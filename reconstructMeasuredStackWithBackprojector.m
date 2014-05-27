%The user must supply inFilename and outFilename

% plaque_name = '777';
% plaque_date = '2013_12_18';

inFilename = ['~/Dropbox/CT_MATHEMATICS/plaque_images/plaque_' plaque_name '_' plaque_date '/plaque_' plaque_name '_measuredDataPhC.mat'];
outFilename = ['~/Dropbox/CT_MATHEMATICS/plaque_images/plaque_' plaque_name '_' plaque_date '/plaque_' plaque_name '_reconstructedPhC.mat'];


if(~exist('outFilename','var'))
    error('outFilename must be supplied');
end

load(inFilename, 'energies', 'scanParameters', 'geometryParameters', 'registeredPValues', 'nBins', 'nSlices', 'measured360Degrees')

%geometryParameters.angles=geometryParameters.angles+17;
%Sometimes it is necessary to add a constant angle here to make geometryParameters.angles lie close to 0.

nDetectorElements =length(geometryParameters.detectorElementPositions);
nProjectionPositions = length(geometryParameters.projectionPositions);
nAngles=size(geometryParameters.angles,2);

%registeredPvalues(registeredPValues >10) = 10;
outputSize=[5000,5000];

makeScaledImages=true; %If the reconstructed stack is inconveniently large we can make downsized images too.
scaleFactor = 0.1; %For making scaled images
scaledSize = round(scaleFactor*outputSize);

sliceStart=min(slice_nr);
sliceEnd=max(slice_nr);

%backprojector = @backprojectPartialFanbeam;

%old code:
% %Regions of interest for calculating CNR and quantifying beam hardening artifacts.
% ROISizePixels = round(50*outputSize(1)/450); 
%     
% %The ROI position coordinates are the center coordinates of the square
% %ROI regions with sides given by ROISizePixels.
% %These two should have the same distance to the isocenter.
% ROIBackgroundPos =round([331/450,233/450]*outputSize(1)); %Hole in the plastic 
% ROITargetPos = round([235/450,333/450]*outputSize(1)); %Plastic
%     
% %These two are used for cupping artifact quantization
% ROICenterPos = round([225/450,225/450]*outputSize(1));
% ROIEdgePos = round([68/450,225/450]*outputSize(1));
% 
% ROIParameters = struct('ROICenterPos', ROICenterPos, 'ROICenterSizePixels', ROISizePixels, 'ROIEdgePos', ROIEdgePos, 'ROIEdgeSizePixels', ROISizePixels, 'ROITargetPos', ROITargetPos, 'ROITargetSizePixels', ROISizePixels, 'ROIBackgroundPos', ROIBackgroundPos, 'ROIBackgroundSizePixels', ROISizePixels);

%Make sure that the detector array is correctly placed relative to the
%isocenter. These are the projection positions that are not thrown away.
geometryParameters.projectionPositions = (49*0.4/2.1093)*[-3:3]-0.3;
%geometryParameters.projectionPositions =(-1.671*30:1.671:1.671*30)-27*0.050;
%geometryParameters.projectionPositions =(-1.671*30:1.671:1.671*30)+32*0.050; %27 or 28 would probably be better, at least for image_2011_02_25
%geometryParameters.projectionPositions=(-1.671*30:1.671:1.671*30)+33*0.050;
%geometryParameters.projectionPositions=(-1.671*29.5:1.671:1.671*29.5)+33*0.050;
%geometryParameters.projectionPositions =(-1.671*30:1.671:1.671*30)+15.5*0.050; %for Baby phantom images 2011-04-01; 2011-04-02
%geometryParameters.projectionPositions =(-1.671*30:1.671:1.671*30)+1.25*0.050; %Good for knee images 2011-04-08, 2011-04-09 and 2011-04-11
%Adding a value lowers the left part of the sinogram relative to the right
%part, when the sinogram is plotted twice with a command such as
%imagesc([ interpolatedProjection(:,1:387,1,1) flipud(interpolatedProjection(:,1:387,1,1))])
%(1:387 can be replaced with : if only 180 degrees have been measured.)
%If the right part lies N pixels above the left part, add N/4 times
%the pixel size to geometryParameters.projectionPositions.

%If the final image should be rotated
rotationAngleDeg=233.8; %positive = counter clockwise

reconstructedImages = zeros(outputSize(1),outputSize(2),nBins,sliceEnd-sliceStart+1);
if (makeScaledImages == true)
    scaledImages = zeros(scaledSize(1), scaledSize(2), nBins, sliceEnd-sliceStart+1);
end

if (measured360Degrees ==true)
    angleRangeCoveredDeg=360;
else
    angleRangeCoveredDeg=180;
end

for sliceNo=sliceStart:sliceEnd
    fprintf('Reconstructing slice %d.\n',sliceNo);
    %A single slice version of the geometry parameters, for feeding into
    %partialFanbeamToParallel, containing only the angles for one slice.
    singleSliceGeometryParameters = geometryParameters;
    singleSliceGeometryParameters.angles = singleSliceGeometryParameters.angles(:,:,sliceNo);
    for binNo=1:nBins
        fprintf('  Reconstructing bin %d.\n',binNo);
        [currentInterpolatedProjection currentInterpolationAngles] = partialFanbeamToParallel(registeredPValues(:,:,binNo,sliceNo), singleSliceGeometryParameters, angleRangeCoveredDeg);
        currentInterpolationAngles=currentInterpolationAngles+rotationAngleDeg*pi/180; %If we want a rotated output image
        currentInterpolatedProjection=interp1(1:1:size(currentInterpolatedProjection,1),currentInterpolatedProjection,1:0.1:size(currentInterpolatedProjection,1),'cubic');
        
        if(measured360Degrees ==true)
            nAnglesInAveragedProjection = ceil(length(currentInterpolationAngles)/2);
            if (mod(length(currentInterpolationAngles),2)==0)
                %This exponential average simulates a larger aquisition and
                %should cause less metal problems than an arithmetic average...
                %- I hope.
                % MATS : flipud?
                currentAveragedProjection =-log(0.5*(exp(-currentInterpolatedProjection(:,1:nAnglesInAveragedProjection))+...
                    +flipud(exp(-currentInterpolatedProjection(:,nAnglesInAveragedProjection+1:end)))));
            else
                currentAveragedProjection =[-log(0.5*(exp(-currentInterpolatedProjection(:,1:nAnglesInAveragedProjection-1))+...
                    +flipud(exp(-currentInterpolatedProjection(:,nAnglesInAveragedProjection+1:end))))) currentInterpolatedProjection(:,nAnglesInAveragedProjection)];
            end
            reconstructedImages(:,:,binNo,sliceNo) = backprojectRadon(currentAveragedProjection, currentInterpolationAngles(1:nAnglesInAveragedProjection), outputSize);
        else
            reconstructedImages(:,:,binNo,sliceNo) = backprojectRadon(currentInterpolatedProjection, currentInterpolationAngles, outputSize);
        end
        if((sliceNo == sliceStart) && (binNo==1))
            %We need to know the dimensions of the output from
            %partialfanbeamToParallel before we create interpolatedProjection and interpolationAngles
            interpolatedProjection= zeros([size(currentInterpolatedProjection),nBins,sliceEnd-sliceStart+1]);
            %We allow for different bins having different angles, which is not really
            %needed.
            interpolationAngles= zeros(length(currentInterpolationAngles),nBins,sliceEnd-sliceStart+1);
        end
        interpolatedProjection(:,:,binNo,sliceNo) = currentInterpolatedProjection;
        interpolationAngles(:,binNo,sliceNo) = currentInterpolationAngles;
        if makeScaledImages==true
            scaledImages(:,:,binNo, sliceNo) = imresize(reconstructedImages(:,:,binNo,sliceNo),scaleFactor);
        end
    end
end

%plotWithOptimalWeighting(reconstructedImages, ROIParameters)
%figure(1)
%colormap gray;
%axis off;
%axis image
if (makeScaledImages==true)
    save(outFilename, 'outputSize', 'geometryParameters', 'rotationAngleDeg', 'reconstructedImages', 'interpolatedProjection', 'interpolationAngles', 'scaledImages','inFilename', '-v7.3')
else
    save(outFilename, 'outputSize', 'geometryParameters,', 'rotationAngleDeg', 'reconstructedImages', 'interpolatedProjection', 'interpolationAngles','inFilename', '-v7.3')
end
