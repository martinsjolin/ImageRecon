function [interpolatedProjection interpolationAngles] = partialFanbeamToParallel(estimatedProjection, geometryParameters, angleRangeCoveredDeg)
    %note that this function takes angles in radians as input and converts
    %them into degrees which is required by iradon
    %
    %For the moment, only square images are possible!
    %Interpolate to parallel beam geometry and backproject with iradon.
    %
    %This function expects that the projection data covers eiher 180 or 360,
    %degrees, specified by the optional parameter angleRangeCoveredDeg. (Default:180)
    %
    %It is possible to feed NaN data into this function. I either the angle
    %or estimatedProjection is NaN, then that data point is simply skipped.

    if nargin < 3
        angleRangeCoveredDeg=180;
    end
    if angleRangeCoveredDeg ~= 180 && angleRangeCoveredDeg ~= 360
    error('Incorrect value of angleRangeCoveredDeg')
    end
    DEBUG_MODE = false;
    
    %The triangulation algorithm in triscatteredinterp treats the t and
    %theta variables on equal footing. In order to adjust the triangulation
    %grid, we might choose to rescale one of these dimensions relative to
    %the other. All angles are therefore multiplied with angleExtensionFactor
    %before they are entered into the interpolation. Regardless of what
    %value this factor has, it is still true that the input to the function
    %is in radians and the output is in degrees.
    angleExtensionFactor = 100;
    
    nInterpolationPointsSubdivisions = 20;
    %If the computer is unable to interpolate on all points at the same
    %time, the interpolation step can be subdivided with this parameter.
    
    nDetectorElements =length(geometryParameters.detectorElementPositions);
    nProjectionPositions = length(geometryParameters.projectionPositions);
    sizeOfAngles = size(geometryParameters.angles);
    nAngles = sizeOfAngles(2);
    
    if sizeOfAngles(1) > 1
        projectionDependentAngles = true;
    else
        projectionDependentAngles = false;
    end
    
    %Grid with true sample point coordinates in radon space.
    samplePoints = zeros(nDetectorElements*nProjectionPositions*nAngles,2);
    
    sampledValues = reshape(estimatedProjection', nDetectorElements*nProjectionPositions*nAngles, 1);
    
    for projectionPositionNo = 1:nProjectionPositions
        projectionPositionPx = geometryParameters.projectionPositions(projectionPositionNo)/geometryParameters.mmPerPixel;
        for detectorElementPositionNo=1:nDetectorElements
            
%           %NOT NECESSARY, THEREFORE COMMENT detectorElementPositionPx = geometryParameters.detectorElementPositions(detectorElementPositionNo)/geometryParameters.mmPerPixel;
            
            %The projected line is assumed to go from the _center_ of the
            %focal spot to the detector pixel
            
            %Angle of the projected line measured from the focal
            %spot -> detector center line
            
            alpha = atan(geometryParameters.detectorElementPositions(detectorElementPositionNo)/...
                (geometryParameters.sourceToIsocenterDistance+geometryParameters.isocenterToDetectorDistance));
%             
%             alpha = atan(detectorElementPositionPx/...
%                 (geometryParameters.sourceToIsocenterDistance+geometryParameters.isocenterToDetectorDistance)*geometryParameters.mmPerPixel)

            %Unit normal of projected line in a coordinate system where the focal spot lies on the positive y axis
            unitNormal = [cos(alpha), sin(alpha)];
            %Calculate orthogonal distance to the line as a scalar product.
            t = unitNormal*[projectionPositionPx; geometryParameters.sourceToIsocenterDistance/geometryParameters.mmPerPixel];
            lineIndexStart = ((projectionPositionNo-1)*nDetectorElements+(detectorElementPositionNo-1))*nAngles+1;
            samplePoints(lineIndexStart:lineIndexStart+nAngles-1, 1) = t;
            if (projectionDependentAngles == false)
                samplePoints(lineIndexStart:lineIndexStart+nAngles-1, 2) = geometryParameters.angles+alpha;
            else
                samplePoints(lineIndexStart:lineIndexStart+nAngles-1, 2) = geometryParameters.angles(projectionPositionNo,:)+alpha;
            end
        end
    end
    
    %Interpolate at points uniformly spaced in t and theta.
    %Important: this has to be done before we add any wrapped points.
    tMeasuredMin = min(samplePoints(:,1));
    tMeasuredMax = max(samplePoints(:,1));
    %Note: If we have measured data assymetrically in t (i.e. tMeasuredMax != tMeasuredMin),
    %we pick out the largest symmetric set in t around t=0 to interpolate
    %in -  because we need a symmetric set around 0 in order to
    %reconstruct the picture
    tMax = min(abs(tMeasuredMax), abs(tMeasuredMin));
    tMin = -tMax;

    %The number of t values should be chosen to that the step in t is the
    %same as for the original samplings.
    
    if(nProjectionPositions*nDetectorElements == 1)
        tStepPreliminary=1; %Avoids divide by zero
    else
        tStepPreliminary = (tMeasuredMax-tMeasuredMin)/(nProjectionPositions*nDetectorElements-1);
        %NOT tMax-tMin. If we crop the sinogram, we get fewer
        %interpolation points but with the same spacing as before.
    end
    
    %Find a step size approximately the same size as tStepPreliminary
    %but chosen so that an integer number of steps fit in between
    %tMeasuredMin and tMeasuredMax.
    nTStepsInInterpolation = round((tMax-tMin)/tStepPreliminary);
    tStep =(tMax-tMin)/nTStepsInInterpolation;
    
    %These are in radians
    angleMin = 0;
    angleStep = (angleRangeCoveredDeg*pi/180)/nAngles;
    angleMax = (angleRangeCoveredDeg*pi/180)-angleStep;
    interpolationAngles = angleMin:angleStep:angleMax;
    
    
    %Add some "wrapped" sample points in the edges (for the minimal and
    %maximal theta angles) so that no interpolation point will fall outside
    %the convex hull of the sample points. Note: after this step,
    %samplePoints and sampledValues are not ordered. (These two vectors have the same orders though.)
    
    %If these extra points are not enough to cover the convex hull, it
    %should be possible to wrap more points quite easily (maybe all points).
    
    nPoints = length(sampledValues);
    
    %If we want to wrap only the first projection line of each rotation
    %wrappedMinThetaPoints = samplePoints(1:nAngles:nPoints,:); 
    wrappedMinThetaPoints = samplePoints; %Wrap everything...

    if (angleRangeCoveredDeg==180)
        wrappedMinThetaPoints(:,1) = -wrappedMinThetaPoints(:,1);
    end
    wrappedMinThetaPoints(:,2) = wrappedMinThetaPoints(:,2)+ angleRangeCoveredDeg*pi/180;
    
    %If we want to wrap only the first projection line of each rotation
    %wrappedMinThetaValues = sampledValues(1:nAngles:nPoints);
    wrappedMinThetaValues = sampledValues;
    
    %If we want to wrap only the first projection line of each rotation
    %wrappedMaxThetaPoints = samplePoints(nAngles:nAngles:nPoints,:);
    wrappedMaxThetaPoints = samplePoints; %Wrap everything...
    if (angleRangeCoveredDeg==180)
        wrappedMaxThetaPoints(:,1) = -wrappedMaxThetaPoints(:,1);
    end
    wrappedMaxThetaPoints(:,2) = wrappedMaxThetaPoints(:,2)-angleRangeCoveredDeg*pi/180;
    %If we want to wrap only the first projection line of each rotation
    %wrappedMaxThetaValues = sampledValues(nAngles:nAngles:nPoints);
    wrappedMaxThetaValues = sampledValues;
    
    samplePoints = [samplePoints; wrappedMaxThetaPoints; wrappedMinThetaPoints];
    sampledValues = [sampledValues; wrappedMaxThetaValues; wrappedMinThetaValues];
    
    %Memory is an issue
    clear wrappedMinThetaPoints wrappedMinThetaValues wrappedMaxThetaPoints wrappedMaxThetaValues
    
    %Take the angle extension factor into account - see a comment at the top
    %of the file.
    samplePoints(:,2)= samplePoints(:,2)*angleExtensionFactor;
    
    %Remove any NaN values that might be present.
    nanIndices = (isnan(samplePoints(:,2)) | isnan(sampledValues));
    samplePoints(nanIndices,:) = [];
    sampledValues(nanIndices) = [];
    
    interpolant = TriScatteredInterp(samplePoints, sampledValues, 'linear');
    %linear interpolation gives better-looking pictures than the other
    %kinds of interpolation (natural neighbor and nearest neighbor)
    %natural neighbor makes matlab crash for the knee dataset!
    %(segmentation fault)
    
    %Debug code - invoke if a subset of the sample and interpolation grids
    %should be plotted later
    %This code will always be skipped except is the DEBUG_MODE vairable is
    %set to truewhile stepping through the code
    if(DEBUG_MODE ==true)
        %This code allows the user to specify the boundaries of the plotted
        %region with pixel positions in the pre-interpolation sinogram.
        %This assumes that t increases as we move down in the sinogram
        %and that theta increases as we move right in the sinogram.
        %These formulas are approximate of course.
        gridPlotProjectionMin=2033;
        gridPlotProjectionMax=2123;
        gridPlotAngleMin = 194;
        gridPlotAngleMax = 211;
        gridPlotTMin = tMeasuredMin+(tMeasuredMax-tMeasuredMin)*(gridPlotProjectionMin-1)/(size(estimatedProjection,1)-1);
        gridPlotTMax = tMeasuredMin+(tMeasuredMax-tMeasuredMin)*(gridPlotProjectionMax-1)/(size(estimatedProjection,1)-1);
        gridPlotThetaUnextendedMin = min(min(geometryParameters.angles))+(max(max(geometryParameters.angles))-min(min(geometryParameters.angles)))*(gridPlotAngleMin-1)/(size(geometryParameters.angles,2)-1);
        gridPlotThetaUnextendedMax = min(min(geometryParameters.angles))+(max(max(geometryParameters.angles))-min(min(geometryParameters.angles)))*(gridPlotAngleMax-1)/(size(geometryParameters.angles,2)-1);
        gridPlotThetaExtendedMin = gridPlotThetaUnextendedMin*angleExtensionFactor;
        gridPlotThetaExtendedMax = gridPlotThetaUnextendedMax*angleExtensionFactor;

        gridPlotTNumberOfSteps = 100;
        gridPlotThetaNumberOfSteps = 100;

        %Hard to vectorize without wasting memory(?) 
        samplePointsToPlot = [];
        sampledValuesToPlot = [];
        for sampleIndex = 1:size(samplePoints,1)
            if(samplePoints(sampleIndex,1)>gridPlotTMin && samplePoints(sampleIndex,1)<gridPlotTMax...
                && samplePoints(sampleIndex,2) > gridPlotThetaExtendedMin && samplePoints(sampleIndex,2) < gridPlotThetaExtendedMax)
                samplePointsToPlot = [samplePointsToPlot; samplePoints(sampleIndex,:)];
                sampledValuesToPlot = [sampledValuesToPlot; sampledValues(sampleIndex)];
            end    
        end
        plot(samplePointsToPlot(:,1),samplePointsToPlot(:,2)/angleExtensionFactor,'o')
        xlabel('T')
        ylabel('theta (rad)')
    end
    
    %We sometimes want to put a breakpoint here
    clear samplePoints sampledValues; %Memory is an issue.
    
    [T,theta] = meshgrid(tMin:tStep:tMax, interpolationAngles);
    interpolationPoints = [reshape(T, numel(T), 1),reshape(theta, numel(theta), 1)];
    clear T theta; %Memory is an issue.
    
    %Take the angle extension factor into account - see a comment at the top
    %of the file.
    interpolationPoints(:,2)= interpolationPoints(:,2)*angleExtensionFactor;
    
    %This code will always be skipped except is the DEBUG_MODE vairable is
    %set to truewhile stepping through the code
    
    if(DEBUG_MODE ==true)
        %Debug code - invoke if a subset of the sample and interpolation grids
        %should be plotted
        %Hard to vectorize without wasting memory(?)
        interpolationPointsToPlot = [];
        interpolationValuesToPlot = [];
        for interpolationIndex = 1:size(interpolationPoints,1);
            if(interpolationPoints(interpolationIndex,1)>gridPlotTMin && interpolationPoints(interpolationIndex,1)<gridPlotTMax...
                    && interpolationPoints(interpolationIndex,2) > gridPlotThetaExtendedMin && interpolationPoints(interpolationIndex,2) < gridPlotThetaExtendedMax)
                interpolationPointsToPlot = [interpolationPointsToPlot; interpolationPoints(interpolationIndex,:)];
            end
        end
        interpolatedValuesToPlot = interpolant(interpolationPointsToPlot);
        
        gridPlotTStep = (gridPlotTMax-gridPlotTMin)/gridPlotTNumberOfSteps;
        gridPlotThetaExtendedStep = (gridPlotThetaExtendedMax-gridPlotThetaExtendedMin)/gridPlotThetaNumberOfSteps;
        tGridPlot = gridPlotTMin:gridPlotTStep:gridPlotTMax;
        thetaExtendedGridPlot = gridPlotThetaExtendedMin:gridPlotThetaExtendedStep:gridPlotThetaExtendedMax;
        [TGrid ThetaGrid] = meshgrid(tGridPlot, thetaExtendedGridPlot);
        gridValues = interpolant(TGrid,ThetaGrid);
        surfc(TGrid,ThetaGrid/angleExtensionFactor, gridValues)
        shading interp
        hold on;
        plot3(samplePointsToPlot(:,1),samplePointsToPlot(:,2)/angleExtensionFactor,sampledValuesToPlot,'o')
        plot3(interpolationPointsToPlot(:,1),interpolationPointsToPlot(:,2)/angleExtensionFactor,interpolatedValuesToPlot,'r+')
        xlabel('T')
        ylabel('theta (rad)')
        zlabel('p')
        %Also make a 2D plot with just the sample points
        figure;
        plot(samplePointsToPlot(:,1),samplePointsToPlot(:,2)/angleExtensionFactor,'o')
        hold on
        xlabel('T')
        ylabel('theta (rad)')
        plot(interpolationPointsToPlot(:,1),interpolationPointsToPlot(:,2)/angleExtensionFactor,'r+')
    end
    
    %We sometimes want to put a breakpoint here
    nInterpolationPoints  = size(interpolationPoints,1);
    %test code
    %nInterpolationPointsSubdivisions=nInterpolationPoints;
    %subdivisionStep=1;
    %End of test code
    subdivisionStep = ceil(nInterpolationPoints/nInterpolationPointsSubdivisions);
    subdivisionIndices = 1:subdivisionStep:nInterpolationPoints;
    if (subdivisionIndices(end) ~= nInterpolationPoints+1)
        subdivisionIndices = [subdivisionIndices, nInterpolationPoints+1];
    end
    
    interpolatedValues =  zeros(nInterpolationPoints,1);
    for subdivisionNo=1:nInterpolationPointsSubdivisions
        subdivisionStartIndex = subdivisionIndices(subdivisionNo);
        subdivisionEndIndex = subdivisionIndices(subdivisionNo+1)-1;
        
        interpolatedValues(subdivisionStartIndex:subdivisionEndIndex,:) = interpolant(interpolationPoints(subdivisionStartIndex:subdivisionEndIndex,:));
    end
    
    %Keep in mind that the number of t values where the image is
    %interpolated is the number of _steps_ in t plus 1.
    interpolatedProjection = reshape(interpolatedValues, length(interpolationAngles), (nTStepsInInterpolation+1));
    interpolatedProjection = interpolatedProjection';
end