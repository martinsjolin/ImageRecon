%Constants that must be checked each time:

loadMeasurementRawDataFromFile = false; %Mostly for faster debugging
loadCalibrationParametersFromFile = false;
generateCalibrationMeasurementsFiles = true;
%For debugging. If generateCalibrationMeasurementsFile is true
%and loadRawDataFromFile is false, a file calibration_raw_data is generated
%which contains the number of registered counts in all calibration
%measurements.
purePhotonCounting = true;
nFaultyFramesBeforeData = 1;% TASK: confirm that the 

%If true, the number of counts in the different bins are summed before
%anything else is done, both in the image and in the calibration.
%If this is true, the number of bins is automatically set to 1.
%Also, the badbins list is changed.

calculateCalibrationParametersFromAverage=true;
%Only applies to the logarithmic curve fit calibration
%If this is true, calibration measurements are averaged before curve
%fitting. In this way, problems with N being 0 are avoided.
%If false, each separate measurement are fed into the curve fitting.

measured360Degrees = true; %false means 180 degrees are measured
if(measured360Degrees == true)
    coveredAngleRad = 2*pi;
else
    coveredAngleRad = pi;
end

calculateCalibrationErrors = false; %Only makes sense for the logarithmic curve fitting method.
%Calculates error bars on the calibration curves. This probably slows down
%the program quite much.
 
%This code was replaced by selective aggregation 2011-04-11
% %If we want to reduce the statistics we can select one of the angles out of
% %angleAggregationFactor consecutive ones instead of taking the mean
% selectOneAngleInsteadOfAggregating = false;


% plaque_name = '774';
% plaque_date = '2013_12_18';
% 
sliceIndexStart= min(slice_nr)-1; %Normally starts counting from 0 since this is the numbering of the raw data files
sliceIndexEnd = max(slice_nr)-1;

nProjectionPositions = 7;

nSlices=sliceIndexEnd-sliceIndexStart+1;
baseFolder = ['~/Dropbox/CT_MATHEMATICS/plaque_images/plaque_' plaque_name '_' plaque_date '/']; %Should end with a slash.
countsFilenameString = [baseFolder 'raw_data/TomographicData/SliceNR_%d_TransPos_%d_DAC00Repeat00DAC0ScanSelected_asic%d.asic'];
angleFilenameString = [baseFolder 'raw_data/AngleDecoder/SliceNR_%d_TransPos_%d_AngleDecoder.ad'];

outFilenameString= [baseFolder 'plaque_' plaque_name '_'];
%This string should contain a path and the beginning of the filename, but
%not the extension .mat. Several files starting with the above string will
%be created. It normally ends with underscore.
intermediateFilenameString=outFilenameString;
%This is normally the same as outFilenameString
if purePhotonCounting == true
    intermediateFilenameRawDataString = strcat(intermediateFilenameString, 'RawDataPhC_%d.mat');
else
    intermediateFilenameRawDataString = strcat(intermediateFilenameString, 'RawData_%d.mat');
end
%The number in intermediateFilenameRawDataString should be the slice index.
if(purePhotonCounting==true)
    intermediateFilenameCalibrationParametersString =strcat(intermediateFilenameString,'CalibrationParametersPhC.mat');
else
    intermediateFilenameCalibrationParametersString =strcat(intermediateFilenameString,'CalibrationParameters.mat');
end
if(purePhotonCounting==true)
    outFilenameMeasuredDataString =strcat(outFilenameString,'measuredDataPhC.mat');
else
    outFilenameMeasuredDataString =strcat(outFilenameString,'measuredData.mat');
end
if(purePhotonCounting==true)
    intermediateFilenameCalibrationRawDataString = strcat(intermediateFilenameString,'calibrationRawDataPhC_%d.mat');
else
    intermediateFilenameCalibrationRawDataString = strcat(intermediateFilenameString,'calibrationRawData_%d.mat');
end
skippedLinesInBeginningOfAnglesFile = 2; %These are header lines - not angle measurements
%exposureTimeSeconds= 1.3e-3; %Old code. This is assumed to be the same for calibration and measurements.


badBins= [];
%preBadBins =[10:64:3904, 49:64:3904];
%preBadBins =[2:64:3904, 10:64:3904, 17:64:3904, 20:64:3904, 43:64:3904, 49:64:3904, 55:64:3904];
%Used to cunstruct badBins
%badBins =[preBadBins;ones(1,length(preBadBins))]';
%badBins =[[12:64:512];ones(1,8)*6]';
%badBins =[[29:64:512, 12:64:512];ones(1,16)*6]';
if purePhotonCounting == true
    sizeOfBadBins=size(badBins);
    badBins(:,2) = ones(sizeOfBadBins(1),1);
    clear sizeOfBadBins;
end
%Important: adjust this manually each time a picture is taken. Each row
%contains the detector, bin coordinates of measurements that should be
%replaced by interpolation.
%If one bad pixel gives rise to several bad lines in the sinogram because
%of the detector translation, these lines require separate entries in this
%table. Also, the detector coordinates apply in the coordinates of the
%_cropped_ sinogram, not the original one.
%All angle measurements for that detector and bin number are discarded.
%Note the current inplementation cannot handle more than one pixel wide
%incorrect stripes.

firstUsedProjectionPosition = 1;
lastUsedProjectionPosition = nProjectionPositions;
nReducedProjectionPositions = lastUsedProjectionPosition-firstUsedProjectionPosition+1;

angleAggregationFactor = 1; %This number of measurements are aggregated in the angular direction
%The option anglesUsedPerAggregatedBlock allows one to use a reduced set of data, giving an
%image corresponding to a lower dose than that corresponding to the full
%dataset. Setting anglesUsedPerAggregatedBlock equal to
%numberOfAnglesToBeUsed.
anglesUsedPerAggregatedBlock = angleAggregationFactor; %Keep all angle measurements
%If not all angles are kept, keep a contiguous block
includedElementsInAggregation = zeros(1,angleAggregationFactor);
includedElementsInAggregation(1:anglesUsedPerAggregatedBlock)=1;
fractionAnglesKept = anglesUsedPerAggregatedBlock/angleAggregationFactor;
%fractionAnglesKept is the relevant fraction to multiply with when
%calculating the dose. Note that it does not take into account the angles
%discarded because the scan covered more than the required angular range
%(180 or 360 degrees). It only takes account of the angles thrown away
%by the setting of anglesUsedPerAggregatedBlock.

%Calibration contants
calibrationIndexStart = sliceIndexStart;
%calibrationIndexEnd = 0; %For flatfield
calibrationIndexEnd = sliceIndexEnd+1; %For wedge calibrator
nCalibrations = calibrationIndexEnd-calibrationIndexStart+1; %One may have several calibration sequences during one scan in order to correct for temperature drifts.

useMeasuredCalibrationThicknessAsReference = false;
%If true, measured values of the calibrator thickness are used as reference
%values for the calibration curves. If false, The average measured count
%number over all detector elements is used as reference.
if (useMeasuredCalibrationThicknessAsReference == false)
    nCalibrationThicknesses=11; %Including air scan. For platic wedge
    %nCalibrationThicknesses=1; %For flatfield
else
    %This part is only needed when we use measured calibrator thickneses as
    %calibration reference.
    %The _first_ of the calibration thicknesses must be 0 (= unattenuated
    %beam).
    calibrationWedgeAngleRad = 26.6154 *pi/180;%The angle between those two sides of the triangle _which the x-ray beam passes through_!
    calibrationVerticalDistancesFromCornermm = 4.64*(1:10);
    %The values in calibrationVerticalDistancesFromCorner are the distances
    %from the calibration corner (the angle angle between those two sides of
    %the triangle _which the x-ray beam passes through_) along the vertical
    %side of the wedge to the calibration measurement positions. %Note that the
    %length of calibrationVerticalDistancesFromCorner is one less than the
    %number of calibration thicknesses, since the air scan is not included
    %here.
    calibrationThicknessesmm = [0 tan(calibrationWedgeAngleRad)*calibrationVerticalDistancesFromCornermm];
    %Note: It is not so important to know the absolute length that the beam has
    %passed through, it is the quotients between the different lengths that is
    %important.
    nCalibrationThicknesses = length(calibrationThicknessesmm); %Including air scan
    calibrationMaterial = 'pmma'; %No it's POM actually but we'll use pmma values until I have attenuation values for POM
    muCalibrator = linAtt(energies,calibrationMaterial);
end
calibrationFilenameString =[baseFolder 'raw_data/CalibrationData/CalibrationNR_%d_Step_%d_DAC00Repeat00DAC0ScanSelected_asic%d.asic'];
nCalibrationMeasurements=1200; %Number of measurements for each thickness of the calibration phantom

% %----------------parameters for simple calibration---------------------------
% 
% calibrationFileNameString = '2011_02_28_0_%d_%d_calib.asic';
% 
% countsAddedInPhotonStarvedBins = 0.5; %in order to avoid Nan in the image, we can add a low number to each photon starved bin.
% nParameters = 1;% Just the average count rate
% nCalibrationThicknesses = 1;
% nCalibrationMeasurements = 2000;
% nCalibrations=nProjectionPositions;
% calibrationMaterial = 'pmma'; %No it's POM actually but we'll use pmma values until I have attenuation values for POM
% %Used to create a dummy variable below
% %--------------End of parameters for ugly calibration----------------------

%in order to avoid Nan in the image, we can add a low number to each photon
%starved bin.
%Since this is done after aggregation (in which the mean of the aggregated
%measurements is calculated), the minimum nonzero value that may be
%registered is 1/angleAggregationFactor.
countsAddedInPhotonStarvedBins = 0.5*1/angleAggregationFactor;

%Energies in keV, used to calculate tube spectrum, detection efficiency and
%also the p values in calibration
energyMin = 1;
energyMax = 140; %This should be larger than the kVp of the scan
energyStep = 1; 

energies = energyMin:energyStep:energyMax; 

%Geometry parameters:
focalSpotPoints = 0;
magnificationFactor = 2.128; %Rough value, too many digits included here
detectorElementPositions = -9.8:0.4:9.8;
projectionPositions = (49*0.4/magnificationFactor)*[-3:3]; 
%Distances from the line connecting the center of the detector and the
%source to the origo on the object stage.

%Note: this is before sinogram cropping
%projectionPositions = -1.671*29.5:1.671:1.671*29.5; %Note: this is before sinogram cropping
%NOTE as it is now, projectionPositions is redefined in
%reconstructMeasuredKneeWithBackprojectorAndPlot. Make sure to keep
%the two definitions compatible!


totalDistanceSourceToDetector = 950; %mm, From a measurement of SID=45 cm (with measuring tape) and the magnification factor.
%The above two are only used to calculate the following two parameters:
sourceToIsocenterDistance = totalDistanceSourceToDetector/magnificationFactor;
isocenterToDetectorDistance = totalDistanceSourceToDetector*(1-1/magnificationFactor);
%scaleFactor = 9/188; %Not sure if this matters... %Don't think this is necessary
mmPerPixel = 0.4/magnificationFactor; %Seems to be unused

%These parameters don't actually do anything
%Scan parameters: (more scan parameters are calculated in the end of the script)
kVp = 80; %Note that this affects the calculation of p value from calibrator thickness.
%tubeSpectrum = computeStandardTubeSpectrum(energies, kVp, 1); %The normalization is unimportant here
%detectionEfficiency = computeStandardDetectionEfficiency(energies);
tubeSpectrum = NaN*ones(size(energies)); %WORKAROUND, since tubeSpectrum is not used anyway
detectionEfficiency = NaN*ones(size(energies)); %WORKAROUND, since detectionEfficiency is not used anyway

%nDetectors = length(detectorElementPositions)*length(projectionPositions);
%This comes later, after projectionPositions has been modified.
%nAngles comes later
detectorBlockSize = 0; %No block division of the detector array
%noBowtieFilter = ones(nDetectors, length(energies));
%This also comes later.

%calibrationMethod = 'linearTransform';
calibrationMethod = 'logarithmicCurveFit';
if (strcmp(calibrationMethod,'logarithmicCurveFit')==true)
    nParameters = 3; %Number of calibration parameters
end
%This is used for the calculation of p values in the calibration curve
%fitting.

%Constants that are normally not changed
nDetectorElements = length(detectorElementPositions);
nAngles = 8186; %Move to a place with frequently changed settings
nAsics=5;
nChannelsPerAsic=160;
nOriginalBins=8; %The number of bins in the original data files. Should never be anything else than 8.
if purePhotonCounting == true
    nBins = 1;
else
    nBins = nOriginalBins;
end
%nExtraAnglesInTheEnd is the number of extra angle values, after the ones belonging to the measurement, in the rotation encoder output file
nExtraAnglesInTheEnd = 4; %When using 5 asics, otherwise 20.
%----------------------End of settings part--------------------------------

packagedDate= datestr(now);

diode_layout = load('via_channels.mat');

if(loadMeasurementRawDataFromFile == false)
    for sliceNo = 1:nSlices
        disp(sprintf('Loading slice %d', sliceNo));
        sliceIndex=sliceIndexStart+sliceNo-1;
        angles = zeros(nProjectionPositions, nAngles);
        registered = zeros(nProjectionPositions*nDetectorElements, nAngles, nBins);

        %Let the positive x axis be directed from the phantom position towards the
        %person sitting at the steering computer. Then the pixels of the detector
        %are arranged in the wrong order (decreasing x) in the data files. The
        %translation stage with the knee moves in the negative x direction, so
        %in a coordinate system where the position of the knee is fixed, the
        %fan locations get higher x values for higher projection numbers.


        for i = 1:nProjectionPositions
            %Note that the projection positions are tacken in the right order -
            %it was wrong to take them backwards as I did before.
            projectionPositionIndex = i-1; %This is always equal to
            %the number in the corresponding filename. Note: it assumes values
            %between 0 and nProjectionPositions-1.
            angleFilename=sprintf(angleFilenameString, sliceIndex, projectionPositionIndex);

            disp(sprintf('  Loading projection position %d', projectionPositionIndex +1));

            %The command importdata cannot be used for the counts since there is
            %extra information after the count values.
            %d = importdata(filename, '\t', 0);
            %M = d.data;

            countsInAsicChannels = zeros(nChannelsPerAsic*nAsics, nAngles, nOriginalBins);
            for asicNo=1:nAsics
                countsFilename=sprintf(countsFilenameString, sliceIndex, projectionPositionIndex,asicNo-1);
                countsInAsicChannels((asicNo-1)*nChannelsPerAsic+1:asicNo*nChannelsPerAsic,:,:) = loadFramesFromFile(countsFilename, nAngles);
            end
            %TODO mark measurements with OF as bad
            countsInDelsSummedOverDepthSegments = zeros(nDetectorElements,nAngles,nOriginalBins);
            %For now, just sum over depth segments.
            for delNo = 1:nDetectorElements %Could maybe be vectorized, although it's not very easy.
                countsInDelsSummedOverDepthSegments(delNo,:,:)=sum(countsInAsicChannels(diode_layout.diode_channel(:,delNo),:,:),1);
            end

            d_ang = importdata(angleFilename, '\t', skippedLinesInBeginningOfAnglesFile);

            angles(i,:)=(d_ang.data(1:nAngles,3))'; 
            
            detectorStart = (i-1)*nDetectorElements+1;
            detectorEnd = (i-1)*nDetectorElements+nDetectorElements;

            %The orientation of the elements within one detector is the
            %right one and does not need flipping like before.
            %TODO check if this is right
            
            if purePhotonCounting == false
                %registered(detectorStart:detectorEnd,:,:) = flipdim(countsInDelsSummedOverDepthSegments,1); %Old, for mammography diode
                registered(detectorStart:detectorEnd,:,:) = countsInDelsSummedOverDepthSegments;
            else
                %registered(detectorStart:detectorEnd,:,:) = flipdim(sum(countsInDelsSummedOverDepthSegments,3),1); %Old, for mammography diode
                registered(detectorStart:detectorEnd,:,:) = sum(countsInDelsSummedOverDepthSegments,3);
            end

            %Since the rotation stage is going back and forth, we need to flip the
            %order of both the angles and the counts for every other translation
            %position.

            if mod(projectionPositionIndex,2) == 1 %This is projection number 2,4,6...
                registered(detectorStart:detectorEnd,:,:) = flipdim(registered(detectorStart:detectorEnd,:,:),2);
                angles(i,:) = fliplr(angles(i,:));
            end

        end

        %The unit measured by the rotation sensor is 0.1 degrees.
        angles = angles*pi/(180*10);

        save(sprintf(intermediateFilenameRawDataString,sliceIndex), 'registered', 'angles')
        clear registered angles;
    end
end

if(loadCalibrationParametersFromFile == true)
    load(intermediateFilenameCalibrationParametersString, 'calibrationParameters')
else
%Load the calibration data
    if strcmp(calibrationMethod,'flatfield')==true
        calibrationParameters = zeros(nCalibrations, nDetectorElements, nBins); %Will contain N0, the number of counts in an air scan, for each del and bin.
    elseif strcmp(calibrationMethod,'logarithmicCurveFit')==true
        calibrationParameters = zeros(nCalibrations, nDetectorElements, nBins, nParameters);
        if calculateCalibrationErrors==true
            calibrationHalfWidths = zeros(nCalibrations, nDetectorElements, nBins, nCalibrationThicknesses);
        end
    elseif(strcmp(calibrationMethod,'affineTransform')==true)
        calibrationParameters=struct('AMatrices', zeros(nCalibrations,nDetectorElements,nBins,nBins),...
            'bVectors', zeros(nCalibrations, nDetectorElements, nBins), 'unattenuatedCounts', zeros(nCalibrations,nBins));
        %UnattenuatedCounts is used for log normalization, i.e. it is N0 in
        %p = -log(N/N0). It is the same for all detector elements.
    elseif(strcmp(calibrationMethod,'linearTransform')==true)
        calibrationParameters=struct('AMatrices', zeros(nCalibrations,nDetectorElements,nBins,nBins),...
            'unattenuatedCounts', zeros(nCalibrations,nBins));
        %UnattenuatedCounts is used for log normalization, i.e. it is N0 in
        %p = -log(N/N0). It is the same for all detector elements.
    else
        error(sprintf('Unknown calibration method: %s'),calibrationMethod);
    end
    %These shouldn't be needed anymore
    %if purePhotonCounting == false
    %    normFactor = zeros(nCalibrations, nDetectorElements, nBins, nCalibrationThicknesses);
    %    countRate = zeros(nCalibrations, nDetectorElements, nBins, nCalibrationThicknesses);
    %else
    %    normFactor = zeros(nCalibrations, nDetectorElements,   1, nCalibrationThicknesses);
    %    countRate = zeros(nCalibrations, nDetectorElements, 1, nCalibrationThicknesses);        
    %end

    for i=1:nCalibrations
        calibrationRegistered = zeros(nDetectorElements, nCalibrationMeasurements, nBins, nCalibrationThicknesses);
        calibrationRegisteredFaultyData = zeros(nDetectorElements, nCalibrationMeasurements, nBins, nCalibrationThicknesses);
        
        %Read the calibration files  starting with the last calibration
        %measurement.
        %After this is done, calibrationParameters(1,:,:,:) will contain
        %calibration data corresponding to the top part of the sinogram and
        %calibrationParameters(1,:,:,:) will contain data corresponding to
        %the bottom part.
        
        calibrationIndex = i+calibrationIndexStart-1; %This is always equal to
        %the number in the corresponding filename.
        
        fprintf('Loading data for calibration sequence %d.\n', i);   
        for j = 1:nCalibrationThicknesses
            %+----------------------Important note------------------------+
            %|                                                            |
            %|          The calibration thickness increases (i.e.         |
            %|         the count rate decreases) when j increases.        |
            %|         In the very first images we took it was the        |
            %|                       other way around!                    |
            %+-------------------End of important note--------------------+
            
            %load raw data into matrix M;
            calibrationThicknessIndex = j-1; %This is always equal to
            %the number in the corresponding filename. Note: it assumes values
            %between 0 and nCalibrationThicknesses-1.

            fprintf('  Loading data for calibration thickness %d.\n',j);
            %This code is wasteful with calibratin data in the sense that
            %it discards the data (in a certain frame)
            %from all depth segments if overflow occurs in at least one of them.
            countsInAsicChannels = zeros(nChannelsPerAsic*nAsics, nCalibrationMeasurements, nOriginalBins);
            overflowInAsicChannels = zeros(nChannelsPerAsic*nAsics, nCalibrationMeasurements);
            for asicNo=1:nAsics
                
                
                if strcmp(calibration_number,'all')
                    calibrationFilename=sprintf(calibrationFilenameString,calibrationIndex,calibrationThicknessIndex,asicNo-1);
                else % Condition for using only one calibration
                    disp(['calibration number ' num2str(calibration_number)])
                    calibrationFilename=sprintf(calibrationFilenameString,calibration_number,calibrationThicknessIndex,asicNo-1);
                end

                [currentCounts, currentMeasurementsWithOverflow] = loadFramesFromFile(calibrationFilename, nCalibrationMeasurements);
                countsInAsicChannels((asicNo-1)*nChannelsPerAsic+1:asicNo*nChannelsPerAsic,:,:) = currentCounts;
                overflowInAsicChannels((asicNo-1)*nChannelsPerAsic+1:asicNo*nChannelsPerAsic,:) = currentMeasurementsWithOverflow;
            end
            countsInDelsSummedOverDepthSegments = zeros(nDetectorElements,nCalibrationMeasurements,nOriginalBins);
            overflowInDelsSummedOverDepthSegments = zeros(nDetectorElements,nCalibrationMeasurements); %Will be 1 if overflow has happened in any of the depth segments and 0 otherwise

            %For now, just sum over depth segments.
            for delNo = 1:nDetectorElements %Could maybe be vectorized, although it's not very easy.
                countsInDelsSummedOverDepthSegments(delNo,:,:)=sum(countsInAsicChannels(diode_layout.diode_channel(:,delNo),:,:),1);
                overflowInDelsSummedOverDepthSegments(delNo,:)=(sum(overflowInAsicChannels(diode_layout.diode_channel(:,delNo),:),1)>0); %Will assume values 1 and 0
            end

            calibrationRegisteredFaultyData(:,:,:,j)=repmat(overflowInDelsSummedOverDepthSegments,[1 1 nBins]);
            if purePhotonCounting == false
                %calibrationRegistered(:,:,:,j) = flipdim(countsInDelsSummedOverDepthSegments,1);
                calibrationRegistered(:,:,:,j) = countsInDelsSummedOverDepthSegments; %TEST, added 2013-09-26
            else
                %calibrationRegistered(:,:,:,j) = flipdim(sum(countsInDelsSummedOverDepthSegments,3),1);
                calibrationRegistered(:,:,:,j) = sum(countsInDelsSummedOverDepthSegments,3); %TEST, added 2013-09-26
            end

            %Old code
            %Is this possible to vectorize?
            %for detectorElementNo = 1:nDetectorElements
            %    for binNo = 1:nBins
            %        currentCalibrationRegistered = calibrationRegistered(detectorElementNo,find(~calibrationRegisteredFaultyData(detectorElementNo,:,binNo,j)),binNo,j);
            %        countRate(i,detectorElementNo,binNo,j)=mean(currentCalibrationRegistered)/exposureTimeSeconds;
            %    end
            %end
            
            nBadMeasurements = sum(sum(sum(calibrationRegisteredFaultyData(:,:,:,j))))/nBins;
            
            fprintf('  Discarding %d measurements with overflow (out of %d*%d = %d).\n', nBadMeasurements-nFaultyFramesBeforeData,nDetectorElements,nCalibrationMeasurements-nFaultyFramesBeforeData ,nDetectorElements*(nCalibrationMeasurements-nFaultyFramesBeforeData));
            %-nFaultyFramesBeforeData for the first few measurements which are always bad.
            
        end
        
        %Old code
        %countRateDetectorAverage = sum(countRate,2)/nDetectorElements;
        %for detectorElementNo = 1:nDetectorElements
            %Note: CountRateDetectorAverage still hs 4 dimensions even though
            %he detectorElement dimension is singleton. This is intentional, for
            %squeezing countRateDetectorAverage may remove one dimension too
            %much if there is only one calibration.
            %normFactor(i,detectorElementNo,:,:) = squeeze(countRate(i,detectorElementNo,:,:)./countRateDetectorAverage(i,:,:,:));
        %end
        %cpsUnattenuated(i,:) = countRateDetectorAverage(:,1);
        
        if (generateCalibrationMeasurementsFiles == true)
            calibrationRegisteredLegend = 'detector element, calibration measurement, bin, calibration thickness number';
            save(sprintf(intermediateFilenameCalibrationRawDataString,calibrationIndex),'calibrationRegistered','calibrationRegisteredFaultyData','calibrationRegisteredLegend', 'packagedDate');
        end
        
        %This turned out to be unnecessary
        %unattenuatedMeanRegisteredInCalibration = zeros(nBins,nCalibrations);
        %Loop needed since we must discard faulty data.
        %for binNo = 1:nBins
        %    %Same value in all columns since we only make one unattenuated
        %    %measurement
        %    detectorElementMeanRegistered = zeros(1,nDetectorElements);
        %    for detectorElementNo = 1:nDetectorElements
        %        %DON'T forget that the unattenuated beam is the last
        %        %measurement.
        %        detectorElementMeanRegistered(detectorElement) = mean(calibrationRegistered(detectorElementNo,:,binNo,end));
        %    end
        %    unattenuatedMeanRegisteredInCalibration(binNo,:)=mean(detectorElementMeanRegistered);
        %    %Note: all columns get the same value.
        %end
        
        if (strcmp(calibrationMethod,'flatfield')==true)
            calibrationParameters(i,:,:)=permute(sum(calibrationRegistered.*(calibrationRegisteredFaultyData==0),2)./sum(calibrationRegisteredFaultyData==0,2),[1 3 4 2]);
            %Note: for the flatfield calibration nCalibrationThicknesses
            %should be 1, so that the right hand side of this assignment has only two indices.
        elseif (strcmp(calibrationMethod,'logarithmicCurveFit')==true)
            %Calculate the p values with a monochromatic model
            %(Or would a polychromatic one be better?)
            if(useMeasuredCalibrationThicknessAsReference==false)
                pValuesInCalibration = zeros(nBins,nCalibrationThicknesses);
                pValuesInCalibration(:,1)=0; %Air scans
                for binNo=1:nBins
                    meanForAirScan=meanWithFaultyDataPoints(calibrationRegistered(:,:,binNo,1),...
                            calibrationRegisteredFaultyData(:,:,binNo,1));
                    for calibrationThicknessNo = 2:nCalibrationThicknesses
                        %We skip number 1 since the air scan is already dealt with.
                        currentMean=meanWithFaultyDataPoints(calibrationRegistered(:,:,binNo,calibrationThicknessNo),...
                            calibrationRegisteredFaultyData(:,:,binNo,calibrationThicknessNo));
                        pValuesInCalibration(binNo,calibrationThicknessNo)=-log(currentMean/meanForAirScan);                   
                    end
                end
            else
                %Effective mu value for use in monochromatic model
                averageMu=sum(tubeSpectrum.*muCalibrator)/sum(tubeSpectrum);
                pValuesInCalibration = zeros(nBins,nCalibrationThicknesses);
                for calibrationThicknessNo = 1:nCalibrationThicknesses
                    %Monochromatic model
                    pValuesInCalibration(:,calibrationThicknessNo) = averageMu*calibrationThicknessesmm(calibrationThicknessNo)/10;

                    %For polychromatic model:
                    %pValuesInCalibration(reverseCalibrationThicknessNumber) = ...
                    %     -log(sum(exp(-muCalibrator*calibrationThicknessesmm(reverseCalibrationThicknessNumber))*tubeSpectrum)/sum(tubeSpectrum));
                end
            end
            if (calculateCalibrationParametersFromAverage==true)
                %createCalibrationParametersSectraStyle was originally intended for
                %use with single unbinned measurements, but if some of those
                %measurements are 0, the output value will be NaN. Therefore one can
                %calculate the mean of all the measurements and feed them into
                %createCalibrationParametersSectraStyle. Think about: is there a
                %better way of handling this (= a way which gives better accuracy)?
                calibrationRegisteredMean = zeros(nDetectorElements, 1, nBins, nCalibrationThicknesses);
                for detectorElementNo=1:nDetectorElements
                    for binNo=1:nBins
                        for calibrationThicknessNo=1:nCalibrationThicknesses
                            calibrationRegisteredMean(detectorElementNo,1,binNo,calibrationThicknessNo) =...
                                meanWithFaultyDataPoints(calibrationRegistered(detectorElementNo,:,binNo,calibrationThicknessNo),...
                                calibrationRegisteredFaultyData(detectorElementNo,:,binNo,calibrationThicknessNo));
                        end
                    end
                end
                if calculateCalibrationErrors==true
                    [calibrationParameters(i,:,:,:) calibrationHalfWidths(i,:,:,:)]= createCalibrationParametersSectraStyle(calibrationRegisteredMean, zeros(size(calibrationRegisteredMean)), pValuesInCalibration);
                else
                    calibrationParameters(i,:,:,:) = createCalibrationParametersSectraStyle(calibrationRegisteredMean, zeros(size(calibrationRegisteredMean)), pValuesInCalibration);
                end    
            else
                if calculateCalibrationErrors==true
                    [calibrationParameters(i,:,:,:) calibrationHalfWidths(i,:,:,:)]= createCalibrationParametersSectraStyle(calibrationRegistered, calibrationRegisteredFaultyData, pValuesInCalibration);
                else
                    calibrationParameters(i,:,:,:) = createCalibrationParametersSectraStyle(calibrationRegistered, calibrationRegisteredFaultyData, pValuesInCalibration);
                end
            end
        
        elseif(strcmp(calibrationMethod, 'affineTransform') ==true)
            [AMatricesForThisCalibration bVectorsForThisCalibration unattenuatedCountsForThisCalibration] = createCalibrationParametersAffineTransform(calibrationRegistered, calibrationRegisteredFaultyData);
            calibrationParameters.AMatrices(i,:,:,:)= AMatricesForThisCalibration;
            calibrationParameters.bVectors(i,:,:)= bVectorsForThisCalibration;
            calibrationParameters.unattenuatedCounts(i,:)= unattenuatedCountsForThisCalibration;
        elseif(strcmp(calibrationMethod, 'linearTransform') ==true)
            referenceDetectorNo=1; %Move this somewhere else after testing.
            [AMatricesForThisCalibration unattenuatedCountsForThisCalibration] = createCalibrationParametersLinearTransform(calibrationRegistered, calibrationRegisteredFaultyData, false, squeeze(calibrationRegistered(referenceDetectorNo,:,:,:)), squeeze(calibrationRegisteredFaultyData(referenceDetectorNo,:,:,:)) );
            calibrationParameters.AMatrices(i,:,:,:)= AMatricesForThisCalibration;
            calibrationParameters.unattenuatedCounts(i,:)= unattenuatedCountsForThisCalibration;
        else
            error(sprintf('Unknown calibration method: %s'),calibrationMethod);
        end
        clear calibrationRegistered    
    end

    if calculateCalibrationErrors==true
        save(intermediateFilenameCalibrationParametersString,'calibrationParameters','calibrationHalfWidths');
    else
        save(intermediateFilenameCalibrationParametersString,'calibrationParameters');
    end
end

%Calculate the number of angles necessary to cover the required angle range
%for all slices and all projection positions
numberOfAnglesToBeUsed = 0;
for sliceNo = 1:nSlices
    sliceIndex=sliceIndexStart+sliceNo-1;
    load(sprintf(intermediateFilenameRawDataString,sliceIndex), 'angles');
    for i = 1:nProjectionPositions
        endAngle = angles(i,1)+coveredAngleRad;
        numberOfAnglesBetweenStartAndEnd = find(angles(i,:) < endAngle,1, 'last');
        if numberOfAnglesBetweenStartAndEnd>numberOfAnglesToBeUsed
            numberOfAnglesToBeUsed=numberOfAnglesBetweenStartAndEnd;
        end
    end
end
nAggregatedAngles=round(numberOfAnglesToBeUsed/angleAggregationFactor);
numberOfAnglesToBeUsed = nAggregatedAngles*angleAggregationFactor;
%numberOfAnglesToBeUsed is now divisible by angleAggregationFactor

angleStack = zeros(nReducedProjectionPositions,nAggregatedAngles,nSlices);

for sliceNo = 1:nSlices
    %Load raw data from files, either old files or the ones created by the
    %above part of the program
    sliceIndex=sliceIndexStart+sliceNo-1;
    load (sprintf(intermediateFilenameRawDataString, sliceIndex), 'angles', 'registered')

    %Reduce the dataset
    %Note that the first few frames for each projection position are junk.
    reducedRegistered = zeros(nDetectorElements*nProjectionPositions,nAngles-nFaultyFramesBeforeData,nBins);
    for i = 1:nProjectionPositions
        detectorStart = (i-1)*nDetectorElements+1;
        detectorEnd = detectorStart+nDetectorElements-1;
        if mod(i,2) == 0 %This is projection number 2,4,6... where the angle decreases with time
            reducedRegistered(detectorStart:detectorEnd,:,:)=registered(detectorStart:detectorEnd,1:end-nFaultyFramesBeforeData,:);
        else
            reducedRegistered(detectorStart:detectorEnd,:,:)=registered(detectorStart:detectorEnd,nFaultyFramesBeforeData+1:end,:);
        end
    end
    registered=reducedRegistered;
    clear reducedRegistered;
    %There is no offset beween angles and frames - that was a misconception
    %nFramesDataLagsBehindAngles should be 0.
    reducedAngles =zeros(nProjectionPositions,nAngles-nFaultyFramesBeforeData);
    for i = 1:nProjectionPositions
        if mod(i,2) == 0 %This is projection number 2,4,6... where the angle decreases with time
            reducedAngles(i,:)=angles(i,1+nFramesDataLagsBehindAngles:end-nFaultyFramesBeforeData+nFramesDataLagsBehindAngles);
        else
            reducedAngles(i,:)=angles(i,nFaultyFramesBeforeData+1-nFramesDataLagsBehindAngles:end-nFramesDataLagsBehindAngles);
        end
    end
    angles = reducedAngles;
    clear reducedAngles
    
%   Old code: 
%     %For the first slice: calculate how many angles should be used
%     %For subsequent slices the number of angles is taken to be the same as
%     %for the first slice.
%     if (sliceNo ==1)
%         startAngle = max(angles(:,1));
%         endAngle = startAngle+coveredAngleRad;
%         startAngleIndex=zeros(1,nProjectionPositions);
%         endAngleIndex=zeros(1,nProjectionPositions);
%         for i = 1:nProjectionPositions
%             startAngleIndex(i) = find(angles(i,:) >= startAngle,1, 'first');
%             endAngleIndex(i) = find(angles(i,:) < endAngle,1, 'last');
%             %If the required angular range is not covered between startAngle and
%             %angles(i,end), we must wrap some angles later. For now, just extend endAngle
%             %to something which lies outside the angles matrix 
%             if(angles(i,end) < endAngle)
%                 endAngleIndex(i) = startAngleIndex(i)+round((endAngle-startAngle)/...
%                     (angles(i,end)-angles(i,startAngleIndex(i)))*...
%                     (size(angles,2)-startAngleIndex(i)));
%             end
%             nAnglesBetweenStartAndEnd = endAngleIndex-startAngleIndex+1; %BUGGY!!!
%         end
%         nAggregatedAngles = round(mean(nAnglesBetweenStartAndEnd)/angleAggregationFactor);
%         numberOfAnglesToBeUsed = nAggregatedAngles*angleAggregationFactor;
%         %When processing the first slice, also create the matrix which will
%         %store the flatfielded slices.
%         registeredPValues = zeros(nDetectorElements*nReducedProjectionPositions,nAggregatedAngles,nBins,nSlices);
%         angleStack = zeros(nReducedProjectionPositions,nAggregatedAngles,nSlices);
%     else
%         %Just calculate the starting angle indices
%         startAngle = max(angles(:,1));
%         startAngleIndex=zeros(1,nProjectionPositions);
%         for i = 1:nProjectionPositions
%             startAngleIndex(i) = find(angles(i,:) >= startAngle,1, 'first');
%         end
%     end
%     
%     reducedAngles = zeros(nProjectionPositions,numberOfAnglesToBeUsed);
%     for i = 1:nProjectionPositions
%         if(startAngleIndex(i)+numberOfAnglesToBeUsed-1 <= nAngles-nFaultyFramesBeforeData)
%             reducedAngles(i,:) = angles(i,startAngleIndex(i):startAngleIndex(i)+numberOfAnglesToBeUsed-1);
%         else
%             %If there is not a common range of 180 (or 360) degrees covered for all
%             %translation positions, we might need to wrap a few measurements
%             %from < 0 degrees to just below 180 (or 360) degrees.
%             nWrappedAngles = startAngleIndex(i)+numberOfAnglesToBeUsed-1 - (nAngles-nFaultyFramesBeforeData);
%             if (measured360Degrees == true)
%                 reducedAngles(i,:) = [angles(i,startAngleIndex(i):end) coveredAngleRad+(angles(i,(startAngleIndex(i)-nWrappedAngles):(startAngleIndex(i)-1)))];
%             else
%                 %Don't add pi to the wrapped angles if half a rotation is
%                 %sampled. It is better to keep the angle even though
%                 %registeredPValues will look a bit strange.
%                 reducedAngles(i,:) = [angles(i,startAngleIndex(i):end) (angles(i,(startAngleIndex(i)-nWrappedAngles):(startAngleIndex(i)-1)))];
%             end
%         end
%     end
%     angles = reducedAngles;
%     clear reducedAngles;
% 
%     reducedRegistered = zeros(nProjectionPositions*nDetectorElements,numberOfAnglesToBeUsed,nBins);
%     for i = 1:nProjectionPositions
%         detectorStart = (i-1)*nDetectorElements+1;
%         detectorEnd = detectorStart+nDetectorElements-1;
%         if(startAngleIndex(i)+numberOfAnglesToBeUsed-1 <= nAngles-nFaultyFramesBeforeData)    
%             reducedRegistered(detectorStart:detectorEnd,:,:) = registered(detectorStart:detectorEnd,startAngleIndex(i):startAngleIndex(i)+numberOfAnglesToBeUsed-1,:);
%         else
%             %If there is not a common range of 180 (or 360) degrees covered for all
%             %translation positions, we might need to wrap a few measurements
%             %from < 0 degrees to just below 180 (or 360) degrees.
%             nWrappedAngles = startAngleIndex(i)+numberOfAnglesToBeUsed-1 - (nAngles-nFaultyFramesBeforeData);
%             reducedRegistered(detectorStart:detectorEnd,:,:) = [registered(detectorStart:detectorEnd,startAngleIndex(i):end,:) registered(detectorStart:detectorEnd,(startAngleIndex(i)-nWrappedAngles):(startAngleIndex(i)-1),:)];
%         end
%     end
%     registered = reducedRegistered;
%     clear reducedRegistered;

    %Remove surplus angles. Note that some of the remaining angles are set
    %to NaN below. The code lines here are meant to reduce the set of
    %angles to the maximum number needed in any slice and any projection
    %position. We do not bother about ensuring that the first angle is the
    %same for all projection positions.
    registered=registered(:,1:numberOfAnglesToBeUsed,:);
    angles = angles(:,1:numberOfAnglesToBeUsed,:);

    %Aggregate
    %Note:we take the mean of each set of samples of size
    %anglesUsedPerAggregatedBlock instead of summing
    registered=aggregateSelectively(registered,1,angleAggregationFactor,includedElementsInAggregation)/anglesUsedPerAggregatedBlock;
    angles = aggregateSelectively(angles,1,angleAggregationFactor,includedElementsInAggregation)/anglesUsedPerAggregatedBlock;
    
    %referenceRegistered = createFlatfieldFromCalibration(registered, nDetectorElements, nProjectionPositions, normFactor, countRate, effectiveExposureTimeSeconds);
    registered = registered((firstUsedProjectionPosition-1)*nDetectorElements+1:lastUsedProjectionPosition*nDetectorElements,:,:);

    %This is not needed now that we don't calibrate separately for each
    %translation.
    %calibrationParameters=calibrationParameters(firstUsedProjectionPosition:lastUsedProjectionPosition+1,:,:,:);

    %referenceRegistered = referenceRegistered((firstUsedProjectionPosition-1)*64+1:lastUsedProjectionPosition*64,:,:);

    %For some of the projection positions, the covered angular range may
    %still be more than the required range. Therefore we set the surplus
    %angles to NaN in order to disable those points.
    for projectonPositionNo=1:nProjectionPositions
        endAngle = angles(projectonPositionNo,1)+coveredAngleRad;
        lastAngleIndexToBeUsed = find(angles(projectonPositionNo,:) < endAngle,1, 'last');
        if (lastAngleIndexToBeUsed < nAggregatedAngles)
            angles(projectonPositionNo,lastAngleIndexToBeUsed+1:end)=NaN;
        end
    end

    angleStack(:,:,sliceNo) = angles;

    %In order to avoid NaNs when we have some photon starved bins:
    registered(registered == 0) =countsAddedInPhotonStarvedBins;

    % %------------Simple calibration--------------------------------------------
    % registeredPValues = zeros(nDetectorElements*nReducedProjectionPositions,nAggregatedAngles,nBins);
    % for projectionPositionNo = 1:nProjectionPositions
    %     calibrationNo = projectionPositionNo;
    %     for detectorElementNo = 1:nDetectorElements
    %         projectionNo = detectorElementNo+(projectionPositionNo-1)*nDetectorElements;
    %         for binNo = 1:nBins
    %             registeredPValues(projectionNo,:,binNo) = -log(registered(projectionNo,:,binNo)/...
    %                 calibrationParameters(calibrationNo,detectorElementNo,binNo,1));
    %         end
    %     end
    % end
    % 
    % %------------End of simple calibration-------------------------------------

    if(sliceNo == 1)
        %We can not declare registeredPValues outside the for loop since we
        %do not know how large it should be until now.
        registeredPValues=zeros(size(registered));
        %Size: nDetectorElements*nReducedProjectionPositions,nAngles,nBins,nSlices
    end
    
    %With the curve-fitting kind of calibration, no referenceRegistered matrix
    %is created. Instead, a registeredPValues matrix containing p=integral mu dx
    %(which is equal to -log(N/N_0) _only_ in the ideal case) is created.

    if strcmp(calibrationMethod,'flatfield')==true
        for projectionPositionNo = 1:nProjectionPositions
            calibrationNo = projectionPositionNo; %Assume one air scan per projection position
            for detectorElementNo = 1:nDetectorElements
                projectionNo = detectorElementNo+(projectionPositionNo-1)*nDetectorElements;
                for binNo = 1:nBins
                    registeredPValues(projectionNo,:,binNo,sliceNo) = -log(registered(projectionNo,:,binNo)/...
                        calibrationParameters(calibrationNo,detectorElementNo,binNo,1));
                end
            end
        end        
    elseif strcmp(calibrationMethod,'logarithmicCurveFit')==true    
        for reducedProjectionPositionNo=1:nReducedProjectionPositions

            %Interpolate between the calibrations before and after the slice as
            %we step through the field of view.

            %actualProjectionPositionNo = reducedProjectionPositionNo+firstUsedProjectionPosition-1;
            %interpolationParameter = (actualProjectionPositionNo-0.5)/nProjectionPositions;

            %Goes from almost 0 to almost 1 as actualProjectionPosition goes from 1
            %to nProjectionPositions.
            %This code is made for the situation where one calibration is made
            %between every slice acquisition as well as before and after the
            %whole image acquisiton. So te number of calibrations sould be one
            %more than the numer of slices.

            %A simpler variant: take the mean of the two calibrations.
            interpolationParameter = 0.5;

            %Some lines here could be mover out of the for loop over
            %reducedProjectionPositionNo, but keeping it inside simplifies
            %if one wants to interpolate continuously as the object moves
            %through the field of view.
            
            currentCalibrationParameters = (1-interpolationParameter)*calibrationParameters(sliceNo,:,:,:)+...
                interpolationParameter*calibrationParameters(sliceNo+1,:,:,:);        
            %We keep the first singleton dimension as squeezing can be dangerous
            %It will disappear in a, b and c

            detectorStart = (reducedProjectionPositionNo-1)*nDetectorElements+1;
            detectorEnd = detectorStart+nDetectorElements-1;
            %We could use the permutation [1 2 3 4] but then we would not get rid
            %of the singleton dimension. We do that in the following way:
            a = permute(repmat(currentCalibrationParameters(:,:,:,1),[1 1 1 nAggregatedAngles]),[2 4 3 1]);
            b = permute(repmat(currentCalibrationParameters(:,:,:,2),[1 1 1 nAggregatedAngles]),[2 4 3 1]);
            c = permute(repmat(currentCalibrationParameters(:,:,:,3),[1 1 1 nAggregatedAngles]),[2 4 3 1]);
            %The aggregated count numbers are on the same order as the unaggregated
            %ones so we don't have to divide by angleAggregationFactor here.
            partialRegisteredPValues1 = -b./(2*a) - sqrt(b.^2./(4*a.^2)-(c./a)+log(registered(detectorStart:detectorEnd,:,:))./a);
            partialRegisteredPValues2 = -b./(2*a) + sqrt(b.^2./(4*a.^2)-(c./a)+log(registered(detectorStart:detectorEnd,:,:))./a);
            %The correct solution is likely to be the one closest to (log N-c)/b
            %which is what one will get if one sets a = 0
            %By observation: it seems that the solution number one is the right
            %one if and only if a is positive. (Possible speedup by exploiting!)
            approximateSolution = (log(registered(detectorStart:detectorEnd,:,:))-c)./b;
            indicesWhereSolution1IsBest = (abs(partialRegisteredPValues1-approximateSolution)<abs(partialRegisteredPValues2-approximateSolution));
            partialRegisteredPValues = zeros(nDetectorElements,nAggregatedAngles,nBins);
            partialRegisteredPValues(indicesWhereSolution1IsBest) = partialRegisteredPValues1(indicesWhereSolution1IsBest);
            partialRegisteredPValues(~indicesWhereSolution1IsBest) = partialRegisteredPValues2(~indicesWhereSolution1IsBest);
            if(numel(find(a==0))>0)
                partialRegisteredPValues(a==0)=approximateSolution;
            end
            %If a is zero, we can't use the formula for the solution of a
            %quadratic equation. Note that the approximate solution is in fact
            %exact in this special case.
            registeredPValues(detectorStart:detectorEnd,:,:,sliceNo) = real(partialRegisteredPValues);
            %Sometimes the square root formula above gives complex numbers. We
            %take the real part in order to avoid error messages later.
        end
    elseif(strcmp(calibrationMethod,'affineTransform')==true)
        %Take the mean of the calibrations before and after the slice
        interpolationParameter = 0.5;
        
        currentAMatrices = permute((1-interpolationParameter)*calibrationParameters.AMatrices(sliceNo,:,:,:)+...
            interpolationParameter*calibrationParameters.AMatrices(sliceNo+1,:,:,:),[2 3 4 1]);
        currentbVectors = permute((1-interpolationParameter)*calibrationParameters.bVectors(sliceNo,:,:)+...
            interpolationParameter*calibrationParameters.bVectors(sliceNo+1,:,:),[2 3 1]);
        currentUnattenuatedCounts = (1-interpolationParameter)*calibrationParameters.unattenuatedCounts(sliceNo,:)+...
            interpolationParameter*calibrationParameters.unattenuatedCounts(sliceNo+1,:);
        currentCalibrationParameters = struct('AMatrices', currentAMatrices, 'bVectors', currentbVectors, 'unattenuatedCounts', currentUnattenuatedCounts);        
                
        %A possible speedup(?) would be to make the calibration for all lines in the
        %sinogram that correspond to the same detector element at the same
        %time.
        reducedProjectionNumber =1;
        for reducedProjectionPositionNo=1:nReducedProjectionPositions
            for detectorElementNo=1:nDetectorElements
                registeredPValues(reducedProjectionNumber,:,:,sliceNo) = -log(affineTransformFunction(squeeze(registered(reducedProjectionNumber,:,:)),currentCalibrationParameters,detectorElementNo)./(ones(size(registered,2),1)*currentUnattenuatedCounts));
                reducedProjectionNumber = reducedProjectionNumber+1;
            end
        end
    elseif(strcmp(calibrationMethod,'linearTransform')==true)
        %Take the mean of the calibrations before and after the slice
        interpolationParameter = 0.5;
        
        currentAMatrices = permute((1-interpolationParameter)*calibrationParameters.AMatrices(sliceNo,:,:,:)+...
            interpolationParameter*calibrationParameters.AMatrices(sliceNo+1,:,:,:),[2 3 4 1]);
        currentUnattenuatedCounts = (1-interpolationParameter)*calibrationParameters.unattenuatedCounts(sliceNo,:)+...
            interpolationParameter*calibrationParameters.unattenuatedCounts(sliceNo+1,:);
        currentCalibrationParameters = struct('AMatrices', currentAMatrices, 'unattenuatedCounts', currentUnattenuatedCounts);        
                
        %A possible speedup(?) would be to make the calibration for all lines in the
        %sinogram that correspond to the same detector element at the same
        %time.
        reducedProjectionNumber =1;
        for reducedProjectionPositionNo=1:nReducedProjectionPositions
            for detectorElementNo=1:nDetectorElements
                registeredPValues(reducedProjectionNumber,:,:,sliceNo) = -log(linearTransformFunction(squeeze(registered(reducedProjectionNumber,:,:)),currentCalibrationParameters,detectorElementNo)./(ones(size(registered,2),1)*currentUnattenuatedCounts));
                reducedProjectionNumber = reducedProjectionNumber+1;
            end
        end
    end
end

%Interpolate in logNormalizedRegistered for all bad bins

sizeOfBadBins=size(badBins);
nBadBins=sizeOfBadBins(1);

for i=1:nBadBins
    registeredPValues(badBins(i,1),:,badBins(i,2),:) = ...
        (registeredPValues(badBins(i,1)-1,:,badBins(i,2),:)+...
        registeredPValues(badBins(i,1)+1,:,badBins(i,2),:))/2;
end

angles=angleStack;

%Scan parameters again
projectionPositions=projectionPositions(firstUsedProjectionPosition:lastUsedProjectionPosition);
nDetectors = length(detectorElementPositions)*length(projectionPositions);
nAngles = nAggregatedAngles;
%No bowtie filter for now:
noBowtieFilter = ones(nDetectors, length(energies));
%Important note: No bowtie filter means ones in all elements!

%Create the parameters which must be bundled with the data when passed to
%the reconstruction algorithm.

geometryParameters = struct('focalSpotPoints', focalSpotPoints,'angles', angles, 'detectorElementPositions', detectorElementPositions, 'projectionPositions', projectionPositions, 'sourceToIsocenterDistance', sourceToIsocenterDistance, 'isocenterToDetectorDistance', isocenterToDetectorDistance, 'mmPerPixel', mmPerPixel);
scanParameters = struct('detectionEfficiency', detectionEfficiency, 'tubeSpectrum', tubeSpectrum, 'nDetectors', nDetectors, 'nAngles', nAngles, 'detectorBlockSize', detectorBlockSize, 'bowtieFilter', noBowtieFilter);

save(outFilenameMeasuredDataString, 'energies', 'scanParameters', 'geometryParameters', 'registeredPValues', 'calibrationParameters', 'nBins', 'nSlices', 'packagedDate', 'fractionAnglesKept', 'measured360Degrees','-v7.3');
