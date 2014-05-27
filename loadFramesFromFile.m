function [counters, measurementsWithOverflow, measurementsWithTE] = loadFramesFromFile(filename, varargin)
    %[counters, measurements_with_overflow, measurementsWithTE] = loadFramesFromFile(filename, nChannels, nFrames, nEnergyBins)
    %Can also be called as:
    %[counters, measurements_with_overflow, measurementsWithTE] = loadFramesFromFile(filename, nChannels, nFrames)
    %[counters, measurements_with_overflow, measurementsWithTE] = loadFramesFromFile(filename, nFrames)
    %[counters, measurements_with_overflow, measurementsWithTE] = loadFramesFromFile(filename)
    %Input: the filename where the data is stored and the number of
    %channels, frames and energy bins in the file.
    %If nFrames is NaN, the function will sense the file length
    %automatically and return as many whole frames and fit in the file. If
    %the file length is not divisible by the frame length, the
    %overshooting bytes in the end are ignored.
    %The function can also be called with fewer parameters, whereby
    %parameters not passed get their default values: nChannels=160, nFrames=NaN (autodetect), nEnergyBins=8
    %Output:
    %counters, contains the counter values. Size: (nChannels, nFrames, nEnergyBins)
    %measurementsWithOverflow, measurementsWithTE: contains a 1 for those
    %channels and frames where an overflow or a TE error, repectively, has
    %occurred. Both have size (nChannels, nFrames)
    %
    %This function does not currently check for parity errors (which takes time)
    
    if(nargin==4)
        nEnergyBins=varargin{3};
    else
        nEnergyBins=8;
    end
    if(nargin>=3)
        nChannels=varargin{1};
        nFrames=varargin{2};
    elseif nargin==2
        nChannels=160;
        nFrames=varargin{1};
    else
        nChannels=160;
        nFrames=NaN; %Means autodetect nmber of frames.
    end
    
    if(nEnergyBins~=8)
        error('This function does not support anything else than 8 energy bins.')
        %Reduces the risk of bugs. I don't see any reason why this should
        %have to be changed in the future.
    end
    
    shiftedLfsrEncodingTable =[1, 128, 64, 160, 208, 104, 52, 26, 141, 198, 227, 241, 120, 188, 94, 47, 23, 11, 133, 194, 225, 240, 248, 252, 254, 255, 127, 63, 159, 79, 39, 19, 9, 132, 66, 161, 80, 40, 148, 202, 229, 242, 249, 124, 190, 95, 175, 87, 171, 85, 170, 213, 234, 117, 58, 29, 14, 7, 131, 193, 96, 48, 24, 140, 70, 163, 81, 168, 212, 106, 53, 154, 205, 102, 51, 153, 76, 166, 211, 233, 244, 250, 253, 126, 191, 223, 239, 247, 123, 61, 158, 207, 103, 179, 217, 236, 118, 187, 221, 238, 119, 59, 157, 78, 167, 83, 169, 84, 42, 149, 74, 165, 82, 41, 20, 138, 69, 34, 145, 72, 164, 210, 105, 180, 90, 45, 22, 139, 197, 98, 49, 152, 204, 230, 115, 57, 156, 206, 231, 243, 121, 60, 30, 143, 199, 99, 177, 216, 108, 54, 27, 13, 134, 67, 33, 16, 136, 68, 162, 209, 232, 116, 186, 93, 174, 215, 235, 245, 122, 189, 222, 111, 183, 219, 237, 246, 251, 125, 62, 31, 15, 135, 195, 97, 176, 88, 44, 150, 203, 101, 178, 89, 172, 214, 107, 181, 218, 109, 182, 91, 173, 86, 43, 21, 10, 5, 130, 65, 32, 144, 200, 228, 114, 185, 220, 110, 55, 155, 77, 38, 147, 73, 36, 146, 201, 100, 50, 25, 12, 6, 3, 129, 192, 224, 112, 184, 92, 46, 151, 75, 37, 18, 137, 196, 226, 113, 56, 28, 142, 71, 35, 17, 8, 4, 2, 0];
    %The value at position i of this table is the LFSR encoding of i-1.
    %The reason for this shift of one step is that we want to be able to
    %encode 0 as well.
    frameLengthBytes=1445;
    skippedDataPositionsInBeginningOfFrame = 2; %Number of the first byte in each frame that contains a counter value (channel 0, lowest bin)
    
    [fileHandle message] = fopen(filename,'r');
    if length(message)>0 %A successful file opening gives an empty string as error message
        error(message)
    end
    rawData = fread(fileHandle,'uint8');
    fclose(fileHandle);
    nBytesInRawData = length(rawData);
    if(isfinite(nFrames))
        if(nBytesInRawData<nFrames*frameLengthBytes)
            nFramesOld=nFrames;
            nFrames = floor(nBytesInRawData/frameLengthBytes);
            warning('Data file is too short. Expected %d bytes (=%d frames) but the raw data file contains %d bytes.\nReturning %d frames.',nFramesOld*frameLengthBytes,nFramesOld,nBytesInRawData,nFrames)
        elseif(nBytesInRawData>nFrames*frameLengthBytes)
            warning('Data file is too long. Expected %d bytes (=%d frames) but the raw data file contains %d bytes.\nReturning %d frames.',nFrames*frameLengthBytes,nFrames,nBytesInRawData,nFrames)
        end
    else %Autodetect number of frames
        nFrames = floor(nBytesInRawData/frameLengthBytes);
        if(nBytesInRawData~=nFrames*frameLengthBytes)
            warning('File does not contain an even nuber of frames. Returning %d bytes (=%d frames) out of the %d bytes in the raw data file.',nFrames*frameLengthBytes,nFrames,nBytesInRawData)
        end
    end
    
    counters = zeros(nChannels,nFrames, nEnergyBins);
    controlBytes = zeros(nChannels,nFrames);
    for currentFrameNo = 1:nFrames
        currentFrameData=rawData((currentFrameNo-1)*frameLengthBytes+1:currentFrameNo*frameLengthBytes);
        for binNo=1:nEnergyBins
            counters(:,currentFrameNo,binNo) = currentFrameData(...
                skippedDataPositionsInBeginningOfFrame+binNo:nEnergyBins+1:skippedDataPositionsInBeginningOfFrame+binNo+(nEnergyBins+1)*(nChannels-1));
            controlBytes(:,currentFrameNo) = currentFrameData(...
                skippedDataPositionsInBeginningOfFrame+nEnergyBins+1:nEnergyBins+1:skippedDataPositionsInBeginningOfFrame+nEnergyBins+1+(nEnergyBins+1)*(nChannels-1));
        end
    end
    %The control bytes were decoded in the FPGA, although they were
    %never encoded in the first place. So we need to encode them again.
    reencodedControlBytes = shiftedLfsrEncodingTable(controlBytes+1);
    measurementsWithOverflow = bitand(reencodedControlBytes,2)~=0;
    measurementsWithTE = bitand(reencodedControlBytes,1)~=0;% non-temperature encoded value, a zero with a one above in the DAC comparators.
end