
function weighted = weight(spectralData, weightFunction)
    %function weighted = weight(spectralData, weightFunction)
    %SpectralData is a stack of nImages images. Size: (nRows,nColumns,nImages)
    %weightFunction is a row or column vector with nImages entries.
    %Returns a linear combination of the stack of images with the weighting
    %coefficients specified in weightFunction.
    %Can be either image based or projection based weighting, depending on 
    %whether spectralData contains intensities or projections.
    %Don't normalize the weight function since we might want to pass in one
    %that sums to 0.
    if length(weightFunction) == 1
        weighted=spectralData;
    else
        sizeOfImage = size(spectralData);
        weighted =zeros(sizeOfImage(1), sizeOfImage(2));
        for i = 1:sizeOfImage(3)
            weighted = weighted + spectralData(:,:,i)*weightFunction(i); 
        end
    end
end
