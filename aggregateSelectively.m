function reducedmatrix= aggregateSelectively(matrix, rowFactor, columnFactor, includedElements)
    %Reduces the number of rows and columns by aggregating the
    %matrix values in blocks of size rowFactor*columnFactor.
    %If the matrix is three dimensional, the third dimension is left
    %unchanged.
    %This version of the function allows the user to select which elements
    %in each block will be included in the sum and which will be left out.
    %This is useful for example for making a lower dose image out of a
    %sampled dataset.
    %includedElements is a matrix of size rowFactor*columnFactor. It
    %contains ones for the elements of each block which will be included in
    %the sum and zeros otherwise.
    
    sizeOfMatrix = size(matrix);
    if(numel(sizeOfMatrix) >2)
        reducedSize = [sizeOfMatrix(1)/rowFactor, sizeOfMatrix(2)/columnFactor, sizeOfMatrix(3)];
        includedElements = repmat(includedElements,[1 1 sizeOfMatrix(3)]);
    else
        reducedSize = [sizeOfMatrix(1)/rowFactor, sizeOfMatrix(2)/columnFactor, 1];
    end
    reducedmatrix = zeros(reducedSize);
    for i=1:reducedSize(1)
        iOriginalStart = (i-1)*rowFactor+1;
        iOriginalEnd = iOriginalStart+rowFactor-1;    
        for j=1:reducedSize(2)
            jOriginalStart = (j-1)*columnFactor+1;
            jOriginalEnd = jOriginalStart+columnFactor-1;
            reducedmatrix(i,j,:) = sum(sum(matrix(iOriginalStart:iOriginalEnd,jOriginalStart:jOriginalEnd,:).*includedElements,2),1);
        end
    end
    if numel(sizeOfMatrix) <= 2
        matrix=squeeze(matrix);
    end
end