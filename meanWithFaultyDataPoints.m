function meanNonfaulty = meanWithFaultyDataPoints(data, faulty)

    %function meanNonfaulty = meanWithFaultyDataPoints(data, faulty)
    %Calcuates the mean of data with the points marked with ones in
    %"faulty" excluded. The mean is taken over all dimensions and is
    %therefore a scalar.
    %
    %data and faulty are n-dimensional matrices which must have the same
    %size. faulty contains zeros and ones. A one means that the
    %corresponding entry in data is faulty.
    
    nDataPoints = numel(data);
    if(numel(faulty)~=nDataPoints)
        error('Error in meanWithFaultyDataPoints: size of data and faulty incompatible.')
    end
    
    dataVector = reshape(data, 1, nDataPoints);
    faultyVector = reshape(faulty, 1, nDataPoints);
    
    nonfaultyData = dataVector(find(~faultyVector));
    meanNonfaulty = mean(nonfaultyData);
end