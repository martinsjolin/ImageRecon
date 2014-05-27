function [calibrationParameters confIntHalfWidths] = createCalibrationParametersSectraStyle(calibrationRegistered, calibrationRegisteredFaultyData, pValuesInCalibration)

    %Possibility for improvement: use weighted least squares.

    %Watch out -  things are changing rapidly so here comes a description
    %of the model and input and output.
    
    %The model: p = f(N)
    %N = registered counts in the current pixel
    %f is given implicitly by N = exp(a*p^2+b*p+c)
    %This means tht there are two possible values for p given N but
    %this should not be  problem in practice - one of the solutions should
    %be way outside the range of values we're interested in.
    %Note that N_0 is is not icluded in the model. Altering N_0 only alters
    %c. However, If there are flux variations one might want to multiply N
    %(i.e. calibrationRegistered) with some factor k taking this into
    %account before passing it as a parameter to this function.
    
    %Calibration is done on a detector element basis - i.e. the calibration
    %factors are assumed to be the same for all translations of the object.
    
    %Input:
    
    %calibrationRegistered contains all measurements for _one_ calibration
    %sequence.
    %The exposure time does not matter here.
    %Dimensions:
    %nDetectorElements, nCalibrationMeasurements, nBins, nCalibrationThicknesses
    
    %calibrationRegisteredFaultyData has the same size as
    %calibrationRegistered and contains ones in positions where the
    %measurement is wrong (e.g. because of overflow), and zeros otherwise.

%     This one turned out to be unnecessary
%     %unattenuatedMeanRegisteredInCalibration is a matrix with dimensions
%     %(nBins,nCalibrationThicknesses) elements. It contains the N_0 values i.e. the
%     %average count rates over all detector elements, for the calibration
%     %measurements. _Caution:_ The entries _all_ refer to unattenuated beam
%     %measurements, for all values of the last idex! The reason why there is
%     %a calibration thickness index is that one may in certain situations
%     %want to measure the unattenuated beam between each value of the
%     %calibrator thickness in order to counter the effect of drift.
     
    %pValuesInCalibration (dimensions: nBins, nCalibrationThicknesses)
    %gives the p Values of the calibration measurements. These can be
    %calculated from measured calibrator thicknesses of from the count
    %numbers themselves.

    %Output:
    %calibrationParameters, a matrix with the parameters a, b and c for
    %all detector elements and all bins.
    %dimensions: nDetectorElements, nBins, nParameters
    %the last index selects between the parameters a (1) b (2) and c (3).

    %Also there is an optional output argument confIntHalfWidths that gives
    %the half-widths of the 95 % confidence intervals of the estimates of p
    %(of the fitted curves) at the points pValuesInCalibration. (This
    %particular choice of p values at which to return the confidence
    %intervals is quite arbitrary - any set of p values could have been
    %chosen.)
    %dimensions: nDetectorElements, nBins, nCalibrationThicknesses
    
    DEBUG_MODE = false; %Can only be set to true when stepping through the code manually in debug mode.
    
    nDetectorElements = size(calibrationRegistered,1);
    nCalibrationMeasurements = size(calibrationRegistered,2);
    nBins = size(calibrationRegistered,3);
    nCalibrationThicknesses = size(calibrationRegistered,4);
    
    nParameters=3;
    
    calibrationParameters = zeros(nDetectorElements, nBins, nParameters);
    if(nargout == 2)
        confIntHalfWidths=zeros(nDetectorElements,nBins,nCalibrationThicknesses);
    end
    
    for detectorElementNo = 1:nDetectorElements
        for binNo = 1:nBins
            p = (repmat(pValuesInCalibration(binNo,:),1,nCalibrationMeasurements))';
            logN=squeeze(log(calibrationRegistered(detectorElementNo,:,binNo,:)));
            faultyDataPoints = reshape((squeeze(calibrationRegisteredFaultyData(detectorElementNo,:,binNo,:)))',nCalibrationMeasurements*nCalibrationThicknesses,1);
            A = [(p(~faultyDataPoints)).^2 p(~faultyDataPoints) ones(length(find(~faultyDataPoints)),1)]; 
            bWithFaultyData = reshape(logN',nCalibrationMeasurements*nCalibrationThicknesses,1);
            b = bWithFaultyData(~faultyDataPoints);
            %tic; %Factor 30 speedup compared to nlinfit
            calibrationParameters(detectorElementNo, binNo,:) = A\b;
            %The above line calculates a, b and c by least squares fitting.
            %toc;
            
            if(nargout == 2)
                %Calculate confidence bounds as well. Note: We calculate the
                %curve by least-squares but the confidence intervals
                %correspond to another curve which is calculated with
                %nlinfit. In practice this does not make any difference
                %since they seem to agree to very high (roundoff?)
                %precision.
                [coeffs, R, J, covb, mse] = nlinfit(p(~faultyDataPoints),b,@f,[0 -1 5]);
                [logNPred, delta] = nlpredci(@f,pValuesInCalibration(binNo,:),coeffs,R,'covariance',covb);
                confIntHalfWidths(detectorElementNo,binNo,:)=delta./abs(2*coeffs(1)*pValuesInCalibration(binNo,:)+coeffs(2))';
            end
            
            %----------------------Test section----------------------------
            if(DEBUG_MODE==true)
                disp('Test code executed in createCalibrationParametersSectraStyle.')
                %A check to see if the confidence intervals are believable
                
                %I'm confused abut the cofidence intervals. The error bars
                %seem to be a litle bit to small. Why?
                hold off;
                for i = 1:100
                    index = 1:floor(length(p)/101):length(p);
                    nonFaultyP=p(~faultyDataPoints);
                    [coeffs, R, J, covb, mse] = nlinfit(nonFaultyP(index(i):index(i+1)),b(index(i):index(i+1)),@f,[0 -1 5]);
                    plot(b(index(i):index(i+1)),nonFaultyP(index(i):index(i+1))+0.5e-2*randn(size(index(i):index(i+1)))','g+')
                    hold on;
                    logNForPlotting = 6:0.05:6.6;
                    for j = 1:length(logNForPlotting)
                        pForPlotting(j) = fsolve(@(p) f(coeffs, p)-logNForPlotting(j),2);
                    end
                    [logNPred, delta] = nlpredci(@f,pForPlotting,coeffs,R,'covariance',covb);
                    if(i==1)
                        errorbar(logNPred,pForPlotting,delta./abs(2*coeffs(1)*pForPlotting+coeffs(2))');
                    else
                        plot(logNPred,pForPlotting,'y')
                    end
                end
            end
            if(DEBUG_MODE==true)
                disp('Test code executed in createCalibrationParametersSectraStyle.')
                hold off;
                plot(b,p(~faultyDataPoints)+0.5e-2*randn(size(p(~faultyDataPoints))),'g+')
                % ^^ with small random offsets for better visual appearance
                plot(b,p(~faultyDataPoints),'g+')
                hold on;
                pForPlotting = pValuesInCalibration;
                disp('a b and c calculated by least squares fitting')
                calibrationParameters(detectorElementNo, binNo,:)
                plot(calibrationParameters(detectorElementNo,binNo,1)*pForPlotting.^2+calibrationParameters(detectorElementNo, binNo,2)*pForPlotting+calibrationParameters(detectorElementNo, binNo,3),pForPlotting,'r')
                xlabel('log(N)')
                ylabel('p')
                %pause;
                
                %Also try curve fitting with nlinfit
                %The starting guess shouldn't matter but check anyway
                tic; [coeffs, R, J, covb, mse] = nlinfit(p(~faultyDataPoints),b,@f,[0 -1 5]); toc;
                disp('a b and c calculated by nlinfit')
                coeffs
                disp('mean square error:')
                mse
                [logNPred, delta] = nlpredci(@f,pForPlotting,coeffs,R,'covariance',covb);
                errorbar(logNPred,pForPlotting,delta./abs(2*coeffs(1)*pForPlotting+coeffs(2))');
            end
            %-----------------End of test section--------------------------
            
        end
    end
end

function res = f(c,p)
    res= c(1)*p.^2+c(2)*p+c(3);
end