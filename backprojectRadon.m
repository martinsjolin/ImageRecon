function res = backprojectRadon(interpolatedProjection, interpolationAngles, outputSize)
    %It doesn't seem like one can specify the output scaling in iradon.
    %Passing an outputSize parameter only changes the field of view.
    
    res = iradon(interpolatedProjection, interpolationAngles*180/pi, 'linear','Ram-Lak', 1, outputSize(1));
%    res = iradon(interpolatedProjection, interpolationAngles*180/pi, 'linear','Ram-Lak', 1, 4000);

    %iradon wants angles in degrees.
    %If one wants to get a scaled output image, remove outputSize in the
    %above call to iradon and uncomment the following line:
    %res = imresize(res, outputSize);
    %This step could be avoided by writing a better backprojection
    %algorithm.
end