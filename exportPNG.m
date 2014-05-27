function exportPNG(image, outFilename, window)
    %function exportPNG(image, outFilename, window)
    %16 bit grayscale png export
    %image = matrix containing te image to be exported
    %outFilename is the output file name
    %window is a 1x2 vector containing the lower and upper limits (in that
    %order) of thegrayscale window.
    %If no window is specified, the window is chosen so that the  image
    %fits precisely in the window
    %RGB images can also be exported. For this purpose, make sure to
    %specify the window, which is then assumed to apply to each color
    %component separately.
    if (nargin<3)
        %If window is not specified: use default window
        window = [min(min(image)) max(max(image))];
    end
    %Convert to uint16 and save as 16-bit grayscale. The conversion from
    %double to uint16 is carried out with the mapping ax+b where a and b are
    %constants. windowMin is mapped onto 0 and windowMax is mapped onto 65535.
    a = 65535/(window(2)-window(1));
    b = -a*window(1);

    uint16Image = uint16(a*image+b);
    imwrite(uint16Image, outFilename,'png');
end