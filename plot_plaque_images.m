


load(['~/Documents/Images/plaque_' num2str(plaque_name) '_' plaque_date '/plaque_' ...
    num2str(plaque_name) '_reconstructedPhC.mat']);

close all

FigHandle = figure;
set(FigHandle, 'Position', [800, 500, 1049, 500]);


subplot(1,2,1)
colormap gray
imagesc(reconstructedImages(rows,columns,1,1))
axis image
subplot(1,2,2)
colormap gray
imagesc(reconstructedImages(rows,columns,1,2))
axis image


% figure
% colormap gray
% imagesc(abs(reconstructedImages(rows,columns,1,2)-reconstructedImages(3000:4000,1500:3000,1,1)))