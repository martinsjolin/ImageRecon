% run_plaque_reconstruction

clf,clear, close all

plaque_name = '763';
plaque_date = '2013_12_17';

slice_nr = 1;
calibration_number = 1;

packageMeasuredDataStack
%%
reconstructMeasuredStackWithBackprojector

rows = 1:5000;
columns = 1:5000;

plot_plaque_images

%%

close all

plaque_name = '764';
plaque_date = '2013_12_18';
calib_number = 0;

load(['~/Documents/Images/plaque_' plaque_name '_' plaque_date '/plaque_' plaque_name '_calibrationRawDataPhC_' num2str(calib_number) '.mat'])

imagesc(squeeze(sum(calibrationRegisteredFaultyData,2))/1200)