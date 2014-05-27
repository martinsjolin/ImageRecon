reconstruct_image_tabletop_setup
Mats P oct 2013

First run packageMeasuredDataStack to read the raw data files and save them as .mat files. The data needed for the reconstruction is called _measuredData.mat.
The other generated files are intermdiate files that could be useful for debugging but are not needed for the reconstruction.
If purePhotonCounting is true, the eight bins are summed into one bin, and the generated files get the ending PhC.
To get the right settings for different images, look at the packageMeasuredDataStack files located on the file server in the same folders as the data from scans with the lab setup.
Note that these files may have to be modified by replacing

tubeSpectrum = computeStandardTubeSpectrum(energies, kVp, 1); %The normalization is unimportant here
detectionEfficiency = computeStandardDetectionEfficiency(energies);

with

tubeSpectrum = NaN*ones(size(energies)); %WORKAROUND, since tubeSpectrum is not used anyway
detectionEfficiency = NaN*ones(size(energies)); %WORKAROUND, since detectionEfficiency is not used anyway

The variables tubeSpectrum and detectionEfficiency are normally not used anywayso I did not include computeStandardTubeSpectrum and computeStandardDetectionEfficiency
in the package.

When packageMeasuredDataStack is finished, reconstruct by running reconstructMeasuredStackWithBackprojector.
When running this, two variables must already be defined in the workspace: inFilename and outFilename
These tell the program which .mat file to load data from (normally ending in _measuredData.mat or _measuredDataPhC.mat) and which .mat file to save the reult in (normally ending in _reconstructed.mat or _reconstructedPhC.mat).

There is a known bug that causes problems if the registered angles in geometryParameters.angles are far from 0 (even though changing them by n*2*pi should not matter).
The workaround for this is to add a constant to geometryParameters.angles in the beginning of reconstructMeasuredStackWithBackprojector, making geometryParameters.angles lie close to 0 (I'm not sure but I think it should be roughly between 0 and 2*pi).

The end result is a variable called reconstructedImages. If purePhotonCounting is true, it is an image. If purePhotonCounting is false, it is a stack of eight images, one for each bin.

The function weight.m is useful to make a linear combination of the bin images with specified weights.
exportPNG.m can be used to export an image matrix to a PNG image file.