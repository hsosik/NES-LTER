%Calibration Curve

%extract mode from SSC-A
filepath = 'E:\Attune_Data\Size_calibration_July2018\ExportedFCS';

filelist = dir(filepath)

for ii in length(filelist)
[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(fullfile(filepath,filename))

