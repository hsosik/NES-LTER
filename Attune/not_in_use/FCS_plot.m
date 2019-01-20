function  [ FCSplot ] = FCS_plot(fcs_path, index_x, index_y)
%'\\sosiknas1\Lab_data\Attune\EN608\ExportedFCS'
%output: figure with PE signal graph

[~,fcshdr,fcsdatscaled] = fca_readfcs(fullfile(fcs_path, fcslist{ii}));

if ~exist('fcs_path', 'var')
    fcs_path = uigetdir(pwd, 'Pick a Directory of FCS files')
else
    if ~exist(fcs_path, 'dir')
        disp('WARNING: Directory not found. Check input path.')
        FCSPEplot = [];
        return
    end
end

fcslist = dir(fullfile(fcs_path, '*.fcs'));
fcslist = {fcslist.name}';
FCSPEplot.x = NaN(size(fcslist));
FCSPEplot.y = FCSPEplot.x;
index_x = index_x;
index_y = index_y;
for ii = 1:length(fcslist)
    if ~rem(ii,10)
        disp([num2str(ii) ' of ' num2str(length(fcslist))])
    end
    [~, fcsdatscaled] = fca_readfcs(fullfile(fcs_path, fcslist{ii}));
    FCSplot.x(ii)= fcsdatscaled(:,index_x);
    FCSplot.y(ii) = fcsdatscaled(index_y);
end

FCSPEplot.filelist = fcslist;

loglog(fcsdatscaled(:,11),fcsdatscaled(:,19),'.')
xlim([10^2  10^6])
ylim([10^2  10^6])
xlabel(fcshdr.par(index_x).name)
ylabel(fcshdr.par(index_y).name)
title('PE Signal for $\Syn$')

end 
