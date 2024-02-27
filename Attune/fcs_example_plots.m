fullfilename = '\\sosiknas1\Lab_data\Attune\cruise_data\20200201_EN649\FCS\NESLTER_EN649_03Feb2020c_Group_day0_Sample(14).fcs';

[~,fcshdr,fcsdatscaled] =fca_readfcs(fullfilename);
t = fcsdatscaled(:,2:10);
n = {fcshdr.par(2:10).name};
t2 = fcsdatscaled(:,11:19);
n2 = {fcshdr.par(11:19).name};

cmap = [1 0 0;  2/3 0 1; 0 0 1; 1 0.5 0; 1 1 0; 0 1 0; 0 0 1; 0 0.5 0 ];
cmap = [1 0 0;  2/3 0 1; 0 0 1; 1 0.5 .5; 1 1 0; 0 1 0; 0 0 1; 0 0.5 0 ];
cmap = [1 0 0;  2/3 0 1; 0 0 1; 2/3 0.5 0; 1 1 0; 0 1 0; 0 0 1; 0 0.5 .5 ]; 
[h,ax,BigAx,hhist,pax] = plotmatrix_scatter(log10(t2),class, '.');
colormap(cmap)
for ii = 1:length(n), ax(1,ii).XLabel.String=n2{ii}; set(ax(1,ii),'xaxislocation', 'top'); end
for ii = 1:length(n), ax(ii,1).YLabel.String=n2{ii}; end

