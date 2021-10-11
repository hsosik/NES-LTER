myreadtable = @(filename)readtable(filename, 'Format', '%s%s%s %f%f%f%f%f %s%s %f%s%f %s%s%s%s%s %f'); %case with 5 tags
metaT =  webread('https://ifcb-data.whoi.edu/api/export_metadata/NESLTER_transect', weboptions('Timeout', 60, 'ContentReader', myreadtable));
mdate = datenum(metaT.sample_time, 'yyyy-mm-dd HH:MM:ss+00:00');
ind = strmatch('underway', metaT.sample_type);
ind2 = strmatch(' ', char(metaT.sample_type));
ind = [ind; ind2]; clear ind2

%uw = load('c:\work\LTER\uw_all');
%uw.uw_sum = [];
%for ii = 1:length(uw.cruises)
%    uw.uw_sum = [uw.uw_sum; datenum(uw.uw_all{ii}.date,'yyyy-mm-dd hh:MM:ss') uw.uw_all{ii}.(uw.vname_uw{ii}{3}) uw.uw_all{ii}.(uw.vname_uw{ii}{4})  uw.uw_all{ii}.(uw.vname_uw{ii}{1}) uw.uw_all{ii}.(uw.vname_uw{ii}{2})];
%end
%d = datevec(uw.uw_sum(:,1));
%yd = uw.uw_sum(:,1) - datenum(d(:,1),1,0);

d = datevec(mdate);
yd = mdate-datenum(d(:,1),1,0);

figure
plot(yd(ind), d(ind,1), 'r.', 'markersize',20)
datetick('x', 'm')
ylim([min(d(ind,1))-.2 max(d(ind,1))+.2])
set(gca, 'ydir', 'rev', 'ytick', min(d(:,1)):max(d(:,1)))
set(gcf, 'position', [550 200 300 340])

