%myreadtable = @(filename)readtable(filename, 'Format', '%s%s%s %f%f%f%f%f %s%s %f%s%f %s%s%s%s%s %f'); %case with 5 tags
%metaT =  webread('https://ifcb-data.whoi.edu/api/export_metadata/NESLTER_transect', weboptions('Timeout', 60, 'ContentReader', myreadtable));
opts = delimitedTextImportOptions("NumVariables", 21);
opts.DataLines = [2 inf];
%case with 5 tags
opts.VariableTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "string", "string", "double", "string", "double", "string", "string", "string", "string", "string", "string", "double", "double"];
opts.VariableNamesLine = 1;
myreadtable = @(filename)readtable(filename, opts);
metaT = webread('https://ifcb-data.whoi.edu/api/export_metadata/NESLTER_transect', weboptions('Timeout', 60, 'ContentReader', myreadtable));
  
mdate = datenum(metaT.sample_time, 'yyyy-mm-dd HH:MM:ss+00:00');
ind = find(strcmp(metaT.sample_type, 'underway') | ismissing(metaT.sample_type));
mdate = mdate(ind);
%FUDGE for AR77 until data uploaded
mdate = [mdate; datenum(2023,10,11:16)'];
%%
%case to include SPIROPA north of 39.5 latitude
opts = delimitedTextImportOptions("NumVariables", 20);
opts.DataLines = [2 inf];
opts.VariableTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "string", "string", "double", "string", "double", "string", "string", "string", "string", "string", "double", "double"];
opts.VariableNamesLine = 1;
myreadtable = @(filename)readtable(filename, opts); %
metaT =  webread('https://ifcb-data.whoi.edu/api/export_metadata/SPIROPA', weboptions('Timeout', 60, 'ContentReader', myreadtable));
mdate2 = datenum(metaT.sample_time, 'yyyy-mm-dd HH:MM:ss+00:00');
ind = (strcmp('underway', metaT.sample_type) & metaT.latitude>=39.5);
mdate = [mdate; mdate2(ind)];
%%

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
plot(yd, d(:,1), 'r.', 'markersize',20)
xlim([1 365])
datetick('x', 'm', 'keeplimits')
ylim([min(d(:,1))-.2 max(d(:,1))+.2])
set(gca, 'ydir', 'rev', 'ytick', min(d(:,1)):max(d(:,1)))
set(gcf, 'position', [550 200 300 340])

unique(metaT.cruise(~ismissing(metaT.cruise)))