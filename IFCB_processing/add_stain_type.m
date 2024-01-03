% ToTag_xlsFile = '\\sosiknas1\IFCB_data\SPIROPA\to_tag\SPIROPA_AR29_to_tag_IFCB14staining_needsSTAININGediting.xls';
ToTag_xlsFile = '\\sosiknas1\IFCB_data\NESLTER_broadscale\to_tag\NESLTER_broadscale_PC1405_to_tag_IFCB014_needsstaining.xls';

totag = readtable(ToTag_xlsFile);
% urlbase = 'https://ifcb-data.whoi.edu/SPIROPA/';
urlbase = 'https://ifcb-data.whoi.edu/NESLTER_broadscale/';

%runType = cell(size(totag,1),1);
for count = 1:size(totag,1)
    h = IFCBxxx_readhdr([urlbase char(totag.Filename{count}) '.hdr']);
    runType{count} = h.runtype;
end

ii = strmatch('NORMAL', runType);
totag.tag2(ii) = {'unstained'};

ii = strmatch('ALT', runType);
totag.tag2(ii) = {'FDA_stained'};

f = regexprep(ToTag_xlsFile, '.xls', '_with_runtype.csv');

writetable(totag, f);

