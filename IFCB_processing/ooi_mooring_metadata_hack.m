t = dir('\\sosiknas1\IFCB_data\ooi\pioneer_mooring\data\\**\*.roi');
p = {t.name}';
meta_data = table;
meta_data.pid = regexprep({t.name}', '.roi', '');
meta_data.sample_time = datetime(IFCB_file2date(meta_data.pid), 'ConvertFrom', 'datenum', 'Format', 'uuuu-MM-dd'' ''HH:mm:ss+00:00');
meta_data.ml_analyzed = IFCB_volume_analyzed(strcat({t.folder}', filesep, regexprep({t.name}', 'roi', 'hdr')));
meta_data.sample_type

bd = dir('\\sosiknas1\IFCB_data\ooi\pioneer_mooring\data\beads\*.roi');
[a,b,c] = intersect(regexprep({b.name}', '.roi', ''),meta_data.pid);
meta_data.sample_type(:) = {'mooring'};
meta_data.sample_type(c) = {'bead'};
meta_data.skip(:) = 0;

save('\\sosiknas1\IFCB_products\ooi\pioneer_mooring\summary\meta_data', 'meta_data')
