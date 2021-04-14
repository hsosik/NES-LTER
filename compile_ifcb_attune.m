% AR29 compile IFCB and Attune data, calculate size structure statistics

clear; close all;

% load metadata

attnpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20190705_TN368\bead_calibrated\class\';
ifcbpath = '\\sosiknas1\IFCB_products\SPIROPA\features\D2019\';

load \\sosiknas1\Lab_data\Attune\cruise_data\20190705_TN368\Summary\Attune_uw_match
attnfiles = Attune_uw_match.Filename(good);
attnvolanalyzed = Attune_uw_match.VolAnalyzed_ml(good);
attntime = datenum(Attune_uw_match.StartDate(good));

myreadtable = @(filename)readtable(filename, 'Format', '%s%s%s %f%f%f%f%f %s%s %f%s%f %s%s%s%s%s %f');
metaT =  webread('https://ifcb-data.whoi.edu/api/export_metadata/SPIROPA', ...
    weboptions('Timeout', 60, 'ContentReader', myreadtable));
%metaid = strcmp(metaT.cruise, 'AR29') & strcmp(metaT.sample_type, 'underway') & metaT.ifcb==127;
metaid = strcmp(metaT.cruise, 'TN368') & strcmp(metaT.sample_type, 'underway');
ifcbfiles = metaT.pid(metaid);
ifcbvolanalyzed = metaT.ml_analyzed(metaid);
ifcbtime = datenum(metaT.sample_time(metaid), 'yyyy-mm-dd HH:MM:SS+00:00');

load \\sosiknas1\IFCB_data\SPIROPA\match_up\SPIROPA_TN368__newdashboard_USEMEuw_match

% pre-allocate data structures

bins = logspace(log10(0.5), log10(150), 128);
ifcbll = 10.1344; % IFCB lower limit and Attune upper limit line up with bin edges
attnul = 15.1824;
counts = NaN(length(ifcbfiles), length(bins)-1); % counts per mL
biovol = counts; % total biovolume per mL
attnvol = NaN(length(ifcbfiles), 1); % total sample volume
ifcbvol = attnvol;

% compile data from IFCB and Attune

for i = 1:length(ifcbfiles)
    
if rem(i,10)==0
    fprintf('%i of %i \n', i, length(ifcbfiles))
end

if isfile([ifcbpath ifcbfiles{i}(1:9) filesep ifcbfiles{i} '_fea_v4.csv'])
    ifcbdata = readmatrix([ifcbpath ifcbfiles{i}(1:9) filesep ifcbfiles{i} '_fea_v4.csv']);
    ifcb_bv = ifcbdata(:,22)/21.2539;
    ifcb_esd = 1.2407*ifcb_bv.^(1/3);
    ifcbvol(i) = ifcbvolanalyzed(i);
    attni = find(attntime > ifcbtime(i)-0.0069 & attntime < ifcbtime(i)+0.0069); % +/- 10 mins
    attn_bv = [];
    attn_esd = [];
    sample_vol = 0;
    for j=1:length(attni)
        if isfile([attnpath strrep(attnfiles{attni(j)}, '.fcs', '.mat')])
            load([attnpath strrep(attnfiles{attni(j)}, '.fcs', '.mat')]);
            volume(class==0 | class==7) = []; % remove junk/noise
            attn_bv = [attn_bv; volume];
            attn_esd = [attn_esd; 1.2407*volume.^(1/3)];
            sample_vol = sample_vol + attnvolanalyzed(attni(j));
        end
    end
    attnvol(i) = sample_vol;
    % compile counts per mL
    for j = 1:length(bins)-1
        if bins(j) < ifcbll
            counts(i,j) = sum(attn_esd>=bins(j)&attn_esd<bins(j+1))/attnvol(i);
        elseif bins(j) >= ifcbll && bins(j) <= attnul
            counts(i,j) = (sum(attn_esd>=bins(j)&attn_esd<bins(j+1)) + ...
                sum(ifcb_esd>=bins(j)&ifcb_esd<bins(j+1)))/(attnvol(i) + ifcbvol(i));
        else
            counts(i,j) = sum(ifcb_esd>=bins(j)&ifcb_esd<bins(j+1))/ifcbvol(i);
        end
    end
    % compile total biovolume per mL
    for j = 1:length(bins)-1
        if bins(j) < ifcbll
            biovol(i,j) = sum(attn_bv(attn_esd>=bins(j)&attn_esd<bins(j+1)))/attnvol(i);
        elseif bins(j) >= ifcbll && bins(j) <= attnul
            biovol(i,j) = (sum(attn_bv(attn_esd>=bins(j)&attn_esd<bins(j+1))) + ...
                sum(ifcb_bv(ifcb_esd>=bins(j)&ifcb_esd<bins(j+1))))/(attnvol(i) + ifcbvol(i));
        else
            biovol(i,j) = sum(ifcb_bv(ifcb_esd>=bins(j)&ifcb_esd<bins(j+1)))/ifcbvol(i);
        end
    end
end

% append metadata rows
matchid = find(strcmp(ifcbfiles(i), IFCB_match_uw_results.pid));
if  exist('metadata')
    metadata = [metadata; IFCB_match_uw_results(matchid, :)];
else
    metadata = IFCB_match_uw_results(matchid, :);
end

end

% remove NaN rows (IFCB files with no corresponding Attune files)

nani = isnan(counts(:,1));
counts(nani,:) = [];
biovol(nani,:) = [];
attnvol(nani) = [];
ifcbvol(nani) = [];
metadata(nani,:) = [];

% compute total biovolume in each size class

pico = NaN(length(biovol), 1);
nano = pico;
micro = pico;
for i=1:length(biovol)
    pico(i) = sum(biovol(i,bins<2));
    nano(i) = sum(biovol(i,bins>=2&bins<=20));
    micro(i) = sum(biovol(i,bins>20&bins<150));
end
total_bv = pico + nano + micro;

% match samples to cruise transect number

run('rb1904_transect_startime')
tnum = NaN(height(metadata), 1);
for i=1:length(tstime)-1
    tnum(metadata.mdate >= tstime(i) & metadata.mdate < tstime(i+1)) = i;
end
metadata.transect = tnum;

save('TN368_compiled_ifcb_attune_data', 'bins', 'counts', 'biovol', 'attnvol', 'ifcbvol', 'metadata', ...
    'pico', 'nano', 'micro', 'total_bv');
