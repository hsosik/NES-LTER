function [] = IFCB_writesb_withFeatures_withClass_NESLTER(dataset,sampletype,yr)
%
% IFCB_writesb_withFeatures_withClass_NESLTER('NESLTER_broadscale','underway',2025)
% IFCB_writesb_withFeatures_withClass_NESLTER('NESLTER_transect','cast',2024)
%
% Here are some USER choices
%dataset = 'NESLTER_transect'; %NESLTER_transect NESLTER_broadscale
%sampletype = 'discrete'; %'underway' 'discrete'
%yr = 2025;

%%
disp('Loading metadata...')
metaT = readtable(['https://ifcb-data.whoi.edu/api/export_metadata/' dataset]);
metaT.datetime = datetime(metaT.sample_time, 'InputFormat','uuuu-MM-dd HH:mm:ss+00:00');
%%
%handle lack of consistent sample depth info in ifcb-dashboard metadata for ships
ship = table;
ship.name =         {'AR' 'AT' 'EN' 'HB' 'RB' 'HRS' 'TN' 'OE' 'GU' 'AE' 'PC'}';
ship.intake_depth = [2     5    5    5    5    3     4.5  4.5  4.5  1.5  5]';
%%
outpath_base = ['\\sosiknas1\IFCB_products\' dataset '\SeaBASS\20220209_Jan2022_NES_2.4\sb_files\'];
feabase = ['\\sosiknas1\IFCB_products\' dataset '\features\'];
classbase = ['\\sosiknas1\IFCB_products\' dataset '\class\v3\20220209_Jan2022_NES_2.4\'];
urlbase = ['http://ifcb-data.whoi.edu/' dataset '/'];
r2rbase = 'https://nes-lter-api.whoi.edu/api/events/';
load(['\\sosiknas1\IFCB_products\' dataset '\SeaBASS\20220209_Jan2022_NES_2.4\automated_assessed_ID_table_NES_2.4.mat']); %created with make_SB_assessed_ID_table.m
%%
if exist('hdr', 'var')
    clear hdr*
end
%ship = 'survey'; %'survey' 'process'

hdr.investigators = 'Heidi_Sosik';
hdr.affiliations = 'Woods_Hole_Oceanographic_Institution';
hdr.contact = 'hsosik@whoi.edu';
hdr.documents = ['IFCB_brief_protocol_Sosik_Jul2026.docx,checklist_IFCB_plankton_and_particles_NES-LTER_' dataset '_R1_Sosik_automated_classification.docx,namespace_ptwg_nonconforming_roi_v1.yml,automated_assessed_ID_table_NES_2.4.txt'];

switch sampletype
    case  'underway' %
        sbdatatypestr = 'flow_thru';
        IFCBtypestr = {'underway'};
    case 'discrete'
        sbdatatypestr = 'bottle'; %default, overwrite with net_tow as appropriate below
        IFCBtypestr = {'cast' 'underway_discrete'};
end
ii = find(~metaT.skip & year(metaT.datetime) == yr & ismember(metaT.sample_type, IFCBtypestr) & ~ismember(metaT.ifcb, [9 14]));
cruise_list = unique(metaT.cruise(ii));
outpath = [outpath_base num2str(yr) filesep];
disp('Loading R2R event logs if needed...')
for tt = 1:length(cruise_list)
    c = cruise_list{tt};
    if ~exist([outpath c], 'dir')
        mkdir([outpath c])
    end
    if strcmp(sampletype, 'discrete') %get the event logs
        r2r.(c) = readtable([r2rbase c '.csv'], 'VariableNamingRule','preserve');
        if iscell(r2r.(c).Cast)
            ind = strcmp(r2r.(c).Instrument, 'CTD911');
            cast_temp = NaN(size(r2r.(c).Cast));
            cast_temp(ind) = str2num(char((r2r.(c).Cast(ind))));
            r2r.(c).Cast = cast_temp;
        end
    end
end

hdr.experiment = 'NES-LTER';
hdr.cruise = 'temp';
hdr.r2r_event = 'NA'; %default for underway
hdr.data_file_name = 'temp';
hdr.calibration_files = 'no_cal_files';
hdr.eventID = 'temp';
hdr.data_type = sbdatatypestr;
hdr.instrument_model = 'temp';
hdr.instrument_manufacturer = 'McLane_Research_Laboratories_Inc';
hdr.data_status = 'final';
hdr.start_date = 'temp';
hdr.end_date = 'temp';
hdr.start_time = 'temp';
hdr.end_time = 'temp';
hdr.north_latitude = 'temp';
hdr.south_latitude = 'temp';
hdr.east_longitude = 'temp';
hdr.west_longitude = 'temp';
hdr.water_depth = 'NA';
hdr.measurement_depth = 'temp';
hdr.volume_sampled_ml = '5';
hdr.volume_imaged_ml = 'temp';
hdr.pixel_per_um = '2.77';
hdr.associatedMedia_source = 'temp';
%hdr.associated_archives = [hdr.experiment '-' hdr.cruise '_' cruise '_IFCB_raw_data_associated.tgz'];
hdr.associated_archives = [hdr.experiment '_' dataset '_' num2str(yr) '_IFCB_raw_data_associated.tgz'];
hdr.associated_archive_types = 'raw';
hdr.length_representation_instrument_varname = 'maxFeretDiameter';
hdr.width_representation_instrument_varname = 'minFeretDiameter';
hdr.missing = '-9999';
hdr.delimiter = 'comma';
hdr.fields = 'associatedMedia,data_provider_category_automated,scientificName_automated,scientificNameID_automated,prediction_score_automated_category,biovolume,area_cross_section,length_representation,width_representation,equivalent_spherical_diameter';
hdr.units = 'none,none,none,none,none,um^3,um^2,um,um,um';

ii = find(~metaT.skip & year(metaT.datetime) == yr & ismember(metaT.sample_type, IFCBtypestr) & ~ismember(metaT.ifcb, [9 14]));

comment1 = ['! ' dataset ' '  num2str(yr)];
comment2 = '! To access each image directly from the associatedMedia string: replace .html with .png';
comment3 = '! v4 ifcb-analysis image products; https://github.com/hsosik/ifcb-analysis';
%comment6 = '! Files updated (R2 version) to include taxonomic classification';

pixel_per_micron = str2num(hdr.pixel_per_um);
hdr_copy = hdr;
%%
for count = 1:length(ii) %1:round(length(ii)/5):length(ii)
    hdr = hdr_copy;
    m = metaT(ii(count),:);
    hdr.data_type = sbdatatypestr; %default
    hdr.cruise = char(m.cruise);
    % switch char(m.sample_type) %all other cases as defaulted to 'flow_thru' or 'bottle'
    %     case 'underway_discrete'
    %         hdr.data_type = 'flow_thru';
    %         hdr = rmfield(hdr,'data_use_warning');
    %     case 'cast'
    %         hdr = rmfield(hdr,'data_use_warning');
    % end
    comment5 = '!';
    %comment6 = '!';
    if strcmp(sampletype, 'discrete') % for this case, check for type and if concentrated'
        comment5 = ['! ' char(m.sample_type)];
        if strcmp(m.sample_type, 'cast')
            comment5 = [comment5 '  ' num2str(m.cast) ' Niskin ' num2str(m.niskin)];
            ind = r2r.(hdr.cruise).Cast==m.cast & strcmp(r2r.(hdr.cruise).Instrument, 'CTD911');
            if ~startsWith(r2r.(hdr.cruise).R2R_Event{ind}, hdr.cruise)
                hdr.r2r_event = [lower(hdr.cruise) '-SE-' r2r.(hdr.cruise).R2R_Event{ind}];
            end
        end
        %      if ~ismissing(m.comment_summary)
        %          comment6 = [comment6 ' ' char(m.comment_summary)];
        %      end
        % %            if strmatch(m.sample_type, 'cast')
        % %                if isnumeric(r2r.Cast(1))
        % %                    ind = find(strcmp(r2r.Instrument, 'CTD911') &  (str2num(char(m.cast)) == r2r.Cast) & ~strcmp(r2r.Action, 'other'));
        % %                else
        % %                    ind = find(strcmp(r2r.Instrument, 'CTD911') &  (str2num(char(m.cast)) == str2double(regexprep(r2r.Cast, 'C', '')) & ~strcmp(r2r.Action, 'other')));
        % %                end
        % %                if isempty(ind)
        % %                    disp(strcat("cast: ", m.cast, " missing in R2R file"))
        % %                    keyboard
        % %                else
        % %                    eventlist = strcat(r2r.R2R_Event(ind), ','); eventlist = strcat(eventlist{:}); eventlist = eventlist(1:end-1);
        % %                    hdr.r2r_event = eventlist;
        % %                 %   %temp = regexprep(r2r.Station(ind), ' ', '');
        % %                 %   %hdr.station = temp{end};
        % %                 %   hdr.station = r2r.Station(ind(end));
        % %                 %   %if ~isequal(temp{:})%isequal(temp(ind(1)),temp(ind(2)))
        % %                 %   if sum(r2r.Station(ind)-r2r.Station(ind(1)))~=0
        % %                 %       disp('mismatched station info')
        % %                 %       disp(r2r.Station(ind))
        % %                 %       disp(hdr.station)
        % %                 %       keyboard
        % %                 %   end
        % %                end
        % %            end
        % % %       end
    end
    if ~strcmp(m.sample_type, 'underway_discrete') && ~strcmp(m.sample_type, 'underway')
        if isnan(m.depth) || m.depth==0
            hdr.measurement_depth = 'NA';
        else
            hdr.measurement_depth = num2str(m.depth);
        end
    elseif strcmp(m.sample_type, 'underway')
        hdr.measurement_depth = num2str(ship.intake_depth(strncmp(ship.name, hdr.cruise,2)));
    end
    hdr.eventID = char(m.pid);
    disp(hdr.eventID)
    hdr.data_file_name = [hdr.experiment '_' dataset '_'  num2str(yr) '_' hdr.cruise '_' sampletype '_IFCB_plankton_and_particles_' hdr.eventID([2:9 11:16]) '_R1.sb'];
    hdr.instrument_model = ['Imaging_FlowCytobot_IFCB' num2str(m.ifcb)];
    hdr.start_date = datestr(datenum(m.sample_time, 'yyyy-mm-dd HH:MM:ss+00:00'), 'yyyymmdd');
    hdr.end_date = hdr.start_date;
    hdr.start_time = [datestr(datenum(m.sample_time, 'yyyy-mm-dd HH:MM:ss+00:00'), 'HH:MM:ss') '[GMT]'];
    hdr.end_time = hdr.start_time;
    if isnan(m.latitude)
        hdr.north_latitude = 'NA';
        hdr.east_longitude = 'NA';
    else
        hdr.north_latitude = [num2str(m.latitude) '[DEG]'];
        hdr.east_longitude = [num2str(m.longitude) '[DEG]'];
    end
    hdr.south_latitude = hdr.north_latitude;
    hdr.west_longitude = hdr.east_longitude;
    hdr.volume_imaged_ml = num2str(m.ml_analyzed);
    hdr.associatedMedia_source = [urlbase char(m.pid) '.html'];
    f = fields(hdr);
    for eventlist = 1:length(f), hdrstr{eventlist} = ['/' f{eventlist} '=' hdr.(f{eventlist})]; end

    if m.trigger_selection == 2
        comment4 = ['! IFCB trigger mode: chlorophyll fluorescence (PMTB)'];
    elseif m.trigger_selection == 3
        comment4 = ['! IFCB trigger mode: chlorophyll fluorescence (PMTB) OR side scattering (PMTA)'];
    else
        comment4 = ['!'];
    end

    %hdrstr = ['/begin_header'; hdrstr(1:end-2)'; '!'; comment1; '!'; comment5; '!'; comment4; '!'; comment2; '!'; comment3; '!'; comment6; '!'; hdrstr(end-1:end)'; '/end_header'];
    hdrstr = ['/begin_header'; hdrstr(1:end-2)'; '!'; comment1; '!'; comment5; '!'; comment4; '!'; comment2; '!'; comment3; '!'; hdrstr(end-1:end)'; '/end_header'];
    outfullfile = [outpath hdr.cruise filesep hdr.data_file_name];
    %writecell(hdrstr, outfullfile, 'FileType', 'text', 'QuoteStrings', false);
    writelines(hdrstr, outfullfile, LineEnding="\n") %use this approach to force file to have on LF (UNIX style)
    clear hdrstr

    feafullname = [feabase hdr.eventID(1:5) filesep hdr.eventID(1:9) filesep hdr.eventID '_fea_v4.csv'];
    classfullname = [classbase hdr.eventID(1:5) filesep hdr.eventID(1:9) filesep hdr.eventID '_class.h5'];
    f = readtable(feafullname);
    c = load_class_scores(classfullname);
    [class_score,class_ind] = max(c.scores');
    class_label_data_provider = c.class_labels(class_ind);
    [~,label_ind] = ismember(class_label_data_provider, assessed_ID_table.data_provider_category_automated);
    outT = table;
    outT.associatedMedia = strcat(urlbase, hdr.eventID, '_', num2str(f.roi_number,'%05.0f'), '.html');
    outT.data_provider_category_automated = class_label_data_provider;
    outT.scientificName_automated = assessed_ID_table.scientificName_automated(label_ind);
    outT.scientificNameID_automated = assessed_ID_table.scientificNameID_automated(label_ind);
    outT.prediction_score_automated_category = class_score';
    outT.biovolume = round(f.Biovolume./(pixel_per_micron^3),3);
    outT.area = round(f.Area/(pixel_per_micron^2),3);
    outT.length = round(f.maxFeretDiameter/pixel_per_micron,3);
    outT.width = round(f.minFeretDiameter/pixel_per_micron,3);
    outT.esd = round((outT.biovolume/4*3/pi).^(1/3)*2,3);
    writetable(outT, outfullfile, 'WriteVariableNames', false, 'FileType', 'text', 'WriteMode', 'append')
    if ispc
            % This is dumb and time consuming but it gets rid of the CRLF and changes to LF only so fcheck doesn't throw warnings
        % maybe some day matlab will add the LineEnding option to writetable
        datalines = readlines(outfullfile);
        fid = fopen(outfullfile, 'w');
        fprintf(fid,'%s\n',datalines{:});
        fclose(fid);
    end
end
end

