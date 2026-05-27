function spec = parse_Bofu_SUNA(suna_file)
%
% PURPOSE: 
%   This function parses a SUNA file
%
% USAGE:
%	spec =  parse_Bofu_SUNA(data file name)
%
% INPUTS:
%	suna_file   = Path to suna text file
%   ctd_file    = Path to ctd text file [SDN T S P]
%   time_offset = offset between clocks (suna time - ctd time) in seconds
%                 difference will get added to CTD time
%
% OUTPUTS:
%   spec =   a structure of data and calibration info
%       spec.SDN        = sample time, Matlab sdn        (array)
%      	spec.DC         = Measures DC intensity, counts    (matrix)
%     	spec.UV_INTEN   = Measured UV intensities, counts        (matrix)
%      	spec.SWDC       = sea water DC intensity, counts (scalar)
%                         NaN if doesn't exist (earlier floats)
%       spec.STP        = Salinity, temperature, Pressure of sample
%       spec.NO3        = SBE calculated nitrate
%       spec.WL         = wavelengths of returned sample spectra
%
%      	spec.spectra_pix_range = pixel registration for the wavelengths in
%                                the calibration file
%       spec.WL_Fit_Range      = wavelength bounds of fit window - used in
%                                SUNA floats to determine spectra_pix_range
%    	spec.pix_fit_win       = pixel range for Multiple linear regression.
%                                Used to subset sample spectra to fit window
%
%     	spec.CalTemp    = calibration temperature this is the temperature
%                         the instrument was calibrated at in the lab.
%    	spec.CalDateStr = lab calibration date string
%     	spec.CalDate    = serial date number for the lab cal date

% TESTING
%suna_file = 'C:\Users\jplant\Documents\MATLAB\Bofu_Zheng\A0000477-SUNA1227.csv';


% ************************************************************************
%PREDIM OUTPUT IF NO FILES EXIST
spec.DC         = [];
spec.UV_INTEN   = [];
spec.SWDC       = [];
spec.STP        = [];
spec.NO3        = []; % SBE Calculated nitrate, uM

spec.CalTemp           = NaN;
spec.CalDate           = NaN;
spec.CalDateStr        = '';
spec.WL_fit_win = [217 240];
spec.spectra_pix_range = NaN; % for no float data all pixels returned
spec.spectra_WL_range = NaN; % for no float data all pixels returned
spec.path = suna_file;
spec.pix_fit_win =[NaN NaN];

% sep_ct    = regexp(spec.path,'\\');
% spec.file = spec.path(max(sep_ct)+1:end);


% ************************************************************************
% PARSE SPEC HEADER LINE AND GET FIT WAVELENGTHS FROM HEADER LINE
% ************************************************************************
fid   = fopen(suna_file);
tline = ' ';
while ischar(tline)
    if regexp(tline,'^Index','once')  % stop at hdr line
        t_tab   = regexp(tline,'\t');
        t_comma = regexp(tline,',');
        if size(t_tab, 2) > 10
            delim   = '\t';
        elseif size(t_comma, 2) > 10
            delim   = ',';
        else
            delim   = ' ';
        end
        spec_hdr = regexp(tline,delim,'split'); % break up into col names
        WL  = regexp(tline,'(?<=UV\()\d+\.\d+','match'); % get fit WL's
    elseif regexp(tline,'^#','once')
        break
    end
    tline = fgetl(fid);
end

% BUILD FORMAT STRING FOR DATA PARSING AND PARSE DATA
rhdr = size(spec_hdr,2);
df = '';
for i = 1:rhdr
    if regexp(spec_hdr{i},'DATE|TIME','once')
        df = [df,'%s'];
    else
        df = [df,'%f'];
    end
end
df = [df,'\r\n'];

% CHANGE DELIM IN TEXTSCAN CALL TO MATCH DATA FORMAT
d = textscan(fid, df,'Delimiter', delim); 
fclose(fid);

% Does file have data
if isempty(d)
    disp(['File exists but no data returned for : ',suna_file]);
    return
end

% MERGE DATE STRINGS AND CONVERT TO MATLAB SDN AND COMBINE
d_tmp = strcat(d{1,2},regexprep(d{1,3}, '(\d+:\d+)',' $1'));
sdn   = datenum(d_tmp,'dd mmm yyyy HH:MM:SS');

% CONDENSE REMAINING DATA
iIND  = find(strcmp('Index', spec_hdr)        == 1);
iN    = find(strcmp('NITRATE_UM', spec_hdr)        == 1);
iDAVG = find(strcmp('DARK_AVG', spec_hdr) == 1,1,'last');
iUV1  = find(strncmp('UV(', spec_hdr, 3)      == 1,1,'first'); % start of intensities
iUV2  = find(strncmp('UV(', spec_hdr, 3)      == 1,1,'last'); % end of intensities
iS = find(strcmp('CTD_SAL', spec_hdr)      == 1); % S of STP
iP = find(strcmp('CTD_DEPTH', spec_hdr)    == 1); % P of STP

keep_index = [iN, iDAVG, iUV1:iUV2, iS:iP];
spec_data = [cell2mat(d(:,iIND)), sdn, cell2mat(d(:,keep_index))];
spec_hdr = [spec_hdr(iIND),'SDN', spec_hdr(keep_index)];

% spec_data = [d{1,1},sdn,d{1,3},d{1,5}];
% 
% spec_hdr = spec_hdr([1,2,4:47,49:52]); % Remove 'Time', 'CTD_TIME' 'DATETAG' and 'TIMETAG2' from hdr
clear d d_tmp df i rhdr sdn tline fid_spec ans

% REDO INDICES
iN    = find(strcmp('NITRATE_UM', spec_hdr)        == 1);
iDAVG = find(strcmp('DARK_AVG', spec_hdr) == 1,1,'last');
iUV1  = find(strncmp('UV(', spec_hdr, 3)      == 1,1,'first'); % start of intensities
iUV2  = find(strncmp('UV(', spec_hdr, 3)      == 1,1,'last'); % end of intensities
iS = find(strcmp('CTD_SAL', spec_hdr)      == 1); % S of STP
iP = find(strcmp('CTD_DEPTH', spec_hdr)    == 1); % P of STP
iSDN  = find(strncmp('SDN', spec_hdr, 5)    == 1);

% ************************************************************************
% ASSIGN DATA TO STRUCTURE
% USING TIMP STAMP FROM SAMPLE SCAN NOT DARK SCANS
% ************************************************************************
% returned sample spectra
spec.spectra_WL_range = [str2double(WL{1}) str2double(WL{end})];

% Take average of LIGHT SCANS between DARK SCANS and discard consecutive
% DARK SCANS
% Dark scans appear to have ABS_254, ABS_350, and DARK_AVG set = 0

isDRK = spec_data(:,iDAVG) == 0;
D_ind = find(isDRK == 1);
% CHECK FOR CONSECUTIVE DARK SCANS AND REMOVE ONE OF THEM
tDUP = diff(D_ind) == 0;
if sum(tDUP) > 0
    disp(['Removing consecutive dark scans (',num2str(sum(tDUP)),')']);
    spec_data(D_ind(tDUP),:) =[];
    D_ind(tDUP) = [];
end
rDRK = size(D_ind,1);


% DETERMINE SAMPLE DATA SIZE & predim outputs
nSAMP = size(spec_data,1) - size(D_ind,1);

DC       = ones(nSAMP, iUV2-iUV1+1)*NaN;
UV_INTEN = spec.DC;
SDN      = ones(nSAMP, 1)*NaN; % Sample time stamp
% STP      = ones(nSAMP, 3)*NaN; % STP
NO3      = ones(nSAMP, 1)*NaN; % NO3 conc already in file

% BUILD UP DC MATRIX TO MATCH SAMPLE MATRIX
line_ct = 1;
for dct = 1:rDRK
    D1 = D_ind(dct); % this is DC for following samples
    
    if dct < rDRK % dct+1 exists
        D2 = D_ind(dct+1)-1; % last sample in block
    elseif D1 == size(spec_data,1) % data end with a DC scan - all done here
        break
    else
        D2 = size(spec_data,1); % data ends with sample scans
    end
    
    tmp  = spec_data(D1+1:D2, iUV1:iUV2);
    rtmp = size(tmp,1);
    
    UV_INTEN(line_ct:line_ct+rtmp-1,:) = tmp;
    DC(line_ct:line_ct+rtmp-1,:)  = ones(rtmp,1)*spec_data(D1, iUV1:iUV2);
    SDN(line_ct:line_ct+rtmp-1,:) = spec_data(D1+1:D2, iSDN);
    NO3(line_ct:line_ct+rtmp-1,:) = spec_data(D1+1:D2, iN);
%     STP(line_ct:line_ct+rtmp-1,:) = spec_data(D1+1:D2, [iS,iT,iP]);
    
    line_ct = line_ct+rtmp;
end

spec.DC       = DC;
spec.UV_INTEN = UV_INTEN;
spec.SDN      = SDN;
% spec.STP      = STP;
spec.NO3      = NO3;



