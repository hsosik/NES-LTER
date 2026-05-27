function out = calc_Bofu_NO3(spec, ncal, Pcor_flag)
% 
% PURPOSE:
%   Calculate nitrate concentration for each SUNA sample.
%   This function calculates nitrate concentration from data supplied as a
%   structure retuned from the function parse_SOSIK_SUNA.m as well as
%   calibration data supplied from parseNO3cal.m. See parse_SOSIK_SUNA.m and 
%   parseNO3cal.m for more details on output structure
%
% USAGE:
%	out = calc_glider_NO3(spec, ncal, Pcor_flag)
%
% INPUTS:
%	spec      =  a structure built by parse_SUNA_tsv.m
%   ncal      =  a structure built by parseNO3cal.m with NO3 cal data
%   Pcor_flag = 1 to pressure correct ESW, otherwise don't
%
% OUTPUTS:
%   out =   a structure of nitrate data and a cell array of column titles 

% HISTORY


% *************************************************************************
% TESTING
% data_dir  = 'C:\Users\jplant\Documents\MATLAB\SUNA_GLIDER\SUNA_files\';
% spec_fn   = 'A0001230-SUNA0642.tsv';
% bin_fn    = [regexp(spec_fn,'\w+', 'once','match'),'.BIN'];
% ctd_fn    = 'm6_ctd4suna.csv';
% cal_fn    = 'SNA0642N.cal';
% 
% suna_path = [data_dir,'ASCII\', spec_fn];
% bin_path  = [data_dir,'DAT\', bin_fn];
% ctd_path  = [data_dir, ctd_fn];
% cal_path  = [data_dir, cal_fn];
% 
% ncal      = parseNO3cal(cal_path);
% spec      = parse_SUNA_tsv(suna_path,ctd_path,bin_path);
% Pcor_flag = 1;

% ************************************************************************

% ************************************************************************
% CONSTANTS
% ************************************************************************
fig_flag = 0; % 1 = show nitrate fit
WL_offset = ncal.WL_offset;  % Optical wavelength offset


% ************************************************************************
% CALCULATE RAW ABSORBANCE 
% CALIBRATION DATA MEASURED FROM ~190 to ~400 nm IN THE LAB. FLOAT CAN'T
% STORE ALL WAVELENGTHS SPECTROPHOTOMETER MEASURES SO ONLY A SUBSET IS
% RETURNED. 
% ************************************************************************
d = spec; clear spec
d.P = d.STP(:,3);
d.T = d.STP(:,2);
d.S = d.STP(:,1);
rmfield(d,'STP');

% MAKE SOME CALIBRATION ADJUSTMENTS BASED ON DATA FLAGS IN CAL STRUCTURE
% 0 based count - pixel registration to wavelength is off by one for some
% instruments. if 217 was pixel 30 it is now pix 31. DOUBLE CHECK THIS!!!
if ncal.pixel_base == 0 
    d.spectra_pix_range = d.spectra_pix_range +1;
    d.pix_fit_win = d.pix_fit_win + 1;
    disp('Pixel registration offset by +1')
end
% CHOOSE REGULAR or SEAWATER DARK CURRENT FOR CALCULATION
DC = d.DC; % default
if ncal.DC_flag == 0
    DC = d.SWDC;
end


% USE 1st >= 217 and LAST <= 240 unless cal file has updated fit window -
% not doing this in this file

% ADJUST FIT WINDOW IF IT HAS BEEN UPDATED IN THE CAL STRUCTURE
if ~isempty(ncal.min_fit_WL)
    d.pix_fit_win(1) = find(ncal.WL >= ncal.min_fit_WL,1,'first');
    d.WL_fit_win(1)  = ncal.WL(d.pix_fit_win(1));
% else
%     d.pix_fit_win(1) = find(ncal.WL >= 217,1,'first'); % MBARI DEFAULTS
end

if ~isempty(ncal.max_fit_WL)
    d.pix_fit_win(2)  = find(ncal.WL <= ncal.max_fit_WL,1,'last');
    d.WL_fit_win(2)   = ncal.WL(d.pix_fit_win(2));
% else
%     d.pix_fit_win(2)  = find(ncal.WL <= 240 ,1,'last'); % MBARI DEFAULTS
end

% CHECK DATA RANGE - IF WLs of returned spectra don't overlap fit window
% the fit window needs to shrink! THIS IS NOT A GOOD SITUATION SO CHECK THE
% DATA!!!
if d.spectra_WL_range(1) > d.WL_fit_win(1)
    fprintf(['Min WL of returned spectra (%0.2f) is greater than Min WL of ', ...
        'fit window (%0.2f).\n'], d.spectra_WL_range(1), d.WL_fit_win(1));
        d.WL_fit_win(1) = d.spectra_WL_range(1);
    fprintf(['Fit window adjusted [%0.2f %0.2f]& NO3 estimate will be ', ...
        'compromised!\n'], d.WL_fit_win);  
end

if d.spectra_WL_range(2) < d.WL_fit_win(2)
    fprintf(['Max WL of returned spectra (%0.2f) is less than Max WL of ', ...
        'fit window (%0.2f).\nFit window adjusted & NO3 estimate will be ', ...
        'compromised!\n'], d.spectra_WL_range(2), d.WL_fit_win(2));
        d.WL_fit_win(2) = d.spectra_WL_range(2);
end

% GET PIXEL FIT WINDOW FOR SUNA FLOATS BASED ON WL FIT RANGE
% SUNA FLOATS DON'T HAVE THIS AS A HEADER LINE
if isnan(sum(d.pix_fit_win)) && isfield(d,'WL_fit_win')
    pfit_low = find(ncal.WL >= d.WL_fit_win(1),1,'first');
    pfit_hi  = find(ncal.WL <= d.WL_fit_win(2),1,'last');
    d.pix_fit_win = [pfit_low pfit_hi];
    clear pfit_low pfit_hi
end


% GET SAMPLE PIXEL RANGE INDEXES 
ind1 = find(ncal.WL >= d.spectra_WL_range(1),1,'first');
ind2 = find(ncal.WL <= d.spectra_WL_range(2),1,'last');
d.spectra_pix_range = [ind1,ind2];
saved_pixels = d.spectra_pix_range(1): d.spectra_pix_range(2); % index
clear ind1 ind2


% SIZE OF SAMPLE INTENSITY SPECTRA MATRIX(samples X wavelengths)
[rows,cols]  = size(d.UV_INTEN);

% SUBSET CALIBRATION DATA OVER PIXEL RANGE IDENTIFIED FROM SAMPLE LINE
% THEN TRANSPOSE TO ROWS & REPEAT ROWS SO SAME SIZE AS UV_INTEN
REF  = ones(rows,1) * ncal.Ref(saved_pixels)';  % Reference intensities
WL   = ones(rows,1) * ncal.WL(saved_pixels)';   % Wavelengths
ESW  = ones(rows,1) * ncal.ESW(saved_pixels)';  % Extinction coefs SW

ENO3 = ncal.ENO3(saved_pixels)'; % Extinction coefs NO3 (array)

if isfield(ncal,'EHS') % NOT USED YET BUT CARRY ALONG
    EHS  = ncal.EHS(saved_pixels)';  % Extinction coefs HS  (array)
end

% SUBTRACT DARK CURRENT INTENSITIES
%UV_INTEN_SW = d.UV_INTEN - DC * ones(1,cols); 
UV_INTEN_SW = d.UV_INTEN - DC; % Seabird *.tsv file DC for each sample WL 05/02/2018 jp
% Get UV reference intensity for pixel range from cal file

% CHECK FOR UV_INTEN_SW <= 0 (DC >= UV_INTEN_SW) & SET TO NaN
tINTEN_SW = UV_INTEN_SW <= 0;
if sum(tINTEN_SW(:)) > 0
    UV_INTEN_SW(tINTEN_SW) = NaN;
    disp(['UV Intensities <= DC found for ', ...
        '. Setting these low intensities to NaN']);
    %clear tINTEN_SW
end
    
% ABSORBANCE SPECTRA FOR ALL SAMPLES IN PROFILE
ABS_SW = -log10(UV_INTEN_SW ./ REF);

% ************************************************************************
% CALCULATE TEMPERATURE CORRECTED ABSORBANCES
% ************************************************************************
A = 1.1500276; % Coefficients from Sakamoto et al., 2009
B = 0.02840;
C = -0.3101349;
D = 0.001222;

% Build matrices for calculations all size of spectra - avoid some loops 
% calculate corrected absorbances for all samples all at once
% ctd_temp  = d.T * ones(1,cols);           % measured sw temperature, C
% ctd_sal   = d.S * ones(1,cols);           % measured sw salinity
% cal_temp  = ones(rows,cols) * ncal.CalTemp;  % temp isus calibrated at, C

ctd_temp  = d.T * ones(1,cols);           % measured sw temperature, C
ctd_sal   = d.S * ones(1,cols);           % measured sw salinity
cal_temp  = ones(rows,cols) * ncal.CalTemp;  % temp isus calibrated at, C

%disp('Using cal temp from cal file not from message file')

% Opt_wl_UV = ones(rows,1) * (isus_cal(:,1))';  % matrix of wavelengths
% E_swa_tcal = ones(rows,1) * (isus_cal(:,2))'; % sw extinction coefs, cal t

% % Refer to Computing DAC Nitrate v0.2
% % Sakamoto temperature correction "f function" at ctd temperature
% F_sw = (A + B .* ctd_temp) .* exp((C + D .* ctd_temp) .* ...     % Eqn. 3
%        (WL - WL_offset));
% 
% % Sakamoto temperature correction "f function" at calibration temperature
% % WL_offset
% F_cal = (A + B .* cal_temp) .* exp((C + D .* cal_temp) .* ...    % Eqn. 3
%         (WL - WL_offset));
% 
% % Calculate in situ sea water extinction coefficients (mostly due to Br)
% ESW_in_situ = ESW .* F_sw ./ F_cal;      % Eqn. 2


% UPDATED ESW TEMPERATURE CORRECTION ALGORITHM TO dLN[ESW]/dT 03/30/2022 JP
% (03/29/2022 JP)EQUITECH PROBES ONLY with 3 anamolously low float experiments
% removed 05/20/2016, 11/21/2016, 10/09/2017 (P CHAMBER EXPERIMENTS?)
Tcorr_coef  = [1.27353e-07 -7.56395e-06 2.91898e-05 1.67660e-03 1.46380e-02];
Tcorr       = polyval(Tcorr_coef,(WL - WL_offset)) .* (ctd_temp - cal_temp);
ESW_in_situ = ESW .* exp(Tcorr);

% Add pressure corrrection if flag is on
if Pcor_flag == 1
    pres_term   = (1- d.P ./ 1000 .* ncal.pres_coef) * ones(1,cols);
    ESW_in_situ = ESW_in_situ .* pres_term;
end

% SW (bromide) absorbances
ABS_Br_tcor      = ESW_in_situ .* ctd_sal; % in situ bromide absorbances

% Calculate predicted absorbance of nitrate + baseline
% The 2nd term here is removing the temperature corrected absorbance
% of seawater from the spectra
ABS_cor = ABS_SW - ABS_Br_tcor  ;    % Eqn. 4

clear ctd_temp cal_temp ctd_sal REF ESW pres_term

% ************************************************************************
% CALCULATE THE NITRATE CONCENTRATION & THE SLOPE AND INTERCEPT OF THE
% BASELINE ABSORBANCE
% This calulation occurs over the "fit window" which is usually a subset
% of the "spectra window" which is a subset of the calibration spectra
% window. The calibration window is the full wavelength range of the
% spectraphotometer
% ************************************************************************

% SUBSET TO FIT WINDOW
t_fit     = saved_pixels >= d.pix_fit_win(1) & ...
            saved_pixels <= d.pix_fit_win(2); % fit wind pixel indicies
      
Fit_WL    = WL(1,t_fit)';
Fit_ENO3  = ENO3(t_fit)';
Ones      = ones(size(Fit_ENO3));

M     = [Fit_ENO3 Ones/100 Fit_WL/1000]; % /1000 numerical stability??
M_INV = pinv(M);
% /100 & /1000 may be in here for numerical stability. Need to /100 and
% /1000 at the end to get back to actual baseline and intercept
colsM = size(M,2);
% Columns may need to be added for more solutions
% # samples x (# fit parameters + QC values)
NO3 = ones(rows, colsM +3)* NaN;

for i = 1:rows
    samp_ABS_cor = (ABS_cor(i,t_fit))'; % NO3 + BL
    
    % NO3 = [ baseline intercept, and slope]
    NO3(i,1:3) = M_INV * samp_ABS_cor;
    NO3(i,2) = NO3(i,2)/100; % baseline intercept
    NO3(i,3) = NO3(i,3)/1000; % baseline slope
    
    ABS_BL  = WL(1,:) .* NO3(i,3) + NO3(i,2); % Baseline absorbance
    ABS_NO3 = ABS_cor(i,:) - ABS_BL; % Calculated nitrate absorbance
    ABS_NO3_EXP = ENO3 .* NO3(i,1);
    FIT_DIF = ABS_cor(i,:) - ABS_BL -ABS_NO3_EXP;
    RMS_ERROR = sqrt((sum(FIT_DIF(t_fit).^2)./sum(t_fit)));
    ind_240  = find(Fit_WL <= 240,1,'last');
    ABS_240  = [Fit_WL(ind_240), samp_ABS_cor(ind_240)]; % WL ~240 & ABS
    
    NO3(i,colsM+1:colsM+3) = [RMS_ERROR ABS_240]; % RMS WL~240 ABS~240
    % ********************************************************************
    % TEST -- PLOT SPECTRA
    % ********************************************************************
    %     if RMS_ERROR > 0.003 || ABS_240(2) > 0.8 % ONLY PLOT IF NOT GOOD
    if fig_flag == 1 && ~isnan(NO3(i,1))
        clf
        figure(101)
        set(gcf, 'Units', 'normalized');
        set(gcf,'Position',[0.17 0.23 0.52 0.55])
        a1 = axes('Position', [0.13 0.12 0.775 0.59]);
        plot(WL(1,:), ABS_SW(i,:), 'k-','LineWidth', 2) % sw spectrum
        hold on
        %x_lim = [min(WL(1,:)) max(WL(1,:))];
        x_lim = [min(Fit_WL(:,1))-5 max(Fit_WL(:,1))+5];
        xlim(x_lim)
        y_lim = ylim;
        set(gca, 'Linewidth',3)
        xlabel('Wavelength, nm', 'FontSize', 14)
        ylabel('Absorbance','FontSize', 14)
        set(gca, 'FontSize', 14)
        
        plot(WL(1,:), ABS_BL, 'g-','LineWidth', 2) % BL
        plot(WL(1,:), ABS_Br_tcor(i,:)+ABS_BL, 'b-','LineWidth', 2) % Br spectrum (+ BL, % Bofu added)
        plot(WL(1,:), ABS_cor(i,:), 'ro','LineWidth', 1) % NO3+BL observed
        plot(WL(1,:), ABS_NO3_EXP +ABS_BL, 'r-','LineWidth', 2) % NO3+BL exp
        
        y_lim = ylim;
        ylim(y_lim)
        
        fit_patch = fill([Fit_WL(1),Fit_WL(1), Fit_WL(end),Fit_WL(end)], ...
            [y_lim(1), y_lim(2), y_lim(2), y_lim(1)],'y');
        set(fit_patch,'FaceAlpha',0.20, 'LineStyle','--','EdgeColor', ...
            [160 160 160]/256, 'LineWidth',2)
        
        %     s1 = {['FIT WINDOW (',num2str(Fit_WL(1),'%4.1f'),' - ', ...
        %           num2str(Fit_WL(end),'%4.1f'),')'],['NO3 = ', num2str(NO3(i,1),...
        %           '%2.1f'), ' �M'], ['RMS error = ', num2str(NO3(i,4),'%2.3e')]};
        
        s1 = {d.file,['Depth:', num2str(d.P(i),'%0.1f')],['FIT WINDOW (',num2str(Fit_WL(1), ...
            '%4.1f'),' - ', num2str(Fit_WL(end),'%4.1f'),')'], ...
            ['NO3 = ', num2str(NO3(i,1),'%2.1f'), ' \muM'], ...
            ['RMS error = ', num2str(NO3(i,4),'%2.3e')]};
        
%         text (Fit_WL(1) + (Fit_WL(end) - Fit_WL(1))/2 , ...
%             y_lim(1) + (y_lim(2) - y_lim(1))* 0.9,...
%             s1,'FontSize', 16, 'HorizontalAlignment','center', ...
%             'VerticalAlignment', 'Top')
        text (Fit_WL(1,1) + (Fit_WL(end,1) - Fit_WL(1,1))* 0.4, ...
            y_lim(1) + (y_lim(2) - y_lim(1))* 0.95,...
            s1,'FontSize', 16, 'HorizontalAlignment','center', ...
            'VerticalAlignment', 'Top')
    
        
        set(gca, 'box', 'on')
        
        
        set(get(get(fit_patch,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Don't show patch
        legend('Sample','Baseline', 'Br+BL', 'NO3+BL obs', 'NO3+BL fit', ...
            'LOcation','NorthEast')
        legend('boxoff')
        
        hold off
        
        % PLOT RESIDUALS
        a2 = axes('Position', [0.13 0.75 0.775 0.20]);
        FIT_DIF_E3 = FIT_DIF*1000;
        max_y = max(abs(FIT_DIF_E3));
        a2_ylim = [-max_y - 0.1*max_y, max_y + 0.1*max_y];
        ylim(a2_ylim)
        xlim(x_lim)
        ylabel({'Obs - Fit','x 10^{3}'},'FontSize', 14)
        set(a2,'Xticklabel','','box', 'on', 'FontSize', 14, 'Linewidth',3);
        hold on
        
        dif_patch = fill([Fit_WL(1),Fit_WL(1), Fit_WL(end),Fit_WL(end)], ...
            [a2_ylim(1), a2_ylim(2), a2_ylim(2), a2_ylim(1)],'y');
        set(dif_patch,'FaceAlpha',0.20, 'LineStyle','--','EdgeColor', ...
            [160 160 160]/256, 'LineWidth',2)
        
        plot(x_lim,[0 0],'--','color','k' ,'LineWidth', 2)
        plot(WL(1,:), FIT_DIF_E3, 'ko','LineWidth', 1, ...
            'MarkerFaceColor','r') % Resid * 1000
        
        hold off
        pause
    end
end

% ************************************************************************
% Create an output matrix and plot the data
% [SDN, DarkCur,Pres,Temp,Sal,NO3,BL_B,BL_M,RMS ERROR,WL~240,ABS~240]
out.hdr = {'SDN' 'AVG DC' 'Pres' 'Temp' 'Sal' 'SBE_NO3, uM/L', ...
           'NO3, uM/L' 'BL_B' 'BL_M' 'RMS ERROR' 'WL~240','ABS~240'};

avg_DC = mean(DC,2,'omitnan');       
NO3 = [d.SDN, avg_DC, d.P, d.T, d.S, d.NO3, NO3];
tnan = isnan(NO3(:,6));
out.data = NO3(~tnan,:); % remove all NaN nitrate lines
out.info.WL_fit_window = d.WL_fit_win;
out.info.spectra_WL_range = d.spectra_WL_range;
out.info.WL_offset = WL_offset;









