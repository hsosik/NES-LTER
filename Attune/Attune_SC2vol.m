function [vol_cubic_micron, func_str] = Attune_SC2vol(SCdat_bdnorm,SCpar, prelim_flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%11 sep 2026, Heidi update with new coefficients from July2026 cal
%14 Nov 2026, SSC-H refit (same 12 culture rows and SSC->GL1 switch as
%   size_calibration_fit_jul2026.m, T. weissflogii kept) with log-space modes:
%   Attune culture modes = KDE peak of log10 signal (were linear-bin histogram
%   modes, 0.03-0.19 log low for broad populations) and Coulter modes per log
%   diameter (compile_cc_histograms.m, CCSummary_saved_14Nov2026).
%   95% CI p1 1.009-1.383, p2 0.879-1.286; rmse 0.261 log10 vol (n = 12);
%   ~x4 either way for a single culture. See classifier_dev/pt_transfer/README.md
%   SSC-A not refit (still linear-bin modes).
switch SCpar
    case 'SSC-A'
        p1 = 1.1108; %2019 v2: p1 = 1.2166;
        p2 = 1.0497;  %2019 v2: p2 = 1.1717;
    case 'SSC-H'
        p1 = 1.1962; %11 Sep 2026: 1.1696; 2019 v2: p1 = 1.3008;
        p2 = 1.0826; %11 Sep 2026: 1.1592; 2019 v2: p2 = 1.3240;
    case 'SSC-H_488_2019'
        % cruises before Aug 2019 (one scatter channel, SSC 488/10): March/April
        % 2019 cultures with OD2 on SSC (linear), 12 cultures, log-space Attune
        % and Coulter modes; 95% CI p1 1.23-1.59, p2 1.03-1.40; rmse 0.233.
        % OD2-on-SSC cruises: apply to SSC-H / 1 um bead directly; no-OD2
        % cruises: decompress first (ssc_decompress). Nov 2026,
        % classifier_dev/pt_transfer/sizecal_2019.m
        p1 = 1.4066;
        p2 = 1.2109;
    %case 'SSC-W'  %NOT GOOD...
    %    p1 = 2.151;
    %    p2 = 0.5907;
end
if prelim_flag == 1 % hack to get rough estimates in 2026 before cal
    p1 = 1.1358;
    p2 = 0.8762;
end
vol_cubic_micron = 10.^(p1*log10(SCdat_bdnorm)+p2);
func_str = ['vol_cubic_microns = 10.(log10(SSCmerge./bead_value)*' num2str(p1) ' + ' num2str(p2) ')'];
end