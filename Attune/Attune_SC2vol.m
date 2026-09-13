function [vol_cubic_micron, func_str] = Attune_SC2vol(SCdat_bdnorm,SCpar, prelim_flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%11 sep 2026, Heidi update with new coefficients from July2026 cal
switch SCpar
    case 'SSC-A'
        p1 = 1.1108; %2019 v2: p1 = 1.2166;
        p2 = 1.0497;  %2019 v2: p2 = 1.1717;
    case 'SSC-H'
        p1 = 1.1696; %2019 v2: p1 = 1.3008;
        p2 = 1.1592; %2019 v2: p2 = 1.3240;
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