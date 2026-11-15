function eb = review_euk_line(F, cls)
% function eb = review_euk_line(F, cls)
% Euk chl-size line [intercept slope] fitted (robustfit) to the v2 Euks of one sample,
% as in plot_class_cytograms; [NaN NaN] with fewer than 50 Euks. The review panel
% "chl minus Euk chl-size line" plots F.BL3 - (eb(1) + eb(2) * F.size).
%
% Nov 2026, classification review prototype

eb = [NaN NaN];
k = cls(:) == 1;
if sum(k) >= 50
    w = warning('off', 'stats:statrobustfit:IterationLimit');
    eb = robustfit(F.size(k), F.BL3(k))';
    warning(w);
end
end
