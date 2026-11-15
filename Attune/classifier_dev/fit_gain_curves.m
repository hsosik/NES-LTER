function G = fit_gain_curves(B, chans)
% function G = fit_gain_curves(B, chans)
% Fit detector response vs voltage from bead runs, assuming a power law
%   log10(signal) = k * log10(hv) + offset(cruise)
% Beads are a fixed target, so within-cruise changes in hv identify k;
% the per-cruise offset absorbs bead lot, laser and alignment differences.
% Fits each bead population separately, then pools those on scale.
%
% Nov 2026, classifier redesign exploration

lo = log10(30); hi = log10(9e5); %stay off the noise floor and saturation
G = table;
for ch = chans
    hv = B.(['hv_' ch{1}]);
    rows = [];
    for b = 1:3
        y = B.(sprintf('bead%d_%s', b, ch{1}));
        ok = isfinite(y) & y > lo & y < hi & isfinite(hv) & hv > 0;
        if sum(ok) < 10, continue, end
        d = table(y(ok), log10(hv(ok)), categorical(B.cruise(ok)), 'VariableNames', {'y' 'lhv' 'cruise'});
        % only cruises with more than one voltage inform the slope
        [g, ~] = findgroups(d.cruise);
        nv = splitapply(@(x) numel(unique(x)), d.lhv, g);
        d = d(nv(g) > 1, :);
        if height(d) < 5, continue, end
        d.cruise = removecats(d.cruise);
        if numel(categories(d.cruise)) > 1
            m = fitlm(d, 'y ~ lhv + cruise');
        else
            m = fitlm(d, 'y ~ lhv');
        end
        k = m.Coefficients{'lhv', {'Estimate' 'SE'}};
        rows = [rows; b k height(d) numel(categories(d.cruise)) numel(unique(d.lhv))]; %#ok<AGROW>
    end
    r = table(string(ch{1}), NaN, NaN, 0, 0, 0, "", 'VariableNames', {'chan' 'k' 'k_se' 'n' 'ncruise' 'nhv' 'detail'});
    if ~isempty(rows)
        w = 1 ./ rows(:, 3).^2; %inverse-variance pooled across bead sizes
        r.k = sum(w .* rows(:, 2)) / sum(w);
        r.k_se = sqrt(1 / sum(w));
        r.n = sum(rows(:, 4));
        r.ncruise = max(rows(:, 5));
        r.nhv = max(rows(:, 6));
        r.detail = strjoin(compose('bead%d k=%.2f±%.2f (n=%d)', rows(:, 1), rows(:, 2), rows(:, 3), rows(:, 4)), '; ');
    end
    G = [G; r]; %#ok<AGROW>
end
end
