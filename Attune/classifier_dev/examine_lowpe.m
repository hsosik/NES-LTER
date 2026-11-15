function R = examine_lowpe(tags, outdir)
% function R = examine_lowpe(tags, outdir)
% Is class 3 (low-PE euk, Syn-like PE spectrum) a real group, or Euks carrying /
% coinciding with Syn? For each scored file: class 3 PE relative to one Syn
% (GL2 - Syn median GL2, log10 Syn equivalents), chl relative to Euks of the same
% size, PE-chl correlation within class 3 and within high-PE euks (4), and per-file
% rates (class 3 per Euk vs Syn per event) to test a coincidence/attachment model.
% Writes outdir/lowpe_<cruise>.png and returns per-file table R.
%
% Nov 2026, classifier redesign prototype

addpath(fullfile(fileparts(mfilename('fullpath')), '..'))
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
cfg = classifier_config();
R = table; E = table;
for t = tags
    a = load(fullfile(outdir, ['agreement_' t{1} '.mat']), 'A'); A = a.A;
    for i = 1:height(A)
        cpath = fullfile(root, char(A.cruise(i))); fn = char(A.file(i));
        bs = load_beadstat(fullfile(cpath, 'bead_calibrated'));
        [d, h] = fca_readfcs(fullfile(cpath, 'FCS', fn));
        ref = bead_reference(bs, A.time(i), cfg);
        [new, info] = assign_class_v2(d, h, ref, cfg, A.pro_measured(i));
        F = info.features;
        sy = new == 2; eu = new == 1; c3 = new == 3;
        if sum(sy) < 50 || sum(eu) < 50, continue, end
        synGL2 = median(F.GL2(sy));
        % Euk chl expected for size (robust line through Euks)
        p = robustfit(F.size(eu), F.BL3(eu));
        dchl = F.BL3 - (p(1) + p(2) * F.size);
        % volume analysed: use event counts per class as concentration proxies within a file
        R = [R; table(string(t{1}), string(fn), height(F), sum(sy), sum(eu), sum(c3), sum(new == 4), sum(new == 8), ...
            median(F.GL2(c3)) - synGL2, median(dchl(c3)), median(dchl(eu)), ...
            'VariableNames', {'cruise' 'file' 'n_events' 'n_syn' 'n_euk' 'n_lowpe' 'n_highpe' 'n_synpart' 'lowpe_PE_vs_syn' 'lowpe_dchl' 'euk_dchl'})]; %#ok<AGROW>
        if any(c3)
            E = [E; table(repmat(string(t{1}), sum(c3), 1), F.GL2(c3) - synGL2, dchl(c3), F.BL3(c3), F.size(c3), F.GL2(c3), ...
                'VariableNames', {'cruise' 'PE_vs_syn' 'dchl' 'BL3' 'size' 'GL2'})]; %#ok<AGROW>
        end
        hp = new == 4;
        if sum(hp) > 20
            E = [E; table(repmat(string(t{1}) + "_hp", sum(hp), 1), F.GL2(hp) - synGL2, dchl(hp), F.BL3(hp), F.size(hp), F.GL2(hp), ...
                'VariableNames', {'cruise' 'PE_vs_syn' 'dchl' 'BL3' 'size' 'GL2'})]; %#ok<AGROW>
        end
    end
    C = E(E.cruise == t{1}, :); H = E(E.cruise == string(t{1}) + "_hp", :); r = R(R.cruise == t{1}, :);
    fprintf('\n%s: low-PE euks %d in %d files\n', t{1}, height(C), height(r))
    fprintf('  PE relative to one Syn (log10): p10/50/90 %s\n', mat2str(prctile(C.PE_vs_syn, [10 50 90]), 2))
    fprintf('  chl relative to Euks of same size: low-PE %s; Euks median %s\n', mat2str(prctile(C.dchl, [10 50 90]), 2), mat2str(median(r.euk_dchl), 2))
    fprintf('  corr(PE, chl): low-PE %.2f', corr(C.GL2, C.BL3))
    if height(H) > 20, fprintf(', high-PE %.2f', corr(H.GL2, H.BL3)), end
    fprintf('\n')
    % coincidence / attachment model: low-PE per Euk should scale with Syn per event
    k = r.n_lowpe > 0;
    if sum(k) > 5
        x = log10(r.n_syn(k) ./ r.n_events(k)); y = log10(r.n_lowpe(k) ./ r.n_euk(k));
        b = polyfit(x, y, 1); cc = corr(x, y);
        fprintf('  across files: log(low-PE/Euk) on log(Syn/event): slope %.2f, r %.2f (coincidence predicts slope ~1)\n', b(1), cc)
    end
end
f = figure('Visible', 'off', 'Position', [20 20 1500 450]); tiledlayout(1, 3, 'TileSpacing', 'compact');
nexttile; hold on
for t = tags
    C = E(E.cruise == t{1}, :); histogram(C.PE_vs_syn, -2:0.1:2, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
end
xline(0, 'k--', 'one Syn'); xline(log10(2), 'k:', 'two Syn'); xlabel('low-PE euk GL2 - Syn median GL2 (log10)'); ylabel('fraction'); legend(tags, 'FontSize', 7); grid on
nexttile; hold on
for t = tags
    C = E(E.cruise == t{1}, :); scatter(C.BL3, C.GL2, 4, 'filled', 'MarkerFaceAlpha', .3);
end
xlabel('chl (BL3)'); ylabel('PE (GL2)'); title('low-PE euks'); grid on
nexttile; hold on
for t = tags
    r = R(R.cruise == t{1} & R.n_lowpe > 0, :);
    scatter(log10(r.n_syn ./ r.n_events), log10(r.n_lowpe ./ r.n_euk), 25, 'filled');
end
xlabel('log10 Syn per event'); ylabel('log10 low-PE euks per Euk'); legend(tags, 'FontSize', 7, 'Location', 'northwest'); grid on
exportgraphics(f, fullfile(outdir, 'lowpe_summary.png'), 'Resolution', 90); close(f)
save(fullfile(outdir, 'lowpe_summary.mat'), 'R', 'E')
end
