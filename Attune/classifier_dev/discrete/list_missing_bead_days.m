function M = list_missing_bead_days()
% function M = list_missing_bead_days()
% Days on which preserved samples were run (Gated_Table date_processed, all
% cruises) without a QC'd lab bead run of each population on the same day.
% Lists the lab bead runs that day (if any) and why they failed QC, and the
% nearest QC'd run. Saves output/discrete/missing_bead_days.csv
%
% Nov 2026, discrete-sample classifier development

here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
root = '\\sosiknas1\lab_data\Attune\cruise_data\';
lab = load(fullfile(od, 'lab_beadstat.mat')).beadstat;
lday = floor(datenum(lab.time));
P = {'05' '1' '6'};
d = dir(fullfile(root, '*_*', 'preserved', 'outputs', 'Gated_Table.mat'));
M = table;
for i = 1:numel(d)
    try, G = load(fullfile(d(i).folder, d(i).name)).gated_table; catch, continue, end
    G.Properties.VariableNames = lower(G.Properties.VariableNames);
    if ~ismember('date_processed', G.Properties.VariableNames), continue, end
    dp = G.date_processed;
    if ~isdatetime(dp), dp = datetime(string(dp), 'InputFormat', 'dd-MMM-yyyy'); end
    day = floor(datenum(dp));
    cr = string(extractBetween(d(i).folder, root, filesep));
    for u = unique(day(isfinite(day)))'
        missing = {};
        for p = P
            if ~any(lday == u & lab.(['QC' p{1}])), missing{end + 1} = strrep(p{1}, '05', '0.5'); end %#ok<AGROW>
        end
        if isempty(missing), continue, end
        that = lab(lday == u, :);
        why = strjoin(unique([that.why05; that.why1; that.why6]), '; ');
        q = find(lab.QC1);
        [dmin, j] = min(abs(lday(q) - u));
        M = [M; table(cr, datetime(u, 'ConvertFrom', 'datenum'), sum(day == u), string(strjoin(missing, ' ')), height(that), ...
            string(why), dmin, lab.time(q(j)), 'VariableNames', ...
            {'cruise' 'run_day' 'nruns' 'missing_pops_um' 'lab_runs_that_day' 'qc_failures' 'days_to_nearest_QC1um' 'nearest_QC1um_run'})]; %#ok<AGROW>
    end
end
writetable(M, fullfile(od, 'missing_bead_days.csv'))
disp(M)
end
