function run_preservation_shift_all(cruises, redo)
% function run_preservation_shift_all(cruises, redo)
% preservation_shift2 on every cruise with hand-gated preserved runs
% (preserved/outputs/Gated_Table.mat) and processed underway files
% (bead_calibrated/FCSfileinfo.mat), Aug 2019 layout onward. Skips cruises that
% already have presshift2_<cruise>.mat unless redo. Nov 2026
if nargin < 1 || isempty(cruises)
    cruises = {'20190820_EN644' '20191005_AR39B' '20200201_EN649' '20200311_AR43' '20200725_EN655' ...
        '20210203_EN661' '20210716_EN668' '20211108_AR61B' '20220216_AT46' '20220729_EN687' '20220806_EN688' ...
        '20230111_EN695' '20230807_EN706' '20240209_EN712' '20240503_EN715' '20240906_EN720' '20250124_EN727' ...
        '20250814_AR92' '20260113_AR99'};
end
if nargin < 2, redo = false; end
here = fileparts(mfilename('fullpath'));
od = fullfile(fileparts(here), 'output', 'discrete');
for c = string(cruises)
    if ~redo && exist(fullfile(od, "presshift2_" + extractAfter(c, "_") + ".mat"), 'file'), continue, end
    try
        S = preservation_shift2(char(c));
        fprintf('%s: %d rows\n', c, height(S))
    catch err
        fprintf('%s failed: %s\n', c, err.message)
        disp(err.stack(1))
    end
end
end
