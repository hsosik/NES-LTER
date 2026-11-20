% score an automatic classifier against the hand gates on preserved runs
% set tag (and optionally cfg, classify) before running; Nov 2026
if ~exist('tag', 'var'), tag = 'test'; end
if ~exist('cfg', 'var'), cfg = []; end
if ~exist('classify', 'var'), classify = []; end
if ~exist('nper', 'var'), nper = 6; end
cruises = {'20200201_EN649' '20200725_EN655' '20210203_EN661' '20210716_EN668' '20211108_AR61B' '20220729_EN687' ...
    '20230111_EN695' '20230807_EN706' '20240503_EN715' '20240906_EN720' '20250124_EN727' '20250814_AR92'};
A = score_discrete(cruises, nper, cfg, tag, classify);
G = groupsummary(A, {'runtype'}, {'median' 'min'}, 'agreement');
disp(G)
S = groupsummary(A, {'cruise' 'runtype'}, 'sum', {'hand_1' 'auto_1' 'hand_2' 'auto_2' 'hand_7' 'auto_7' 'hand_34' 'auto_34' 'auto_8'});
disp(S)
