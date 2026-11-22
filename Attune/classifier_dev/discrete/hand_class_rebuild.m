function [c, added] = hand_class_rebuild(gate_names, gate_assign, extended)
% function [c, added] = hand_class_rebuild(gate_names, gate_assign, extended)
% Hand classes for one preserved run from its Gated_Table gate memberships, with the
% precedence of Process_Preserved_AllSteps_2026 step 4 (codes 1 Euk, 2 Syn,
% 3 bacteria, 4 Pro, 5 low-PE euk, 6 high-PE euk; earlier classes keep their events,
% Pro overwrites bacteria). That step maps gates by exact name only, so workspaces
% whose Syn gate has another name get no hand Syn (HB1907 / EN688 Pro runs: SYN_all;
% MVCO MV418N02, MV419N01: Syn_plus; scan, 22 Nov 2026). extended (default true) uses
% SYN_all, Syn_all, Syn_plus, PEB_SYN, PUB_SYN (first found) as Syn, and nonPE_Euk as
% Euk, but only in runs with no regular Syn / Euk gate: where both exist the others
% are subsets or building blocks (EN617 Syn = Syn_sub1 and Syn_sub2, which cut out a
% faint group). false reproduces the stored classes exactly (5,436 of 5,436 runs).
% Ambiguous names are left out: Pro_potential, Pro_all, Pro_sub1,
% large_bacteria(_plus), huge_SYN, weird_SYN, No_SYN_Euk_Pro (an exclusion gate).
% gate_names  - cell of gate names; gate_assign - gates x events membership
% added       - names from the extended lists found in this run
%
% Nov 2026, discrete-sample classifier development

if nargin < 3, extended = true; end
gate_names = cellstr(string(gate_names));
if size(gate_assign, 1) ~= numel(gate_names), gate_assign = gate_assign'; end
gate_assign = logical(gate_assign);
c = zeros(1, size(gate_assign, 2));
added = {};
L = struct();
L.euk = {'euk' 'Euk'};             L.euk_x = {'nonPE_Euk'};
L.lpe = {'LowP_Euk' 'LowPE_Euk'};  L.hpe = {'HighP_Euk' 'HighPE_Euk'};
L.syn = {'Syn' 'syn' 'SYN'};       L.syn_x = {'SYN_all' 'Syn_all' 'Syn_plus' 'PEB_SYN' 'PUB_SYN'};
L.bac = {'Bacteria' 'bacteria' 'BacPro'};
L.pro = {'pro' 'Pro' 'Euk_sm' 'Proc'};
    function fill(names, code, overwrite)
        for n = names
            j = strcmp(gate_names, n{1});
            if ~any(j), continue, end
            m = any(gate_assign(j, :), 1);
            if overwrite, c(m) = code; else, c(c == 0 & m) = code; end
        end
    end
    function fillx(names, code, regular)
        % extra names only where the run has none of the regular gates for this class
        % (where both exist they are subsets or building blocks of the regular gate,
        % e.g. PEB_SYN / PUB_SYN inside Syn, Syn = Syn_sub1 and Syn_sub2 in EN617)
        if ~extended || any(ismember(regular, gate_names)), return, end
        k = names(ismember(names, gate_names));
        if isempty(k), return, end
        k = k(1); %the first listed (e.g. SYN_all rather than its PEB/PUB parts)
        added = [added k]; %#ok<AGROW>
        fill(k, code, false);
    end
fill(L.euk, 1, false); fillx(L.euk_x, 1, L.euk);
fill(L.lpe, 5, false); fill(L.hpe, 6, false);
fill(L.syn, 2, false); fillx(L.syn_x, 2, L.syn);
fill(L.bac, 3, false);
fill(L.pro, 4, true);
end
