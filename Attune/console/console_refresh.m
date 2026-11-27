function st = console_refresh(st, cruise, sample)
% function st = console_refresh(st, cruise, sample)
% Cruise console 'refresh' button: rescans one cruise and sample type (or every
% cruise when cruise is "") with cruise_status_all, which updates the cache and its
% browser copy, and returns the console state (st.S, st.made) read back from it.
%
% Nov 2026, cruise console

cruise = string(cruise); sample = string(sample);
if cruise == ""
    [~, S] = cruise_status_all(Root=st.root, Show=false);
else
    if sample == "", sample = ["underway" "preserved"]; end
    [~, S] = cruise_status_all(Cruises=cruise, Samples=sample, Root=st.root, Show=false);
end
root = st.root; if root == "", root = "\\sosiknas1\Lab_data\Attune\cruise_data\"; end
o = load(fullfile(root, 'processing_status', 'status.mat'), 'made');
st.S = S; st.made = o.made;
end
