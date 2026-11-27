function console_publish(outdir, S, made)
% function console_publish(outdir, S, made)
% Writes the cruise console page for browsers into outdir: index.html (copy of
% console\cruise_console_page.html), status.js (window.STATUS = ..., which a page
% opened from disk can load) and status.json (same content). In a browser the page is
% read-only: overview, details and links; refresh and runs need cruise_console.
%
% Nov 2026, cruise console

if ~isfolder(outdir), mkdir(outdir); end
txt = status_json(S, made, 'browser');
writetxt(fullfile(outdir, 'status.json'), txt);
writetxt(fullfile(outdir, 'status.js'), ['window.STATUS = ' txt ';' newline]);
copyfile(fullfile(fileparts(mfilename('fullpath')), 'cruise_console_page.html'), fullfile(outdir, 'index.html'), 'f');
end


function writetxt(fn, txt)
fid = fopen(fn, 'w', 'n', 'UTF-8');
if fid < 0, error('console_publish:write', 'cannot write %s', fn); end
fwrite(fid, txt, 'char');
fclose(fid);
end
