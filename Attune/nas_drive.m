function out = nas_drive(action, varargin)
% function out = nas_drive(action, ...)
% The drive letter this computer maps the lab NAS share (\\sosiknas1\Lab_data) to,
% kept as a MATLAB preference (group 'attune', name 'nas_drive') so it carries over
% to later MATLAB sessions. Default Z. Review pages save decisions only when they
% are opened from the mapped drive, so messages and the cruise console show paths
% on that drive.
%
%   L = nas_drive()                 the letter, e.g. 'Z'
%   L = nas_drive('set', L)         store a new letter (A-Z; 'Z:' and 'Z:\' accepted)
%   p = nas_drive('path', f)        f with \\sosiknas1\Lab_data\ replaced by L:\
%   p = nas_drive('path', f, L)     the same with a given letter
%
% Nov 2026, cruise console

if nargin < 1 || isempty(action), action = 'get'; end
switch lower(char(action))
    case 'get'
        out = 'Z';
        if ispref('attune', 'nas_drive')
            L = getpref('attune', 'nas_drive');
            if ischar(L) && isscalar(L) && isletter(L), out = upper(L); end
        end
    case 'set'
        L = upper(regexprep(strtrim(char(varargin{1})), '[:\\/]+$', ''));
        if ~isscalar(L) || ~isletter(L)
            error('nas_drive:letter', 'drive letter must be one letter A-Z, not "%s"', char(varargin{1}))
        end
        setpref('attune', 'nas_drive', L);
        out = L;
    case 'path'
        f = char(varargin{1});
        if numel(varargin) > 1, L = char(varargin{2}); else, L = nas_drive(); end
        unc = '\\sosiknas1\Lab_data\';
        if strncmpi(f, unc, numel(unc)), f = [L ':\' f(numel(unc)+1:end)]; end
        out = f;
    otherwise
        error('nas_drive:action', 'unknown action %s', char(action))
end
end
