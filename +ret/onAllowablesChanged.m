function onAllowablesChanged(src, ~)
%RET.ONALLOWABLESCHANGED Update allowables in state and refresh UI live

S = guidata(src);

tag = string(get(src,'Tag'));
raw = string(get(src,'String'));
val = str2double(raw);

% Accept blank as NaN (treat as N/A)
if strlength(strtrim(raw)) == 0
    val = NaN;
end
if ~isfinite(val)
    val = NaN;
end

% Tag format:
%   - allow_y_<field>  (yield column)
%   - allow_u_<field>  (ultimate column)
%   - allow_c_<field>  (single casing allowable; applies to both yield + ultimate)
parts = split(tag, "_");
if numel(parts) < 3
    guidata(S.fig, S);
    return;
end

which = parts(2); % y, u, or c
field = strjoin(parts(3:end), "_");

if which == "c"
    % Single casing allowables apply to both yield + ultimate
    S.allow.yield.(field) = val;
    S.allow.ult.(field)   = val;
elseif which == "y"
    S.allow.yield.(field) = val;
else
    S.allow.ult.(field) = val;
end

% Refresh all derived results (including margins)
S = ret.updateAll(S);

guidata(S.fig, S);

end
