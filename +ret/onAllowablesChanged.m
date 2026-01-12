function onAllowablesChanged(src, ~)
S = guidata(src);
tag = string(get(src,'Tag'));
val = str2double(get(src,'String'));
if ~isfinite(val), val = NaN; end

% Tag format: allow_y_<field> or allow_u_<field>
parts = split(tag, "_");
if numel(parts) < 3
    guidata(S.fig,S);
    return;
end

which = parts(2); % y or u
field = strjoin(parts(3:end), "_");

if which == "y"
    S.allow.yield.(field) = val;
else
    S.allow.ult.(field) = val;
end

guidata(S.fig,S);
retention_app_V4('update');
end