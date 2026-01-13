function S = readConstraintsFromUI(S)

% Ensure structs exist
if ~isfield(S,'allow') || ~isstruct(S.allow), S.allow = struct(); end
if ~isfield(S.allow,'yield') || ~isstruct(S.allow.yield), S.allow.yield = struct(); end
if ~isfield(S.allow,'ult')   || ~isstruct(S.allow.ult),   S.allow.ult   = struct(); end

readBox = @(h) localReadAllow(h);

% -------- Current UI: casing wide allowables + pins/ret yield+ultimate --------
if isfield(S,'edAllowC_shearOut') && isgraphics(S.edAllowC_shearOut)
    v = readBox(S.edAllowC_shearOut);
    S.allow.yield.shearOut = v; S.allow.ult.shearOut = v;
end
if isfield(S,'edAllowC_netT') && isgraphics(S.edAllowC_netT)
    v = readBox(S.edAllowC_netT);
    S.allow.yield.netTension = v; S.allow.ult.netTension = v;
end
if isfield(S,'edAllowC_bear') && isgraphics(S.edAllowC_bear)
    v = readBox(S.edAllowC_bear);
    S.allow.yield.bearing = v; S.allow.ult.bearing = v;
end
if isfield(S,'edAllowC_pv') && isgraphics(S.edAllowC_pv)
    v = readBox(S.edAllowC_pv);
    S.allow.yield.pressureVessel = v; S.allow.ult.pressureVessel = v;
end

% Pins
if isfield(S,'edAllowY_ps') && isgraphics(S.edAllowY_ps)
    S.allow.yield.pinShear = readBox(S.edAllowY_ps);
end
if isfield(S,'edAllowU_ps') && isgraphics(S.edAllowU_ps)
    S.allow.ult.pinShear = readBox(S.edAllowU_ps);
end

% Ret ring
if isfield(S,'edAllowY_retSO') && isgraphics(S.edAllowY_retSO)
    S.allow.yield.ret_shearOut = readBox(S.edAllowY_retSO);
end
if isfield(S,'edAllowU_retSO') && isgraphics(S.edAllowU_retSO)
    S.allow.ult.ret_shearOut = readBox(S.edAllowU_retSO);
end
if isfield(S,'edAllowY_retNT') && isgraphics(S.edAllowY_retNT)
    S.allow.yield.ret_netTension = readBox(S.edAllowY_retNT);
end
if isfield(S,'edAllowU_retNT') && isgraphics(S.edAllowU_retNT)
    S.allow.ult.ret_netTension = readBox(S.edAllowU_retNT);
end
if isfield(S,'edAllowY_retB') && isgraphics(S.edAllowY_retB)
    S.allow.yield.ret_bearing = readBox(S.edAllowY_retB);
end
if isfield(S,'edAllowU_retB') && isgraphics(S.edAllowU_retB)
    S.allow.ult.ret_bearing = readBox(S.edAllowU_retB);
end

end

% -------- locals --------
function v = localReadAllow(h)
s = strtrim(string(get(h,'String')));
if s == ""
    v = NaN;
else
    v = str2double(s);
    if ~isfinite(v)
        v = NaN;
    end
end
end
