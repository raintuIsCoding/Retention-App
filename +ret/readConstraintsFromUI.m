function S = readConstraintsFromUI(S)

% Ensure structs exist
if ~isfield(S,'allow') || ~isstruct(S.allow), S.allow = struct(); end
if ~isfield(S.allow,'yield') || ~isstruct(S.allow.yield), S.allow.yield = struct(); end
if ~isfield(S.allow,'ult')   || ~isstruct(S.allow.ult),   S.allow.ult   = struct(); end
if ~isfield(S,'targets') || ~isstruct(S.targets), S.targets = struct(); end

% -------- New UI: Allowables (Yield + Ultimate) --------
hasAllowUI = isfield(S,'edAllowY_shearOut') && isgraphics(S.edAllowY_shearOut);

if hasAllowUI
    % helper
    readBox = @(h) localReadAllow(h);

    % casing allowables
    if isfield(S,'edAllowY_shearOut'),       S.allow.yield.shearOut       = readBox(S.edAllowY_shearOut); end
    if isfield(S,'edAllowU_shearOut'),       S.allow.ult.shearOut         = readBox(S.edAllowU_shearOut); end
    if isfield(S,'edAllowY_netT'),           S.allow.yield.netTension     = readBox(S.edAllowY_netT);     end
    if isfield(S,'edAllowU_netT'),           S.allow.ult.netTension       = readBox(S.edAllowU_netT);     end
    if isfield(S,'edAllowY_bear'),           S.allow.yield.bearing        = readBox(S.edAllowY_bear);     end
    if isfield(S,'edAllowU_bear'),           S.allow.ult.bearing          = readBox(S.edAllowU_bear);     end
    if isfield(S,'edAllowY_pv'),             S.allow.yield.pressureVessel = readBox(S.edAllowY_pv);       end
    if isfield(S,'edAllowU_pv'),             S.allow.ult.pressureVessel   = readBox(S.edAllowU_pv);       end

    % pins allowables
    if isfield(S,'edAllowY_ps'),             S.allow.yield.pinShear       = readBox(S.edAllowY_ps);       end
    if isfield(S,'edAllowU_ps'),             S.allow.ult.pinShear         = readBox(S.edAllowU_ps);       end

    % ret ring allowables (may exist)
    if isfield(S,'edAllowY_retSO'),          S.allow.yield.ret_shearOut   = readBox(S.edAllowY_retSO);    end
    if isfield(S,'edAllowU_retSO'),          S.allow.ult.ret_shearOut     = readBox(S.edAllowU_retSO);    end
    if isfield(S,'edAllowY_retNT'),          S.allow.yield.ret_netTension = readBox(S.edAllowY_retNT);    end
    if isfield(S,'edAllowU_retNT'),          S.allow.ult.ret_netTension   = readBox(S.edAllowU_retNT);    end
    if isfield(S,'edAllowY_retB'),           S.allow.yield.ret_bearing    = readBox(S.edAllowY_retB);     end
    if isfield(S,'edAllowU_retB'),           S.allow.ult.ret_bearing      = readBox(S.edAllowU_retB);     end

    % Keep optimization pipeline alive:
    % Map ultimate allowables -> legacy target maxima
    S.targets.shearOut_max        = localFiniteOr(S.allow.ult.shearOut,       S.targets, 'shearOut_max',  Inf);
    S.targets.netTension_max      = localFiniteOr(S.allow.ult.netTension,     S.targets, 'netTension_max',Inf);
    S.targets.pinShear_max        = localFiniteOr(S.allow.ult.pinShear,       S.targets, 'pinShear_max',  Inf);
    S.targets.bearing_max         = localFiniteOr(S.allow.ult.bearing,        S.targets, 'bearing_max',   Inf);
    S.targets.pressure_vessel_max = localFiniteOr(S.allow.ult.pressureVessel, S.targets, 'pressure_vessel_max', Inf);

    return;
end

% -------- Old UI fallback: Target Constraints --------
hasOldTargets = isfield(S,'edTgtShearOut') && isgraphics(S.edTgtShearOut);
if hasOldTargets
    S.targets.shearOut_max        = str2double(get(S.edTgtShearOut, 'String'));
    S.targets.netTension_max      = str2double(get(S.edTgtNetTens,  'String'));
    S.targets.pinShear_max        = str2double(get(S.edTgtPinShear, 'String'));
    S.targets.bearing_max         = str2double(get(S.edTgtBearing,  'String'));
    S.targets.pressure_vessel_max = str2double(get(S.edTgtPressureVessel, 'String'));
end

end

% -------- locals --------
function v = localReadAllow(h)
s = strtrim(string(get(h,'String')));
if s == ""
    v = NaN;
else
    v = str2double(s);
    if ~isfinite(v), v = NaN; end
end
end

function out = localFiniteOr(val, Sstruct, field, fallback)
if isfinite(val)
    out = val;
elseif isfield(Sstruct,field) && isfinite(Sstruct.(field))
    out = Sstruct.(field);
else
    out = fallback;
end
end