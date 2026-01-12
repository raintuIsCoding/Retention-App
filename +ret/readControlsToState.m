function S = readControlsToState(S)

% -------- Basic scalar inputs --------
if isfield(S,'edID') && isgraphics(S.edID)
    S.ID = str2double(get(S.edID,'String'));
end
if isfield(S,'edT') && isgraphics(S.edT)
    S.t  = str2double(get(S.edT,'String'));
end
if isfield(S,'edMEOP') && isgraphics(S.edMEOP)
    S.MEOP_psi = str2double(get(S.edMEOP,'String'));
end

% -------- FOS (new UI) or DF (old UI fallback) --------
if ~isfield(S,'FOS') || ~isstruct(S.FOS), S.FOS = struct(); end

% Defaults if missing (from your margins sheet)
if ~isfield(S.FOS,'casing_y') || ~isfinite(S.FOS.casing_y), S.FOS.casing_y = 2.00; end
if ~isfield(S.FOS,'casing_u') || ~isfinite(S.FOS.casing_u), S.FOS.casing_u = 2.00; end
if ~isfield(S.FOS,'pin_y')    || ~isfinite(S.FOS.pin_y),    S.FOS.pin_y    = 1.75; end
if ~isfield(S.FOS,'pin_u')    || ~isfinite(S.FOS.pin_u),    S.FOS.pin_u    = 2.00; end
if ~isfield(S.FOS,'ret_y')    || ~isfinite(S.FOS.ret_y),    S.FOS.ret_y    = 1.75; end
if ~isfield(S.FOS,'ret_u')    || ~isfinite(S.FOS.ret_u),    S.FOS.ret_u    = 2.00; end

% If new FOS edit boxes exist, read them
hasFOS = isfield(S,'edFOS_cy') && isgraphics(S.edFOS_cy);

if hasFOS
    S.FOS.casing_y = localReadNum(S.edFOS_cy, S.FOS.casing_y);
    S.FOS.casing_u = localReadNum(S.edFOS_cu, S.FOS.casing_u);
    S.FOS.pin_y    = localReadNum(S.edFOS_py, S.FOS.pin_y);
    S.FOS.pin_u    = localReadNum(S.edFOS_pu, S.FOS.pin_u);
    S.FOS.ret_y    = localReadNum(S.edFOS_ry, S.FOS.ret_y);
    S.FOS.ret_u    = localReadNum(S.edFOS_ru, S.FOS.ret_u);
else
    % Old UI fallback: DF fields
    if ~isfield(S,'DF_casing') || ~isfinite(S.DF_casing), S.DF_casing = 2.00; end
    if ~isfield(S,'DF_pin')    || ~isfinite(S.DF_pin),    S.DF_pin    = 1.75; end

    if isfield(S,'edDFcasing') && isgraphics(S.edDFcasing)
        S.DF_casing = localReadNum(S.edDFcasing, S.DF_casing);
    end
    if isfield(S,'edDFpin') && isgraphics(S.edDFpin)
        S.DF_pin = localReadNum(S.edDFpin, S.DF_pin);
    end

    % map into FOS so downstream code can rely on it
    S.FOS.casing_y = S.DF_casing;
    S.FOS.casing_u = S.DF_casing;
    S.FOS.pin_y    = S.DF_pin;
    S.FOS.pin_u    = S.DF_pin;
    S.FOS.ret_y    = S.DF_casing;
    S.FOS.ret_u    = S.DF_casing;
end

% -------- Rows / geometry inputs --------
if isfield(S,'edRows') && isgraphics(S.edRows)
    S.nRows = max(1, round(str2double(get(S.edRows,'String'))));
end
if isfield(S,'edPPR') && isgraphics(S.edPPR)
    S.nPinsPerRow = max(1, round(str2double(get(S.edPPR,'String'))));
end
if isfield(S,'edRowSp') && isgraphics(S.edRowSp)
    S.rowSpacing = str2double(get(S.edRowSp,'String'));
end
if isfield(S,'edFirst') && isgraphics(S.edFirst)
    S.firstRowZ = str2double(get(S.edFirst,'String'));
end

% -------- Pin diameter (dropdown preferred) --------
if isfield(S,'ddPinDia') && isgraphics(S.ddPinDia) && isfield(S,'allowedPinDias')
    idx = get(S.ddPinDia,'Value');
    idx = max(1, min(idx, numel(S.allowedPinDias)));
    S.pinDiaIdx = idx;
    S.pinDia    = S.allowedPinDias(idx);
end

% -------- View Toggle --------
if isfield(S,'btnLengthMode') && isgraphics(S.btnLengthMode)
    isFull = (get(S.btnLengthMode, 'Value') == 1);
    if isFull
        S.lengthMode = "full";
        S.mirrorPins = true;
        set(S.btnLengthMode, 'String', 'Full Casing');
    else
        S.lengthMode = "aftOnly";
        S.mirrorPins = false;
        set(S.btnLengthMode, 'String', 'Retention Detail View');
    end
end

% -------- Pattern Toggle --------
if isfield(S,'btnPattern') && isgraphics(S.btnPattern)
    if get(S.btnPattern, 'Value') == 1
        S.pinPatternMode = "alternating";
        set(S.btnPattern, 'String', "Pattern: Alternating");
    else
        S.pinPatternMode = "progressive";
        set(S.btnPattern, 'String', "Pattern: Progressive");
    end
end

% -------- Safety clamps --------
if ~isfinite(S.ID) || S.ID<=0, S.ID = 6.0; end
if ~isfinite(S.t)  || S.t<=0,  S.t = 0.225; end
if ~isfinite(S.MEOP_psi) || S.MEOP_psi<=0, S.MEOP_psi = 650; end
if ~isfinite(S.rowSpacing) || S.rowSpacing<=0, S.rowSpacing = 0.5; end
if ~isfinite(S.firstRowZ)  || S.firstRowZ<0,   S.firstRowZ  = 0.875; end
if ~isfinite(S.pinDia)     || S.pinDia<=0,     S.pinDia     = 0.25; end

% Push formatting back (only if controls exist)
if isfield(S,'edRows') && isgraphics(S.edRows), set(S.edRows,'String', num2str(S.nRows)); end
if isfield(S,'edPPR')  && isgraphics(S.edPPR),  set(S.edPPR,'String',  num2str(S.nPinsPerRow)); end
if isfield(S,'edMEOP') && isgraphics(S.edMEOP), set(S.edMEOP,'String', num2str(S.MEOP_psi)); end

% Keep allowables + bounds synced (robust versions)
S = ret.readConstraintsFromUI(S);
S = ret.readBoundsFromUI(S);

end

% -------- local --------
function v = localReadNum(h, fallback)
v = str2double(get(h,'String'));
if ~isfinite(v) || v <= 0
    v = fallback;
end
end