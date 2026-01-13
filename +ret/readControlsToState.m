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

% Defaults if missing
% - Casing uses a single FOS
% - Pins and retention ring share the same Yield/Ultimate FOS inputs
if ~isfield(S.FOS,'casing') || ~isfinite(S.FOS.casing), S.FOS.casing = 2.00; end
if ~isfield(S.FOS,'comp_y') || ~isfinite(S.FOS.comp_y), S.FOS.comp_y = 1.75; end
if ~isfield(S.FOS,'comp_u') || ~isfinite(S.FOS.comp_u), S.FOS.comp_u = 2.00; end

% If new FOS edit boxes exist, read them
hasFOS = isfield(S,'edFOS_casing') && isgraphics(S.edFOS_casing);

if hasFOS
    S.FOS.casing = localReadNum(S.edFOS_casing, S.FOS.casing);
    S.FOS.comp_y = localReadNum(S.edFOS_compY,  S.FOS.comp_y);
    S.FOS.comp_u = localReadNum(S.edFOS_compU,  S.FOS.comp_u);
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
    S.FOS.casing = S.DF_casing;
    S.FOS.comp_y = S.DF_pin;
    S.FOS.comp_u = S.DF_pin;
end

% Keep legacy DF fields in sync (used by some displays/constraints)
S.DF_casing = S.FOS.casing;
S.DF_pin    = S.FOS.comp_y;

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
        set(S.btnPattern, 'String', "Alternating");
    else
        S.pinPatternMode = "spiral";
        set(S.btnPattern, 'String', "Spiral");
    end
end

% -------- Ret ring visibility --------
if isfield(S,'btnShowRetRings') && isgraphics(S.btnShowRetRings)
    isShown = (get(S.btnShowRetRings,'Value') == 1);
    S.showRetRings = isShown;
    if isShown
        set(S.btnShowRetRings,'String','Visible');
    else
        set(S.btnShowRetRings,'String','Hidden');
    end
end

% -------- Safety clamps --------
if ~isfinite(S.ID) || S.ID<=0, S.ID = 8.0; end
if ~isfinite(S.t)  || S.t<=0,  S.t = 0.25; end
if ~isfinite(S.MEOP_psi) || S.MEOP_psi<=0, S.MEOP_psi = 850; end
if ~isfinite(S.rowSpacing) || S.rowSpacing<=0, S.rowSpacing = 0.75; end
if ~isfinite(S.firstRowZ)  || S.firstRowZ<0,   S.firstRowZ  = 1.0; end
if ~isfinite(S.pinDia)     || S.pinDia<=0,     S.pinDia     = 0.375; end

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