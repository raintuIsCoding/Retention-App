function S = readControlsToState(S)

S.ID         = str2double(get(S.edID,'String'));
S.t          = str2double(get(S.edT,'String'));
S.MEOP_psi   = str2double(get(S.edMEOP,'String'));
S.DF_casing  = str2double(get(S.edDFcasing,'String'));
S.DF_pin     = str2double(get(S.edDFpin,'String'));
S.nRows      = max(1, round(str2double(get(S.edRows,'String'))));
S.nPinsPerRow= max(1, round(str2double(get(S.edPPR,'String'))));
S.rowSpacing = str2double(get(S.edRowSp,'String'));
S.firstRowZ  = str2double(get(S.edFirst,'String'));

% --- Pin diameter (dropdown/popupmenu preferred) ---
if isfield(S,'ddPinDia') && isgraphics(S.ddPinDia)
    idx = get(S.ddPinDia,'Value');
    idx = max(1, min(idx, numel(S.allowedPinDias)));
    S.pinDiaIdx = idx;
    S.pinDia = S.allowedPinDias(idx);

elseif isfield(S,'edPinDia') && isgraphics(S.edPinDia)
    % legacy editbox version
    S.pinDia = str2double(get(S.edPinDia,'String'));
    [~, idx] = min(abs(S.allowedPinDias - S.pinDia));
    S.pinDiaIdx = idx;

else
    % fallback
    if ~isfield(S,'pinDiaIdx') || isempty(S.pinDiaIdx)
        [~, idx] = min(abs(S.allowedPinDias - S.pinDia));
        S.pinDiaIdx = idx;
    end
end

% Modeled length toggle
isFull = (get(S.btnLengthMode, 'Value') == 1);
if isFull
    S.lengthMode = "full";
    S.L_casing   = 24.0;
    S.mirrorPins = true;
    set(S.btnLengthMode, 'String', "Length: 96 in (both ends)");
else
    S.lengthMode = "aftOnly";
    S.L_casing   = 10.0;
    S.mirrorPins = false;
    set(S.btnLengthMode, 'String', "Length: 10 in (aft only)");
end

% Derived: pin engagement length is wall thickness + retention ring thickness
if ~isfield(S,'retRingThk') || ~isfinite(S.retRingThk) || S.retRingThk<=0
    S.retRingThk = 0.25;
end
S.pinLen = S.t + S.retRingThk;

% Pattern toggle
if get(S.btnPattern, 'Value') == 1
    S.pinPatternMode = "alternating";
    set(S.btnPattern, 'String', "Pattern: Alternating");
else
    S.pinPatternMode = "progressive";
    set(S.btnPattern, 'String', "Pattern: Progressive");
end

% Derived
S.edgeMarginEnd = S.t;

% Clamp bad values
if ~isfinite(S.ID) || S.ID<=0, S.ID = 8.0; end
if ~isfinite(S.t)  || S.t<=0,  S.t = 0.25; end
if ~isfinite(S.MEOP_psi) || S.MEOP_psi<=0, S.MEOP_psi = 900; end
if ~isfinite(S.DF_casing) || S.DF_casing<=0, S.DF_casing = 1.5; end
if ~isfinite(S.DF_pin)    || S.DF_pin<=0,    S.DF_pin    = 2.0; end
if ~isfinite(S.rowSpacing) || S.rowSpacing<=0, S.rowSpacing = 0.75; end
if ~isfinite(S.firstRowZ)  || S.firstRowZ<0,   S.firstRowZ  = 0.75; end
if ~isfinite(S.pinDia)     || S.pinDia<=0,     S.pinDia     = 0.375; end

% Push formatting back
set(S.edRows,'String', num2str(S.nRows));
set(S.edPPR,'String',  num2str(S.nPinsPerRow));
set(S.edMEOP,'String', num2str(S.MEOP_psi));
set(S.edDFcasing,'String', num2str(S.DF_casing));
set(S.edDFpin,'String',    num2str(S.DF_pin));

% Also keep targets/bounds synced
S = ret.readConstraintsFromUI(S);
S = ret.readBoundsFromUI(S);

end