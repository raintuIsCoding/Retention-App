function S = readControlsToState(S)

S.ID         = str2double(get(S.edID,'String'));
S.t          = str2double(get(S.edT,'String'));
S.MEOP_psi   = str2double(get(S.edMEOP,'String'));
S.DF         = str2double(get(S.edDF,'String'));
S.nRows      = max(1, round(str2double(get(S.edRows,'String'))));
S.nPinsPerRow= max(1, round(str2double(get(S.edPPR,'String'))));
S.rowSpacing = str2double(get(S.edRowSp,'String'));
S.firstRowZ  = str2double(get(S.edFirst,'String'));
S.pinDia     = str2double(get(S.edPinDia,'String'));

% Modeled length toggle
isFull = (get(S.btnLengthMode, 'Value') == 1);
if isFull
    S.lengthMode = "full";
    S.L_casing   = 96.0;
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
if ~isfinite(S.DF)       || S.DF<=0,       S.DF = 1.5; end
if ~isfinite(S.rowSpacing) || S.rowSpacing<=0, S.rowSpacing = 0.75; end
if ~isfinite(S.firstRowZ)  || S.firstRowZ<0,   S.firstRowZ  = 0.75; end
if ~isfinite(S.pinDia)     || S.pinDia<=0,     S.pinDia     = 0.375; end

% Push formatting back
set(S.edRows,'String', num2str(S.nRows));
set(S.edPPR,'String',  num2str(S.nPinsPerRow));
set(S.edMEOP,'String', num2str(S.MEOP_psi));
set(S.edDF,'String',   num2str(S.DF));

% Also keep targets/bounds synced
S = ret.readConstraintsFromUI(S);
S = ret.readBoundsFromUI(S);

end