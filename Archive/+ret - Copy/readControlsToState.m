function S = readControlsToState(S)

S.ID         = str2double(get(S.edID,'String'));
S.t          = str2double(get(S.edT,'String'));
S.L_casing   = str2double(get(S.edL,'String'));
S.MEOP_psi   = str2double(get(S.edMEOP,'String'));
S.DF         = str2double(get(S.edDF,'String'));
S.nRows      = max(1, round(str2double(get(S.edRows,'String'))));
S.nPinsPerRow= max(1, round(str2double(get(S.edPPR,'String'))));
S.rowSpacing = str2double(get(S.edRowSp,'String'));
S.firstRowZ  = str2double(get(S.edFirst,'String'));
S.pinDia     = str2double(get(S.edPinDia,'String'));
S.pinLen     = str2double(get(S.edPinLen,'String'));

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
if ~isfinite(S.L_casing) || S.L_casing<=0, S.L_casing = 10.0; end
if ~isfinite(S.MEOP_psi) || S.MEOP_psi<=0, S.MEOP_psi = 900; end
if ~isfinite(S.DF)       || S.DF<=0,       S.DF = 1.5; end
if ~isfinite(S.rowSpacing) || S.rowSpacing<=0, S.rowSpacing = 0.75; end
if ~isfinite(S.firstRowZ)  || S.firstRowZ<0,   S.firstRowZ  = 0.75; end
if ~isfinite(S.pinDia)     || S.pinDia<=0,     S.pinDia     = 0.375; end
if ~isfinite(S.pinLen)     || S.pinLen<=0,     S.pinLen     = 2*S.t; end

% Push formatting back
set(S.edRows,'String', num2str(S.nRows));
set(S.edPPR,'String',  num2str(S.nPinsPerRow));
set(S.edMEOP,'String', num2str(S.MEOP_psi));
set(S.edDF,'String',   num2str(S.DF));

% Also keep targets/bounds synced
S = ret.readConstraintsFromUI(S);
S = ret.readBoundsFromUI(S);

end
