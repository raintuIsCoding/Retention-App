function S = readBoundsFromUI(S)

% Safely read optimization bound edit boxes if they exist
if isfield(S,'edBndTmin') && isgraphics(S.edBndTmin) && ...
   isfield(S,'edBndTmax') && isgraphics(S.edBndTmax)
    tMin = str2double(get(S.edBndTmin, 'String'));
    tMax = str2double(get(S.edBndTmax, 'String'));
    if isfinite(tMin) && isfinite(tMax) && tMin > 0 && tMax >= tMin
        S.optBounds.t = [tMin, tMax];
    end
end

if isfield(S,'edBndRowsMin') && isgraphics(S.edBndRowsMin) && ...
   isfield(S,'edBndRowsMax') && isgraphics(S.edBndRowsMax)
    rowsMin = round(str2double(get(S.edBndRowsMin, 'String')));
    rowsMax = round(str2double(get(S.edBndRowsMax, 'String')));
    if isfinite(rowsMin) && isfinite(rowsMax) && rowsMin >= 1 && rowsMax >= rowsMin
        S.optBounds.nRows = [rowsMin, rowsMax];
    end
end

if isfield(S,'edBndPPRmin') && isgraphics(S.edBndPPRmin) && ...
   isfield(S,'edBndPPRmax') && isgraphics(S.edBndPPRmax)
    pprMin = round(str2double(get(S.edBndPPRmin, 'String')));
    pprMax = round(str2double(get(S.edBndPPRmax, 'String')));
    if isfinite(pprMin) && isfinite(pprMax) && pprMin >= 1 && pprMax >= pprMin
        S.optBounds.nPinsPerRow = [pprMin, pprMax];
    end
end

if isfield(S,'edBndRowSpMin') && isgraphics(S.edBndRowSpMin) && ...
   isfield(S,'edBndRowSpMax') && isgraphics(S.edBndRowSpMax)
    rsMin = str2double(get(S.edBndRowSpMin, 'String'));
    rsMax = str2double(get(S.edBndRowSpMax, 'String'));
    if isfinite(rsMin) && isfinite(rsMax) && rsMin > 0 && rsMax >= rsMin
        S.optBounds.rowSpacing = [rsMin, rsMax];
    end
end

if isfield(S,'edBndFirstMin') && isgraphics(S.edBndFirstMin) && ...
   isfield(S,'edBndFirstMax') && isgraphics(S.edBndFirstMax)
    frMin = str2double(get(S.edBndFirstMin, 'String'));
    frMax = str2double(get(S.edBndFirstMax, 'String'));
    if isfinite(frMin) && isfinite(frMax) && frMin >= 0 && frMax >= frMin
        S.optBounds.firstRowZ = [frMin, frMax];
    end
end

if isfield(S,'edBndPinDiaMin') && isgraphics(S.edBndPinDiaMin) && ...
   isfield(S,'edBndPinDiaMax') && isgraphics(S.edBndPinDiaMax)
    pdMin = str2double(get(S.edBndPinDiaMin, 'String'));
    pdMax = str2double(get(S.edBndPinDiaMax, 'String'));
    if isfinite(pdMin) && isfinite(pdMax) && pdMin > 0 && pdMax >= pdMin
        S.optBounds.pinDia = [pdMin, pdMax];
    end
end

% Min pitch UI was removed; keep backwards compatibility:
% - If the control exists, read it
% - Otherwise leave S.minCircPitchFactor unchanged
if isfield(S,'edMinCircPitch') && isgraphics(S.edMinCircPitch)
    pitchFactor = str2double(get(S.edMinCircPitch, 'String'));
    if isfinite(pitchFactor) && pitchFactor >= 2
        S.minCircPitchFactor = pitchFactor;
    end
end

end