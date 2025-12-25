function S = readBoundsFromUI(S)

tMin = str2double(get(S.edBndTmin, 'String'));
tMax = str2double(get(S.edBndTmax, 'String'));
if isfinite(tMin) && isfinite(tMax) && tMin > 0 && tMax >= tMin
    S.optBounds.t = [tMin, tMax];
end

rowsMin = round(str2double(get(S.edBndRowsMin, 'String')));
rowsMax = round(str2double(get(S.edBndRowsMax, 'String')));
if isfinite(rowsMin) && isfinite(rowsMax) && rowsMin >= 1 && rowsMax >= rowsMin
    S.optBounds.nRows = [rowsMin, rowsMax];
end

pprMin = round(str2double(get(S.edBndPPRmin, 'String')));
pprMax = round(str2double(get(S.edBndPPRmax, 'String')));
if isfinite(pprMin) && isfinite(pprMax) && pprMin >= 3 && pprMax >= pprMin
    S.optBounds.nPinsPerRow = [pprMin, pprMax];
end

rsMin = str2double(get(S.edBndRowSpMin, 'String'));
rsMax = str2double(get(S.edBndRowSpMax, 'String'));
if isfinite(rsMin) && isfinite(rsMax) && rsMin > 0 && rsMax >= rsMin
    S.optBounds.rowSpacing = [rsMin, rsMax];
end

frMin = str2double(get(S.edBndFirstMin, 'String'));
frMax = str2double(get(S.edBndFirstMax, 'String'));
if isfinite(frMin) && isfinite(frMax) && frMin > 0 && frMax >= frMin
    S.optBounds.firstRowZ = [frMin, frMax];
end

pdMin = str2double(get(S.edBndPinDiaMin, 'String'));
pdMax = str2double(get(S.edBndPinDiaMax, 'String'));
if isfinite(pdMin) && isfinite(pdMax) && pdMin > 0 && pdMax >= pdMin
    S.optBounds.pinDia = [pdMin, pdMax];
end

pitchFactor = str2double(get(S.edMinCircPitch, 'String'));
if isfinite(pitchFactor) && pitchFactor >= 2
    S.minCircPitchFactor = pitchFactor;
end

end
