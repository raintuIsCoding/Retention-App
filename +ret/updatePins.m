function S = updatePins(S, r_o)

% How many pin sets are we drawing? (aft only vs both ends)
nEnds = 1 + double(isfield(S,'mirrorPins') && S.mirrorPins);

nUsed = S.nRows * S.nPinsPerRow * nEnds;
if nUsed > S.maxPins
    warning('Pin pool too small: need %d pins but maxPins=%d. Increase S.maxPins.', nUsed, S.maxPins);
    nUsed = S.maxPins;
end

phi   = linspace(0, 2*pi, S.nPinCirc);
r_pin = S.pinDia/2;

Yloc_side = [r_o*ones(1,S.nPinCirc); (r_o+S.pinLen)*ones(1,S.nPinCirc)];
Xloc_side = [r_pin*cos(phi);         r_pin*cos(phi)];
Zloc_side = [r_pin*sin(phi);         r_pin*sin(phi)];

Yloc_cap = (r_o + S.pinLen)*ones(1,S.nPinCirc);
Xloc_cap = r_pin*cos(phi);
Zloc_cap = r_pin*sin(phi);

pitch = 2*pi / S.nPinsPerRow;
idx   = 0;

for iRow = 1:S.nRows
    xRow = S.firstRowZ + (iRow - 1)*S.rowSpacing;

    switch S.pinPatternMode
        case "progressive"
            rowAngleOffset = (iRow - 1) * (pitch / S.nRows);
        case "alternating"
            rowAngleOffset = mod((iRow - 1 + S.altStartPhase), 2) * (pitch/2);
        otherwise
            rowAngleOffset = 0;
    end

    % Axial positions for this row:
    % - always include aft end (xRow)
    % - optionally include forward end (L - xRow)
    xPositions = xRow;
    if isfield(S,'mirrorPins') && S.mirrorPins
        xPositions = [xPositions, (S.L_casing - xRow)];
    end

    for e = 1:numel(xPositions)
        xHere = xPositions(e);

        for j = 1:S.nPinsPerRow
            idx = idx + 1;
            if idx > nUsed
                break;
            end

            th = 2*pi*(j-1)/S.nPinsPerRow + rowAngleOffset;

            % Side surface
            Y_side =  Yloc_side*cos(th) - Zloc_side*sin(th);
            Z_side =  Yloc_side*sin(th) + Zloc_side*cos(th);
            X_side =  Xloc_side + xHere;
            set(S.pinSide(idx), 'XData', X_side, 'YData', Y_side, 'ZData', Z_side, 'Visible','on');

            % Cap
            Y_cap =  Yloc_cap*cos(th) - Zloc_cap*sin(th);
            Z_cap =  Yloc_cap*sin(th) + Zloc_cap*cos(th);
            X_cap =  Xloc_cap + xHere;
            set(S.pinCap(idx), 'XData', X_cap, 'YData', Y_cap, 'ZData', Z_cap, 'Visible','on');
        end

        if idx >= nUsed
            break;
        end
    end

    if idx >= nUsed
        break;
    end
end

% Hide unused pooled graphics objects
for k = (nUsed+1):S.maxPins
    set(S.pinSide(k), 'Visible','off');
    set(S.pinCap(k),  'Visible','off');
end

end