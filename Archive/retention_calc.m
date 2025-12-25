%% Retention-end visualizer (horizontal, pins on right, staggered rows, packing checks)
clear; clc; close all;

%% ---------- USER INPUTS ----------
ID          = 8.0;      % inner diameter of casing (in)
t           = 0.25;     % wall thickness (in)
L_casing    = 10.0;     % length of modeled region (in)

nRows       = 3;        % axial rows of pins
nPinsPerRow = 12;       % pins around circumference per row

rowSpacing  = 0.75;     % axial spacing between rows (in)
firstRowZ   = 0.75;      % axial position of first row from x=0 (in)

pinDia      = 0.375;     % pin diameter (in)
pinLen      = 2 * t;    % how far pins stick out radially (in)

% Geometric constraints
minCircPitchFactor = 2.0;  % min center spacing = factor * pinDia (circumferential)
edgeMarginEnd      = t;    % safety margin from last row to casing end (approx)

% Pattern toggle (per your definitions)
pinPatternMode = "progressive";   % "progressive" or "alternating"
% progressive  = distributed offsets between aft-row pins across all rows
% alternating  = 0, half-pitch, 0, half-pitch, ...
altStartPhase  = 0;               % 0 or 1; flips which rows get half-pitch

%% ---------- DERIVED GEOMETRY ----------
r_i = ID/2;
r_o = r_i + t;

%% ---------- CIRCUMFERENTIAL PACKING CHECK ----------
% Use an approximate radius where the pin centers live
r_center = r_o;  % consistent with placement below

circumference = 2 * pi * r_center;
minPitch      = minCircPitchFactor * pinDia;
maxPinsCirc   = floor(circumference / minPitch);

if nPinsPerRow > maxPinsCirc
    error(['Geometric packing error: nPinsPerRow = %d exceeds the maximum ' ...
           'of %d pins for an ID = %.2f in, t = %.2f in, and pinDia = %.2f in ' ...
           'with a min pitch of %.2f*pinDia.'], ...
           nPinsPerRow, maxPinsCirc, ID, t, pinDia, minCircPitchFactor);
end

%% ---------- AXIAL LENGTH CHECK ----------
lastRowX = firstRowZ + (nRows - 1) * rowSpacing;
if lastRowX + edgeMarginEnd > L_casing
    warning(['Axial pattern extends near or beyond the modeled casing length.\n' ...
             '  lastRowX = %.2f in, L_casing = %.2f in. Consider increasing L_casing\n' ...
             '  or decreasing nRows/rowSpacing.'], lastRowX, L_casing);
end

%% ---------- DRAW CASING (HOLLOW CYLINDER), AXIS ALONG X ----------
nCirc = 80;
theta = linspace(0, 2*pi, nCirc);
x     = [0, L_casing];   % axis is X

[Theta, X] = meshgrid(theta, x);

% Outer surface
Yo = r_o * cos(Theta);
Zo = r_o * sin(Theta);
Xo = X;

% Inner surface
Yi = r_i * cos(Theta);
Zi = r_i * sin(Theta);
Xi = X;

figure; hold on;
set(gcf, 'Renderer', 'opengl');

surf(Xo, Yo, Zo, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','none');
surf(Xi, Yi, Zi, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none');

% Ring wall at x = 0
Xcap0 = [zeros(1, nCirc); zeros(1, nCirc)];
Ycap0 = [r_i*cos(theta);  r_o*cos(theta)];
Zcap0 = [r_i*sin(theta);  r_o*sin(theta)];
surf(Xcap0, Ycap0, Zcap0, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');

% Ring wall at x = L_casing (pinned end)
XcapL = [L_casing*ones(1, nCirc); L_casing*ones(1, nCirc)];
YcapL = [r_i*cos(theta);          r_o*cos(theta)];
ZcapL = [r_i*sin(theta);          r_o*sin(theta)];
surf(XcapL, YcapL, ZcapL, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');

%% ---------- PIN GEOMETRY IN LOCAL COORDS ----------
nPinCirc = 30;
phi = linspace(0, 2*pi, nPinCirc);
r_pin = pinDia/2;

% Cylinder along +Y (pins pointing "right")
Yloc_side = [r_o * ones(1, nPinCirc); (r_o + pinLen)*ones(1, nPinCirc)];
Xloc_side = [r_pin*cos(phi);         r_pin*cos(phi)];
Zloc_side = [r_pin*sin(phi);         r_pin*sin(phi)];

% Tip cap at Y = r_o + pinLen
Yloc_cap = (r_o + pinLen)*ones(1, nPinCirc);
Xloc_cap = r_pin*cos(phi);
Zloc_cap = r_pin*sin(phi);

%% ---------- PLACE PINS ----------
for iRow = 1:nRows

    % Axial position of this row along X
    xRow = firstRowZ + (iRow - 1) * rowSpacing;

    % Angular pitch between adjacent pins in a row
    pitch = 2*pi / nPinsPerRow;

    % Row-to-row angular offset (per your definitions)
    switch pinPatternMode
        case "progressive"
            % Evenly distribute offsets within one pitch interval:
            % 2 rows -> 0, 0.5*pitch
            % 3 rows -> 0, (1/3)*pitch, (2/3)*pitch, etc.
            rowAngleOffset = (iRow - 1) * (pitch / nRows);

        case "alternating"
            % Alternate: 0, half-pitch, 0, half-pitch...
            rowAngleOffset = mod((iRow - 1 + altStartPhase), 2) * (pitch / 2);

        otherwise
            error('pinPatternMode must be "progressive" or "alternating".');
    end

    for j = 1:nPinsPerRow

        % Base angle for this pin + row offset
        th = 2*pi*(j-1)/nPinsPerRow + rowAngleOffset;

        % Rotate pin around X-axis (pins point along +Y in local coords)
        Y_side =  Yloc_side*cos(th) - Zloc_side*sin(th);
        Z_side =  Yloc_side*sin(th) + Zloc_side*cos(th);
        X_side =  Xloc_side + xRow;

        surf(X_side, Y_side, Z_side, ...
            'FaceColor',[0.8 0.2 0.2],'EdgeColor','none');

        % Pin tip cap
        Y_cap =  Yloc_cap*cos(th) - Zloc_cap*sin(th);
        Z_cap =  Yloc_cap*sin(th) + Zloc_cap*cos(th);
        X_cap =  Xloc_cap + xRow;

        fill3(X_cap, Y_cap, Z_cap, [0.8 0.2 0.2], 'EdgeColor', 'none');
    end
end

%% ---------- RENDER / INTERACTION SETTINGS ----------
axis equal;
axis vis3d;
axis off;
set(gca,'Visible','off');

lighting gouraud;
camlight headlight;

view(180, 0);

rotate3d on;
pan off;
zoom off;

%% ---------- OVERLAY METRICS HUD (INPUTS + OUTPUTS) ----------
fig = gcf;

% Compute outputs
totalPins = nRows * nPinsPerRow;

% ---- Inputs text (readable labels) ----
inputsText = sprintf([ ...
    'Inner Diameter: %.2f in\n' ...
    'Wall Thickness: %.3f in\n' ...
    '\n' ...
    'Axial Rows: %d\n' ...
    'Pins per Row: %d\n' ...
    'Row Spacing: %.2f in\n' ...
    'First Row Offset from End: %.2f in\n' ...
    '\n' ...
    'Pin Diameter: %.3f in\n' ...
    'Pin Engagement Length: %.3f in\n' ...
    '\n' ...
    'Pin Pattern: %s' ], ...
    ID, t, ...
    nRows, nPinsPerRow, rowSpacing, firstRowZ, ...
    pinDia, pinLen, pinPatternMode);

% ---- Outputs text (split into two buckets) ----
geomOutText = sprintf([ ...
    'Total Pins: %d' ], totalPins);

stressOutText = sprintf([ ...
    '(none yet)' ]);

% ---- Main HUD container (top-left; extended downward) ----
hud = uipanel(fig, ...
    'Units', 'normalized', ...
    'Position', [0.02 0.38 0.34 0.60], ... % extended downwards vs before
    'Title', 'Metrics', ...
    'BackgroundColor', [1 1 1], ...
    'ForegroundColor', [0 0 0], ...
    'BorderType', 'line', ...
    'HighlightColor', [0 0 0]);

% ---- Inputs panel ----
pIn = uipanel(hud, ...
    'Units', 'normalized', ...
    'Position', [0.05 0.52 0.90 0.46], ...  % upper half-ish
    'Title', 'Inputs', ...
    'ForegroundColor', [0 0 0], ...
    'BackgroundColor', [1 1 1]);

uicontrol(pIn, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.96], ...
    'String', inputsText, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1], ...
    'ForegroundColor', [0 0 0], ...
    'FontName', 'Consolas', ...
    'FontSize', 10);

% ---- Outputs container panel ----
pOut = uipanel(hud, ...
    'Units', 'normalized', ...
    'Position', [0.05 0.05 0.90 0.44], ... % lower half-ish
    'Title', 'Outputs', ...
    'ForegroundColor', [0 0 0], ...
    'BackgroundColor', [1 1 1]);

% Geometry outputs subpanel
pOutGeom = uipanel(pOut, ...
    'Units', 'normalized', ...
    'Position', [0.02 0.52 0.96 0.46], ...
    'Title', 'Geometry', ...
    'ForegroundColor', [0 0 0], ...
    'BackgroundColor', [1 1 1]);

uicontrol(pOutGeom, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.96], ...
    'String', geomOutText, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1], ...
    'ForegroundColor', [0 0 0], ...
    'FontName', 'Consolas', ...
    'FontSize', 10);

% Force / Stress outputs subpanel
pOutStress = uipanel(pOut, ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.46], ...
    'Title', 'Force / Stress', ...
    'ForegroundColor', [0 0 0], ...
    'BackgroundColor', [1 1 1]);

uicontrol(pOutStress, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.96], ...
    'String', stressOutText, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1], ...
    'ForegroundColor', [0 0 0], ...
    'FontName', 'Consolas', ...
    'FontSize', 10);

hold off;