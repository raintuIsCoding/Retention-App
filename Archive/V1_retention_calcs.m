%% Retention-end visualizer
% - Horizontal tube, pins on right
% - Spiral pin pattern: no two pins share the same angle (no pins stacked)
% - Circumferential packing check and axial length warning

clear; clc; close all;

%% ---------- USER INPUTS ----------
ID          = 8.0;      % inner diameter of casing (in)
t           = 0.30;     % wall thickness (in)
L_casing    = 10.0;     % length of modeled region (in)

nRows       = 3;        % axial rows of pins
nPinsPerRow = 12;       % pins around circumference per row

rowSpacing  = 0.75;     % axial spacing between rows (in)
firstRowZ   = 1.0;      % axial position of first row from x=0 (in)

pinDia      = 0.25;     % pin diameter (in)
pinLen      = 2 * t;    % how far pins stick out radially (in)

% Geometric constraints
minCircPitchFactor = 2.0;  % min center spacing = factor * pinDia (circumferential)
edgeMarginEnd      = t;    % safety margin from last row to casing end (approx)

%% ---------- DERIVED GEOMETRY ----------
r_i = ID/2;
r_o = r_i + t;

%% ---------- CIRCUMFERENTIAL PACKING CHECK ----------
% Approximate radius where pin centers live
r_center     = r_o;
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

%% ---------- SPIRAL ANGLE LOGIC (NO PINS ON SAME AXIS) ----------
N_totalPins = nRows * nPinsPerRow;
angle_step  = 2 * pi / N_totalPins;   % each pin gets a unique angle

%% ---------- PLACE PINS ----------
for iRow = 1:nRows
    
    % Axial position of this row along X
    xRow = firstRowZ + (iRow - 1) * rowSpacing;
    
    for j = 1:nPinsPerRow
        
        % Global pin index (0-based)
        k = (iRow - 1) * nPinsPerRow + (j - 1);
        
        % Unique angle for each pin -> spiral pattern
        th = k * angle_step;   % no two k produce the same angle in [0, 2π)
        
        % Rotate pin around X-axis
        Y_side =  Yloc_side*cos(th) - Zloc_side*sin(th);
        Z_side =  Yloc_side*sin(th) + Zloc_side*cos(th);
        X_side =  Xloc_side + xRow;
        
        surf(X_side, Y_side, Z_side, ...
            'FaceColor',[0.8 0.2 0.2],'EdgeColor','none');
        
        % Tip cap
        Y_cap =  Yloc_cap*cos(th) - Zloc_cap*sin(th);
        Z_cap =  Yloc_cap*sin(th) + Zloc_cap*cos(th);
        X_cap =  Xloc_cap + xRow;
        
        fill3(X_cap, Y_cap, Z_cap, [0.8 0.2 0.2], 'EdgeColor', 'none');
    end
end

%% ---------- RENDER SETTINGS ----------
axis equal;
axis off;
set(gca,'Visible','off');
lighting gouraud;
camlight headlight;

% Your preferred view
view(180, 0);

hold off;
