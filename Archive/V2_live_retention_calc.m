%% Retention-end visualizer with LIVE UI + live geometry updates (FULL FIXED)
% Notes:
% - MEOP + Design Factor are added as live inputs
% - Force/Stress outputs now show Design Pressure, Axial End Force, Avg Force per Pin
% - Pin pattern toggle button included (progressive vs alternating)
% - Geometry and pins update by updating existing graphics objects (fast)

clear; clc; close all;

%% ---------- DEFAULT INPUTS ----------
S.ID          = 8.0;      % inner diameter of casing (in)
S.t           = 0.25;     % wall thickness (in)
S.L_casing    = 10.0;     % length of modeled region (in)

S.nRows       = 3;        % axial rows of pins
S.nPinsPerRow = 12;       % pins around circumference per row

S.rowSpacing  = 0.5;     % axial spacing between rows (in)
S.firstRowZ   = 0.75;     % axial position of first row from x=0 (in)

S.pinDia      = 0.375;    % pin diameter (in)
S.pinLen      = 2 * S.t;  % how far pins stick out radially (in)

% Loads inputs
S.MEOP_psi   = 900;   % psi
S.DF         = 1.5;   % design factor (MEOP multiplier)

% Geometric constraints
S.minCircPitchFactor = 2.0;  % min center spacing = factor * pinDia (circumferential)
S.edgeMarginEnd      = S.t;  % safety margin from last row to casing end (approx)

% Pattern toggle (per your definitions)
S.pinPatternMode = "progressive";   % "progressive" or "alternating"
S.altStartPhase  = 0;               % 0 or 1

%% ---------- RENDER CONSTANTS ----------
S.nCirc    = 150;  % casing mesh resolution
S.nPinCirc = 30;   % pin mesh resolution
S.maxPins  = 400;  % pool size (increase if you expect more than this)

%% ---------- FIGURE / AXES ----------
S.fig = figure('Renderer','opengl');
S.ax = axes( ...
    'Parent', S.fig, ...
    'Units', 'normalized', ...
    'Position', [0.50 0.10 0.45 0.80]);  % shifted right, ~2/3 scale

hold(S.ax, 'on');

axis(S.ax, 'equal');
axis(S.ax, 'vis3d');
axis(S.ax, 'off');
set(S.ax, 'Visible','off');

lighting(S.ax, 'gouraud');
camlight(S.ax, 'headlight');
view(S.ax, 180, 0);

rotate3d(S.fig, 'on');
pan(S.fig, 'off');
zoom(S.fig, 'off');

%% ---------- DRAW CASING (create handles once) ----------
theta = linspace(0, 2*pi, S.nCirc);
x2    = [0, 1]; % placeholder; will be set in updateCasingSurfaces()
[Theta, X] = meshgrid(theta, x2);

S.hCasingOuter = surf(S.ax, X, X, X, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','none');
S.hCasingInner = surf(S.ax, X, X, X, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none');
S.hCap0        = surf(S.ax, X, X, X, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none');
S.hCapL        = surf(S.ax, X, X, X, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none');

%% ---------- CREATE PIN GRAPHICS POOL (handles once) ----------
% Each pin has: a side surf + a cap patch
S.pinSide = gobjects(S.maxPins,1);
S.pinCap  = gobjects(S.maxPins,1);

for k = 1:S.maxPins
    S.pinSide(k) = surf(S.ax, nan(2,S.nPinCirc), nan(2,S.nPinCirc), nan(2,S.nPinCirc), ...
        'FaceColor',[0.8 0.2 0.2], 'EdgeColor','none', 'Visible','off');
    S.pinCap(k) = patch(S.ax, nan(1,S.nPinCirc), nan(1,S.nPinCirc), nan(1,S.nPinCirc), ...
        [0.8 0.2 0.2], 'EdgeColor','none', 'Visible','off');
end

%% ---------- HUD + LIVE INPUT CONTROLS ----------
S = buildHUDandControls(S);

%% ---------- INITIAL UPDATE ----------
S = updateAll(S);

% Store state
guidata(S.fig, S);

%% ---------- CALLBACK ----------
function onAnyInputChanged(src, ~)
    S = guidata(src);
    S = readControlsToState(S);
    S = updateAll(S);
    guidata(S.fig, S);
end

%% ---------- BUILD HUD / CONTROLS ----------
function S = buildHUDandControls(S)
    fig = S.fig;

    % Main container panel (top-left; extended downward, wider)
    S.hud = uipanel(fig, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.05 0.42 0.90], ...
        'Title', 'Metrics', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'BorderType', 'line', ...
        'HighlightColor', [0 0 0]);

    % Inputs panel (interactive)
    S.pIn = uipanel(S.hud, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.52 0.90 0.46], ...
        'Title', 'Inputs', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);

    % Outputs container
    S.pOut = uipanel(S.hud, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.05 0.90 0.44], ...
        'Title', 'Outputs', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);

    % Geometry outputs subpanel
    S.pOutGeom = uipanel(S.pOut, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.02 0.47 0.96], ...
        'Title', 'Geometry', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);

    S.txtGeom = uicontrol(S.pOutGeom, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.02 0.96 0.96], ...
        'String', '', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'FontName', 'Consolas', ...
        'FontSize', 10);

    % Force / Stress outputs subpanel
    S.pOutStress = uipanel(S.pOut, ...
        'Units', 'normalized', ...
        'Position', [0.51 0.02 0.47 0.96], ...
        'Title', 'Force / Stress', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);

    S.txtStress = uicontrol(S.pOutStress, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.02 0.96 0.96], ...
        'String', '', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'FontName', 'Consolas', ...
        'FontSize', 10);

    % ---- Controls layout inside Inputs panel ----
    % Normalized positions
    yTop = 0.90; dy = 0.075;
    xL  = 0.04; wL = 0.56;
    xC  = 0.62; wC = 0.34;
    h   = 0.07;

    % Helper to create label
    mkLabel = @(txt, yy) uicontrol(S.pIn, 'Style','text', 'Units','normalized', ...
        'Position',[xL yy wL h], 'String',txt, 'HorizontalAlignment','left', ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], 'FontName','Consolas', 'FontSize',9);

    % Helper to create edit
    mkEdit = @(val, yy, tag) uicontrol(S.pIn, 'Style','edit', 'Units','normalized', ...
        'Position',[xC yy wC h], 'String',num2str(val), 'Tag',tag, ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], ...
        'Callback',@onAnyInputChanged);

    % Build controls (added MEOP + DF, shifted indices)
    mkLabel('Inner Diameter (in)', yTop);                  S.edID      = mkEdit(S.ID,         yTop,        'ID');
    mkLabel('Wall Thickness (in)', yTop-dy);               S.edT       = mkEdit(S.t,          yTop-dy,     't');
    mkLabel('Modeled Length (in)', yTop-2*dy);             S.edL       = mkEdit(S.L_casing,   yTop-2*dy,   'L');

    mkLabel('MEOP (psi)', yTop-3*dy);                      S.edMEOP    = mkEdit(S.MEOP_psi,   yTop-3*dy,   'MEOP_psi');
    mkLabel('Design Factor (xMEOP)', yTop-4*dy);           S.edDF      = mkEdit(S.DF,         yTop-4*dy,   'DF');

    mkLabel('Axial Rows', yTop-5*dy);                      S.edRows    = mkEdit(S.nRows,      yTop-5*dy,   'nRows');
    mkLabel('Pins per Row', yTop-6*dy);                    S.edPPR     = mkEdit(S.nPinsPerRow,yTop-6*dy,   'nPinsPerRow');

    mkLabel('Row Spacing (in)', yTop-7*dy);                S.edRowSp   = mkEdit(S.rowSpacing, yTop-7*dy,   'rowSpacing');
    mkLabel('First Row Offset from End (in)', yTop-8*dy);  S.edFirst   = mkEdit(S.firstRowZ,  yTop-8*dy,   'firstRowZ');

    mkLabel('Pin Diameter (in)', yTop-9*dy);               S.edPinDia  = mkEdit(S.pinDia,     yTop-9*dy,   'pinDia');
    mkLabel('Pin Engagement Length (in)', yTop-10*dy);     S.edPinLen  = mkEdit(S.pinLen,     yTop-10*dy,  'pinLen');

    % Pattern toggle button
    mkLabel('Pin Pattern', yTop-11*dy);

    isAlt = (S.pinPatternMode == "alternating");
    btnStr = "Pattern: Progressive";
    if isAlt, btnStr = "Pattern: Alternating"; end

    S.btnPattern = uicontrol(S.pIn, 'Style','togglebutton', 'Units','normalized', ...
        'Position',[xC (yTop-11*dy) wC h], ...
        'String', btnStr, ...
        'Value', isAlt, ...
        'BackgroundColor',[0.95 0.95 0.95], ...
        'ForegroundColor',[0 0 0], ...
        'FontName','Consolas', ...
        'FontSize', 9, ...
        'Callback',@onAnyInputChanged);
end

%% ---------- READ CONTROLS -> STATE ----------
function S = readControlsToState(S)
    % Numeric reads
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

    % Pattern (toggle button)
    if get(S.btnPattern, 'Value') == 1
        S.pinPatternMode = "alternating";
    else
        S.pinPatternMode = "progressive";
    end

    if S.pinPatternMode == "alternating"
        set(S.btnPattern, 'String', "Pattern: Alternating");
    else
        set(S.btnPattern, 'String', "Pattern: Progressive");
    end

    % Derived constraints
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

    % Push formatting back to UI
    set(S.edRows,'String', num2str(S.nRows));
    set(S.edPPR,'String',  num2str(S.nPinsPerRow));
    set(S.edMEOP,'String', num2str(S.MEOP_psi));
    set(S.edDF,'String',   num2str(S.DF));
end

%% ---------- UPDATE EVERYTHING ----------
function S = updateAll(S)
    % Derived geometry
    r_i = S.ID/2;
    r_o = r_i + S.t;

    % Packing check
    r_center = r_o;
    circumference = 2*pi*r_center;
    minPitch = S.minCircPitchFactor * S.pinDia;
    maxPinsCirc = floor(circumference / minPitch);
    packingOK = (S.nPinsPerRow <= maxPinsCirc);

    % Axial check
    lastRowX = S.firstRowZ + (S.nRows - 1)*S.rowSpacing;
    axialOK = (lastRowX + S.edgeMarginEnd <= S.L_casing);

    % Update casing + pins
    S = updateCasingSurfaces(S, r_i, r_o);
    S = updatePins(S, r_o);

    % Density
    retention_ring_thickness = 0.25; % Set by mandrel thickness
    Density_retention_ring = 0.0975; % lb/in^3 for 6061 aluminum

    % Geometry outputs
    totalPins = S.nRows * S.nPinsPerRow;
    Volume_pins = S.nPinsPerRow*S.nRows*pi*retention_ring_thickness*(S.pinDia/2)^2; % to subtract from retention ring
    Volume_retention_ring = (pi/4)*((S.ID^2) - (S.ID - 2*retention_ring_thickness)^2)*(S.edgeMarginEnd*2 + S.rowSpacing*(S.nRows-1)) - (Volume_pins);
    Mass_retention_ring = Volume_retention_ring*Density_retention_ring*2; % lbs, x2 for both ends

    geomOutText = sprintf([ ...
        'Total Pins: %d\n' ...
        'Last Row X: %.2f in\n' ...
        'Max Pins/Row (min pitch): %d\n' ...
        'Packing OK: %s\n' ...
        'Axial OK: %s\n' ...
        '\n' ...
        'Retention Ring Mass: %.2f lb\n' ], ...
        totalPins, lastRowX, maxPinsCirc, tfStr(packingOK), tfStr(axialOK), Mass_retention_ring);

    set(S.txtGeom, 'String', geomOutText);

    %%%%%%%%%%%%%%%% -------- AREA -------- %%%%%%%%%%%%%%%%
    n = 1:S.nRows;
    A_bore   = pi * (S.ID/2)^2;       % in^2
    A_ShearOut = sum((S.rowSpacing*(n-1) + S.firstRowZ) * S.t * 2 * S.nPinsPerRow); % Sum all shearout areas
    A_tension = pi*0.25*((S.ID + 2*S.t)^2 - S.ID^2) - S.pinDia*S.nPinsPerRow*S.nRows*S.t;

    %%%%%%%%%%%%%%%% ----- Load calcs ----- %%%%%%%%%%%%%%%%
    p_design = S.MEOP_psi * S.DF;     % psi
    F_axial  = p_design * A_bore;     % lbf
    F_perPin = F_axial / totalPins;   % lbf average
    Stress_ShearOut = F_axial/A_ShearOut/1000; % KSI - Notes: Axial stress / Net shear out area (2 * distance to edge * casing thickness). May break down at larger pin diameters
    Net_Tension = (F_axial/A_tension)/1000;    % KSI - Assume total number of pins contribute to net tension since close relative proximity requires a conservative approximation. However, if the pin rows are spaced further apart such that the pins in seperate rows are further apart than they are from pins in the same row, then this conservative assumption is redundant. This case is not implemented but could be in the future.
    Pin_Shear = (F_perPin/(pi*0.25*(S.pinDia)^2))/1000; % KSI - Single plane shear on pin
    Bearing = (F_perPin/(S.pinDia*S.t))/1000; % KSI - Bearing per pin
    Hoop = (p_design*(((S.ID+2*S.t)^2 + (S.ID)^2)/((S.ID+2*S.t)^2 - (S.ID)^2)))/1000; % KSI - Shrimply get shmogged by THICC wall assumption
    Axial = Hoop/2; % KSI - get conserved nerd

    %%%%%%%%%%%%%%%% ----- FOS calcs ----- %%%%%%%%%%%%%%%%
    Pin_Shear_Strength = 43.5; % KSI
    Pin_Shear_Strength_Ultimate = 75; % KSI
    Pin_Shear_FOS = (Pin_Shear_Strength/Pin_Shear) - 2; % Needs to be positive
    

    stressOutText = sprintf([ ...
        'Design Pressure: %.1f psi\n' ...
        'Axial End Force: %.0f lbf\n' ...
        'Avg Force per Pin: %.1f lbf\n' ...
        'Shear Out Stress: %.4f KSI\n' ...
        'Net Tension: %.4f KSI\n' ...
        'Pin Shear: %.4f KSI\n' ...
        'Bearing: %.4f KSI\n' ...
        'Hoop: %.4f KSI\n' ...
        'Axial: %.4f KSI\n' ...
        '\n' ...
        'Pin Shear FOS: %.4f \n' ], ...
        p_design, F_axial, F_perPin, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);

    set(S.txtStress, 'String', stressOutText);

    drawnow limitrate;
end

function s = tfStr(tf)
    if tf, s = 'YES'; else, s = 'NO'; end
end

%% ---------- UPDATE CASING SURF HANDLES ----------
function S = updateCasingSurfaces(S, r_i, r_o)
    theta = linspace(0, 2*pi, S.nCirc);
    x     = [0, S.L_casing];
    [Theta, X] = meshgrid(theta, x);

    Yo = r_o*cos(Theta);
    Zo = r_o*sin(Theta);

    Yi = r_i*cos(Theta);
    Zi = r_i*sin(Theta);

    set(S.hCasingOuter, 'XData', X, 'YData', Yo, 'ZData', Zo);
    set(S.hCasingInner, 'XData', X, 'YData', Yi, 'ZData', Zi);

    % Cap at x=0
    Xcap0 = [zeros(1,S.nCirc); zeros(1,S.nCirc)];
    Ycap0 = [r_i*cos(theta);   r_o*cos(theta)];
    Zcap0 = [r_i*sin(theta);   r_o*sin(theta)];
    set(S.hCap0, 'XData', Xcap0, 'YData', Ycap0, 'ZData', Zcap0);

    % Cap at x=L
    XcapL = [S.L_casing*ones(1,S.nCirc); S.L_casing*ones(1,S.nCirc)];
    YcapL = [r_i*cos(theta);            r_o*cos(theta)];
    ZcapL = [r_i*sin(theta);            r_o*sin(theta)];
    set(S.hCapL, 'XData', XcapL, 'YData', YcapL, 'ZData', ZcapL);
end

%% ---------- UPDATE PIN SURF/PATCH HANDLES ----------
function S = updatePins(S, r_o)
    nUsed = S.nRows * S.nPinsPerRow;

    if nUsed > S.maxPins
        warning('Pin pool too small: need %d pins but maxPins=%d. Increase S.maxPins.', nUsed, S.maxPins);
        nUsed = S.maxPins;
    end

    % Local pin mesh (updated live for pinDia/pinLen)
    phi = linspace(0, 2*pi, S.nPinCirc);
    r_pin = S.pinDia/2;

    % Cylinder along +Y
    Yloc_side = [r_o*ones(1,S.nPinCirc); (r_o+S.pinLen)*ones(1,S.nPinCirc)];
    Xloc_side = [r_pin*cos(phi);         r_pin*cos(phi)];
    Zloc_side = [r_pin*sin(phi);         r_pin*sin(phi)];

    % Tip cap
    Yloc_cap = (r_o + S.pinLen)*ones(1,S.nPinCirc);
    Xloc_cap = r_pin*cos(phi);
    Zloc_cap = r_pin*sin(phi);

    pitch = 2*pi / S.nPinsPerRow;

    idx = 0;
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

        for j = 1:S.nPinsPerRow
            idx = idx + 1;
            if idx > nUsed, break; end

            th = 2*pi*(j-1)/S.nPinsPerRow + rowAngleOffset;

            % rotate around X-axis
            Y_side =  Yloc_side*cos(th) - Zloc_side*sin(th);
            Z_side =  Yloc_side*sin(th) + Zloc_side*cos(th);
            X_side =  Xloc_side + xRow;

            set(S.pinSide(idx), 'XData', X_side, 'YData', Y_side, 'ZData', Z_side, 'Visible','on');

            Y_cap =  Yloc_cap*cos(th) - Zloc_cap*sin(th);
            Z_cap =  Yloc_cap*sin(th) + Zloc_cap*cos(th);
            X_cap =  Xloc_cap + xRow;

            set(S.pinCap(idx), 'XData', X_cap, 'YData', Y_cap, 'ZData', Z_cap, 'Visible','on');
        end
    end

    % Hide unused pins
    for k = (nUsed+1):S.maxPins
        set(S.pinSide(k), 'Visible','off');
        set(S.pinCap(k),  'Visible','off');
    end
end