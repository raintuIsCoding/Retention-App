%% Retention-end visualizer with LIVE UI + live geometry updates
clear; clc; close all;

%% ---------- DEFAULT INPUTS ----------
S.ID          = 8.0;      % inner diameter of casing (in)
S.t           = 0.25;     % wall thickness (in)
S.L_casing    = 10.0;     % length of modeled region (in)

S.nRows       = 3;        % axial rows of pins
S.nPinsPerRow = 12;       % pins around circumference per row

S.rowSpacing  = 0.75;     % axial spacing between rows (in)
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

mkLabel('MEOP (psi)', yTop-3*dy);              S.edMEOP = mkEdit(S.MEOP_psi, yTop-3*dy, 'MEOP_psi');
mkLabel('Design Factor (xMEOP)', yTop-4*dy);   S.edDF   = mkEdit(S.DF,       yTop-4*dy, 'DF');

S.MEOP_psi = str2double(get(S.edMEOP,'String'));
S.DF       = str2double(get(S.edDF,'String'));

if ~isfinite(S.MEOP_psi) || S.MEOP_psi <= 0, S.MEOP_psi = 900; end
if ~isfinite(S.DF)       || S.DF <= 0,       S.DF = 1.5; end

set(S.edMEOP,'String', num2str(S.MEOP_psi));
set(S.edDF,'String',   num2str(S.DF));

%% ---------- RENDER CONSTANTS ----------
S.nCirc    = 150;   % casing mesh resolution
S.nPinCirc = 30;   % pin mesh resolution
S.maxPins  = 400;  % pool size (increase if you expect more than this)

%% ---------- FIGURE / AXES ----------
S.fig = figure('Renderer','opengl');
S.ax  = axes('Parent', S.fig);
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
x2    = [0, 1]; % placeholder; will be set in updateCasing()
[Theta, X] = meshgrid(theta, x2);

S.hCasingOuter = surf(S.ax, X, X, X, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','none');
S.hCasingInner = surf(S.ax, X, X, X, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none');
S.hCap0        = surf(S.ax, X, X, X, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none');
S.hCapL        = surf(S.ax, X, X, X, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none');

%% ---------- CREATE PIN GRAPHICS POOL (handles once) ----------
% Each pin has: a side surf + a cap patch (fill3)
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

%% ---------- NESTED CALLBACK ----------
function onAnyInputChanged(src, ~)
    S = guidata(src); %#ok<NASGU>
    S = readControlsToState(S);
    S = updateAll(S);
    guidata(S.fig, S);
end

%% ---------- BUILD HUD / CONTROLS ----------
function S = buildHUDandControls(S)
    fig = S.fig;

    % Main container panel (top-left; extended downward)
    S.hud = uipanel(fig, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.05 0.42 0.90], ...
        'Title', 'Metrics', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'BorderType', 'line', ...
        'HighlightColor', [0 0 0]);

    % Inputs panel (now interactive controls)
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
        'Position', [0.02 0.52 0.96 0.46], ...
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
        'Position', [0.02 0.02 0.96 0.46], ...
        'Title', 'Force / Stress', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);

    S.txtStress = uicontrol(S.pOutStress, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.02 0.96 0.96], ...
        'String', '(none yet)', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'FontName', 'Consolas', ...
        'FontSize', 10);

    % ---- Controls layout inside Inputs panel ----
    % A simple 2-column form: label left, control right
    % Normalized positions
    yTop = 0.88; dy = 0.085;
    xL  = 0.04; wL = 0.56;
    xC  = 0.62; wC = 0.34;
    h   = 0.09;

    % Helper to create label
    mkLabel = @(txt, yy) uicontrol(S.pIn, 'Style','text', 'Units','normalized', ...
        'Position',[xL yy wL h], 'String',txt, 'HorizontalAlignment','left', ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], 'FontName','Consolas', 'FontSize',9);

    % Helper to create edit
    mkEdit = @(val, yy, tag) uicontrol(S.pIn, 'Style','edit', 'Units','normalized', ...
        'Position',[xC yy wC h], 'String',num2str(val), 'Tag',tag, ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], ...
        'Callback',@onAnyInputChanged);

    % Build controls
    mkLabel('Inner Diameter (in)', yTop);          S.edID        = mkEdit(S.ID,        yTop,      'ID');
    mkLabel('Wall Thickness (in)', yTop-dy);       S.edT         = mkEdit(S.t,         yTop-dy,   't');
    mkLabel('Modeled Length (in)', yTop-2*dy);     S.edL         = mkEdit(S.L_casing,  yTop-2*dy, 'L');

    mkLabel('Axial Rows', yTop-3*dy);              S.edRows      = mkEdit(S.nRows,     yTop-3*dy, 'nRows');
    mkLabel('Pins per Row', yTop-4*dy);            S.edPPR       = mkEdit(S.nPinsPerRow,yTop-4*dy,'nPinsPerRow');

    mkLabel('Row Spacing (in)', yTop-5*dy);        S.edRowSp     = mkEdit(S.rowSpacing,yTop-5*dy,'rowSpacing');
    mkLabel('First Row Offset from End (in)', yTop-6*dy); S.edFirst = mkEdit(S.firstRowZ,yTop-6*dy,'firstRowZ');

    mkLabel('Pin Diameter (in)', yTop-7*dy);       S.edPinDia    = mkEdit(S.pinDia,    yTop-7*dy,'pinDia');

    % Pin length is tied to thickness by default; still editable
    mkLabel('Pin Engagement Length (in)', yTop-8*dy); S.edPinLen  = mkEdit(S.pinLen,   yTop-8*dy,'pinLen');

    % Pattern toggle button
    mkLabel('Pin Pattern', yTop-9*dy);
    
    isAlt = (S.pinPatternMode == "alternating");
    btnStr = "Pattern: Progressive";
    if isAlt, btnStr = "Pattern: Alternating"; end
    
    S.btnPattern = uicontrol(S.pIn, 'Style','togglebutton', 'Units','normalized', ...
        'Position',[xC (yTop-9*dy) wC h], ...
        'String', btnStr, ...
        'Value', isAlt, ... % 0 = progressive, 1 = alternating
        'BackgroundColor',[0.95 0.95 0.95], ...
        'ForegroundColor',[0 0 0], ...
        'FontName','Consolas', ...
        'FontSize', 9, ...
        'Callback',@onAnyInputChanged);

    % Store state in figure
    guidata(S.fig, S);
end

%% ---------- READ CONTROLS -> STATE ----------
function S = readControlsToState(S)
    % Numeric reads (with basic safety)
    S.ID        = str2double(get(S.edID,'String'));
    S.t         = str2double(get(S.edT,'String'));
    S.L_casing  = str2double(get(S.edL,'String'));
    S.nRows     = max(1, round(str2double(get(S.edRows,'String'))));
    S.nPinsPerRow = max(1, round(str2double(get(S.edPPR,'String'))));
    S.rowSpacing  = str2double(get(S.edRowSp,'String'));
    S.firstRowZ   = str2double(get(S.edFirst,'String'));
    S.pinDia      = str2double(get(S.edPinDia,'String'));
    S.pinLen      = str2double(get(S.edPinLen,'String'));

    % Pattern (toggle button)
    if get(S.btnPattern, 'Value') == 1
        S.pinPatternMode = "alternating";
    else
        S.pinPatternMode = "progressive";
    end
    
    % Keep the button label synced
    if S.pinPatternMode == "alternating"
        set(S.btnPattern, 'String', "Pattern: Alternating");
    else
        set(S.btnPattern, 'String', "Pattern: Progressive");
    end

    % Keep a few derived constraints consistent (optional)
    S.edgeMarginEnd = S.t;

    % Clamp obvious bad values (avoid NaNs blowing up graphics)
    if ~isfinite(S.ID) || S.ID<=0, S.ID = 8; end
    if ~isfinite(S.t)  || S.t<=0,  S.t = 0.25; end
    if ~isfinite(S.L_casing) || S.L_casing<=0, S.L_casing = 10; end
    if ~isfinite(S.rowSpacing) || S.rowSpacing<=0, S.rowSpacing = 0.75; end
    if ~isfinite(S.firstRowZ) || S.firstRowZ<0, S.firstRowZ = 0.75; end
    if ~isfinite(S.pinDia) || S.pinDia<=0, S.pinDia = 0.375; end
    if ~isfinite(S.pinLen) || S.pinLen<=0, S.pinLen = 2*S.t; end

    % Push back integer formatting to UI
    set(S.edRows,'String', num2str(S.nRows));
    set(S.edPPR,'String', num2str(S.nPinsPerRow));
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

    % Update casing geometry
    S = updateCasingSurfaces(S, r_i, r_o);

    % Update pin geometry (side + cap)
    S = updatePins(S, r_o);

    % Update outputs text
    totalPins = S.nRows * S.nPinsPerRow;

    % ----- Load calcs -----
    p_design = S.MEOP_psi * S.DF;          % psi
    A_bore   = pi * (S.ID/2)^2;            % in^2 (pressure area)
    F_axial  = p_design * A_bore;          % lbf  (since psi * in^2 = lbf)
    
    totalPins = S.nRows * S.nPinsPerRow;
    F_perPin  = F_axial / totalPins;       % lbf average

    geomOutText = sprintf([ ...
        'Total Pins: %d\n' ...
        'Last Row X: %.2f in\n' ...
        'Max Pins/Row (min pitch): %d\n' ...
        'Packing OK: %s\n' ...
        'Axial OK: %s' ], ...
        totalPins, lastRowX, maxPinsCirc, tfStr(packingOK), tfStr(axialOK));

    set(S.txtGeom, 'String', geomOutText);

    stressOutText = sprintf([ ...
        'Design Pressure: %.1f psi\n' ...
        'Axial End Force: %.0f lbf\n' ...
        'Avg Force per Pin: %.1f lbf' ], ...
        p_design, F_axial, F_perPin);
    
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

    % Outer surface
    Yo = r_o*cos(Theta);
    Zo = r_o*sin(Theta);

    % Inner surface
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
        % Hard stop
        warning('Pin pool too small: need %d pins but maxPins=%d. Increase S.maxPins.', nUsed, S.maxPins);
        nUsed = S.maxPins;
    end

    % Local pin mesh
    phi = linspace(0, 2*pi, S.nPinCirc);
    r_pin = S.pinDia/2;

    % Cylinder along +Y
    Yloc_side = [r_o*ones(1,S.nPinCirc); (r_o+S.pinLen)*ones(1,S.nPinCirc)];
    Xloc_side = [r_pin*cos(phi);         r_pin*cos(phi)];
    Zloc_side = [r_pin*sin(phi);         r_pin*sin(phi)];

    % Tip cap at Y=r_o+pinLen
    Yloc_cap = (r_o + S.pinLen)*ones(1,S.nPinCirc);
    Xloc_cap = r_pin*cos(phi);
    Zloc_cap = r_pin*sin(phi);

    pitch = 2*pi / S.nPinsPerRow;

    idx = 0;
    for iRow = 1:S.nRows
        xRow = S.firstRowZ + (iRow - 1)*S.rowSpacing;

        % row offset
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
