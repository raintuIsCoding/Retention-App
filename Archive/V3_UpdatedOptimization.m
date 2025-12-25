%% Retention-end visualizer with LIVE UI + live geometry updates + OPTIMIZATION
% Notes:
% - MEOP + Design Factor are added as live inputs
% - Force/Stress outputs now show Design Pressure, Axial End Force, Avg Force per Pin
% - Pin pattern toggle button included (progressive vs alternating)
% - Geometry and pins update by updating existing graphics objects (fast)
% - OPTIMIZATION SECTION ADDED: Finds optimal parameters to meet stress constraints
%   while minimizing weight and pin extension into casing
clear; clc; close all;
%% ---------- DEFAULT INPUTS ----------
S.ID          = 8.0;      % inner diameter of casing (in)
S.t           = 0.25;     % wall thickness (in)
S.L_casing    = 10.0;     % length of modeled region (in)
S.nRows       = 3;        % axial rows of pins
S.nPinsPerRow = 12;       % pins around circumference per row
S.rowSpacing  = 0.5;      % axial spacing between rows (in)
S.firstRowZ   = 0.75;     % axial position of first row from x=0 (in)
S.pinDia      = 0.375;    % pin diameter (in)
S.pinLen      = 2 * S.t;  % how far pins stick out radially (in)
% Loads inputs (FIXED for this design)
S.MEOP_psi   = 850;   % psi
S.DF         = 1.5;   % design factor (MEOP multiplier)
% Geometric constraints
S.minCircPitchFactor = 2.0;  % min center spacing = factor * pinDia (circumferential)
S.edgeMarginEnd      = S.t;  % safety margin from last row to casing end (approx)
% Pattern toggle (per your definitions)
S.pinPatternMode = "progressive";   % "progressive" or "alternating"
S.altStartPhase  = 0;               % 0 or 1
%% ---------- TARGET STRESS LIMITS (from previous successful rocket) ----------
% These are the values from your previous rocket that worked well
% Goal: achieve these numbers or lower while minimizing total weight
S.targets.shearOut_max    = 1.837831702;    % KSI
S.targets.netTension_max  = 11.50766592;    % KSI
S.targets.pinShear_max    = 18.72;          % KSI
S.targets.bearing_max     = 18.37831702;    % KSI
S.targets.hoop_max        = 15.1125;        % KSI
S.targets.axial_max       = 7.55625;        % KSI
S.targets.pinShearFOS_min = 0.323717949;    % Minimum FOS (must be >= this)
%% ---------- CASING PARAMETERS ----------
S.casingLength = 96;              % Total casing length (in) - for weight calc
S.density_CF = 0.0535;            % Carbon fiber density (lb/in³) - Based off of MIII casing measurement
S.density_Al = 0.0975;            % 6061 Aluminum density (lb/in³)
%% ---------- GEOMETRIC CONSTRAINTS ----------
S.minCircPitchFactor = 3.0;         % min center spacing = factor * pinDia (INCREASED for safety)
S.minAxialPitchFactor = 2.5;        % min axial spacing = factor * pinDia
S.edgeMarginEnd = S.t;              % safety margin from last row to casing end
%% ---------- OPTIMIZATION BOUNDS ----------
% Define the allowable ranges for optimization variables
% Format: [min, max]
S.optBounds.t           = [0.2, 0.75];      % Wall thickness (in) - need thicker for hoop stress
S.optBounds.nRows       = [1, 5];           % Number of axial rows
S.optBounds.nPinsPerRow = [6, 30];          % Pins per row (limited by pitch factor)
S.optBounds.rowSpacing  = [0.5, 2.0];       % Row spacing (in)
S.optBounds.firstRowZ   = [0.75, 2.5];      % First row offset (in)
S.optBounds.pinDia      = [0.25, 0.5];      % Pin diameter (in)
%% ---------- RENDER CONSTANTS ----------
S.nCirc    = 150;  % casing mesh resolution
S.nPinCirc = 30;   % pin mesh resolution
S.maxPins  = 400;  % pool size (increase if you expect more than this)
%% ---------- FIGURE / AXES ----------
S.fig = figure('Renderer','opengl', 'Position', [100, 100, 1400, 800]);
S.ax = axes( ...
    'Parent', S.fig, ...
    'Units', 'normalized', ...
    'Position', [0.55 0.10 0.40 0.80]);  % shifted right for optimization panel
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
%% ---------- CONSTRAINT CHANGED CALLBACK ----------
function onConstraintChanged(src, ~)
    S = guidata(src);
    S = readConstraintsFromUI(S);
    S = updateAll(S);
    guidata(S.fig, S);
end
%% ---------- BOUNDS CHANGED CALLBACK ----------
function onBoundsChanged(src, ~)
    S = guidata(src);
    S = readBoundsFromUI(S);
    S = readConstraintsFromUI(S);
    guidata(S.fig, S);
end
%% ---------- READ CONSTRAINTS FROM UI ----------
function S = readConstraintsFromUI(S)
    S.targets.shearOut_max    = str2double(get(S.edTgtShearOut, 'String'));
    S.targets.netTension_max  = str2double(get(S.edTgtNetTens, 'String'));
    S.targets.pinShear_max    = str2double(get(S.edTgtPinShear, 'String'));
    S.targets.bearing_max     = str2double(get(S.edTgtBearing, 'String'));
    S.targets.hoop_max        = str2double(get(S.edTgtHoop, 'String'));
    S.targets.axial_max       = str2double(get(S.edTgtAxial, 'String'));
    S.targets.pinShearFOS_min = str2double(get(S.edTgtFOS, 'String'));
end
%% ---------- READ BOUNDS FROM UI ----------
function S = readBoundsFromUI(S)
    % Wall thickness bounds
    tMin = str2double(get(S.edBndTmin, 'String'));
    tMax = str2double(get(S.edBndTmax, 'String'));
    % Changed condition to >= to allow fixing values (Max == Min)
    if isfinite(tMin) && isfinite(tMax) && tMin > 0 && tMax >= tMin
        S.optBounds.t = [tMin, tMax];
    end
    
    % Axial rows bounds
    rowsMin = round(str2double(get(S.edBndRowsMin, 'String')));
    rowsMax = round(str2double(get(S.edBndRowsMax, 'String')));
    % Changed condition to >= 
    if isfinite(rowsMin) && isfinite(rowsMax) && rowsMin >= 1 && rowsMax >= rowsMin
        S.optBounds.nRows = [rowsMin, rowsMax];
    end
    
    % Pins per row bounds
    pprMin = round(str2double(get(S.edBndPPRmin, 'String')));
    pprMax = round(str2double(get(S.edBndPPRmax, 'String')));
    % Changed condition to >=
    if isfinite(pprMin) && isfinite(pprMax) && pprMin >= 3 && pprMax >= pprMin
        S.optBounds.nPinsPerRow = [pprMin, pprMax];
    end
    
    % Row spacing bounds
    rsMin = str2double(get(S.edBndRowSpMin, 'String'));
    rsMax = str2double(get(S.edBndRowSpMax, 'String'));
    % Changed condition to >=
    if isfinite(rsMin) && isfinite(rsMax) && rsMin > 0 && rsMax >= rsMin
        S.optBounds.rowSpacing = [rsMin, rsMax];
    end
    
    % First row offset bounds
    frMin = str2double(get(S.edBndFirstMin, 'String'));
    frMax = str2double(get(S.edBndFirstMax, 'String'));
    % Changed condition to >=
    if isfinite(frMin) && isfinite(frMax) && frMin > 0 && frMax >= frMin
        S.optBounds.firstRowZ = [frMin, frMax];
    end
    
    % Pin diameter bounds
    pdMin = str2double(get(S.edBndPinDiaMin, 'String'));
    pdMax = str2double(get(S.edBndPinDiaMax, 'String'));
    % Changed condition to >=
    if isfinite(pdMin) && isfinite(pdMax) && pdMin > 0 && pdMax >= pdMin
        S.optBounds.pinDia = [pdMin, pdMax];
    end
    
    % Min pitch factor
    pitchFactor = str2double(get(S.edMinCircPitch, 'String'));
    if isfinite(pitchFactor) && pitchFactor >= 2
        S.minCircPitchFactor = pitchFactor;
    end
end
%% ---------- OPTIMIZATION CALLBACK ----------
function onOptimizeClicked(src, ~)
    S = guidata(src);
    S = readControlsToState(S);
    S = readConstraintsFromUI(S);
    S = readBoundsFromUI(S);
    
    % Update status
    set(S.txtOptStatus, 'String', 'Optimizing... Please wait.');
    drawnow;
    
    % Run optimization
    [S, success] = runOptimization(S);
    
    if success
        % Update UI with optimized values
        set(S.edT, 'String', num2str(S.t, '%.4f'));
        set(S.edRows, 'String', num2str(S.nRows));
        set(S.edPPR, 'String', num2str(S.nPinsPerRow));
        set(S.edRowSp, 'String', num2str(S.rowSpacing, '%.4f'));
        set(S.edFirst, 'String', num2str(S.firstRowZ, '%.4f'));
        set(S.edPinDia, 'String', num2str(S.pinDia, '%.4f'));
        set(S.edPinLen, 'String', num2str(S.pinLen, '%.4f'));
        
        set(S.txtOptStatus, 'String', 'Optimization Complete!');
    else
        set(S.txtOptStatus, 'String', 'Optimization Failed - Try adjusting bounds');
    end
    
    S = updateAll(S);
    guidata(S.fig, S);
end
%% ---------- RUN OPTIMIZATION ----------
function [S, success] = runOptimization(S)
    success = false;
    
    % Fixed parameters
    ID = S.ID;
    MEOP_psi = S.MEOP_psi;
    DF = S.DF;
    minCircPitchFactor = S.minCircPitchFactor;
    minAxialPitchFactor = S.minAxialPitchFactor;
    L_casing = S.L_casing;
    casingLength = S.casingLength;
    density_CF = S.density_CF;
    density_Al = S.density_Al;
    targets = S.targets;
    bounds = S.optBounds;
    
    % Optimization variables: [t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinDia]
    
    % Initial guess (current values)
    x0 = [S.t, S.nRows, S.nPinsPerRow, S.rowSpacing, S.firstRowZ, S.pinDia];
    
    % Bounds
    lb = [bounds.t(1), bounds.nRows(1), bounds.nPinsPerRow(1), ...
          bounds.rowSpacing(1), bounds.firstRowZ(1), bounds.pinDia(1)];
    ub = [bounds.t(2), bounds.nRows(2), bounds.nPinsPerRow(2), ...
          bounds.rowSpacing(2), bounds.firstRowZ(2), bounds.pinDia(2)];
    
    % Integer constraints for nRows and nPinsPerRow
    intcon = [2, 3];
    
    % Objective function: Minimize total system weight
    objFun = @(x) objectiveFunction(x, ID, MEOP_psi, DF, casingLength, density_CF, density_Al);
    
    % Nonlinear constraints
    nonlconFun = @(x) nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, ...
                                           minCircPitchFactor, minAxialPitchFactor, targets);
    
    % Optimization options
    options = optimoptions('ga', ...
        'Display', 'iter', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 200, ...
        'FunctionTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6, ...
        'UseParallel', false, ...
        'PlotFcn', []);
    
    try
        % Use genetic algorithm for mixed-integer optimization
        [xOpt, fval, exitflag] = ga(objFun, 6, [], [], [], [], lb, ub, nonlconFun, intcon, options);
        
        if exitflag > 0
            % Extract optimized values
            S.t = xOpt(1);
            S.nRows = round(xOpt(2));
            S.nPinsPerRow = round(xOpt(3));
            S.rowSpacing = xOpt(4);
            S.firstRowZ = xOpt(5);
            S.pinDia = xOpt(6);
            S.pinLen = 2 * S.t;  % Pin length = 2x wall thickness
            S.edgeMarginEnd = S.t;
            
            success = true;
            fprintf('\n=== OPTIMIZATION RESULTS ===\n');
            fprintf('Wall Thickness: %.4f in\n', S.t);
            fprintf('Axial Rows: %d\n', S.nRows);
            fprintf('Pins per Row: %d\n', S.nPinsPerRow);
            fprintf('Row Spacing: %.4f in\n', S.rowSpacing);
            fprintf('First Row Offset: %.4f in\n', S.firstRowZ);
            fprintf('Pin Diameter: %.4f in\n', S.pinDia);
            fprintf('Objective Value: %.4f\n', fval);
        else
            fprintf('Optimization did not converge. Exit flag: %d\n', exitflag);
            % Try grid search as fallback
            [S, success] = gridSearchOptimization(S);
        end
    catch ME
        fprintf('Optimization error: %s\n', ME.message);
        fprintf('Attempting grid search fallback...\n');
        [S, success] = gridSearchOptimization(S);
    end
end
%% ---------- GRID SEARCH OPTIMIZATION (Fallback) ----------
function [S, success] = gridSearchOptimization(S)
    success = false;
    
    fprintf('\n=== RUNNING GRID SEARCH OPTIMIZATION ===\n');
    
    % Fixed parameters
    ID = S.ID;
    MEOP_psi = S.MEOP_psi;
    DF = S.DF;
    minCircPitchFactor = S.minCircPitchFactor;
    minAxialPitchFactor = S.minAxialPitchFactor;
    L_casing = S.L_casing;
    casingLength = S.casingLength;
    density_CF = S.density_CF;
    density_Al = S.density_Al;
    targets = S.targets;
    bounds = S.optBounds;
    
    % Calculate approximate number of pins needed
    p_design = MEOP_psi * DF;
    A_bore = pi * (ID/2)^2;
    F_axial = p_design * A_bore;
    
    fprintf('Design axial force: %.0f lbf\n', F_axial);
    
    % Smart grid
    if bounds.t(2) > bounds.t(1)
        t_vals = linspace(bounds.t(1), bounds.t(2), 8);
    else
        t_vals = bounds.t(1);
    end
    
    nRows_vals = bounds.nRows(1):bounds.nRows(2);
    
    if bounds.nPinsPerRow(2) > bounds.nPinsPerRow(1)
        nPinsPerRow_vals = bounds.nPinsPerRow(1):2:bounds.nPinsPerRow(2);
    else
        nPinsPerRow_vals = bounds.nPinsPerRow(1);
    end
    
    if bounds.rowSpacing(2) > bounds.rowSpacing(1)
        rowSpacing_vals = linspace(bounds.rowSpacing(1), bounds.rowSpacing(2), 6);
    else
        rowSpacing_vals = bounds.rowSpacing(1);
    end
    
    if bounds.firstRowZ(2) > bounds.firstRowZ(1)
        firstRowZ_vals = linspace(bounds.firstRowZ(1), bounds.firstRowZ(2), 5);
    else
        firstRowZ_vals = bounds.firstRowZ(1);
    end
    
    if bounds.pinDia(2) > bounds.pinDia(1)
        pinDia_vals = linspace(bounds.pinDia(1), bounds.pinDia(2), 6);
    else
        pinDia_vals = bounds.pinDia(1);
    end
    
    bestObj = inf;
    bestParams = [];
    nEval = 0;
    nFeasible = 0;
    
    totalCombinations = length(t_vals) * length(nRows_vals) * length(nPinsPerRow_vals) * ...
                        length(rowSpacing_vals) * length(firstRowZ_vals) * length(pinDia_vals);
    fprintf('Total combinations to evaluate: %d\n', totalCombinations);
    
    % Store all feasible solutions for reporting
    feasibleSolutions = [];
    
    for t = t_vals
        for nRows = nRows_vals
            for nPinsPerRow = nPinsPerRow_vals
                for rowSpacing = rowSpacing_vals
                    for firstRowZ = firstRowZ_vals
                        for pinDia = pinDia_vals
                            nEval = nEval + 1;
                            
                            % Check constraints
                            [c, ceq] = nonlinearConstraints([t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinDia], ...
                                                            ID, MEOP_psi, DF, L_casing, minCircPitchFactor, minAxialPitchFactor, targets);
                            
                            % All constraints must be <= 0
                            if all(c <= 0)
                                nFeasible = nFeasible + 1;
                                
                                % Calculate objective (total system weight)
                                obj = objectiveFunction([t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinDia], ...
                                                        ID, MEOP_psi, DF, casingLength, density_CF, density_Al);
                                
                                % Store this solution
                                sol.t = t;
                                sol.nRows = nRows;
                                sol.nPinsPerRow = nPinsPerRow;
                                sol.rowSpacing = rowSpacing;
                                sol.firstRowZ = firstRowZ;
                                sol.pinDia = pinDia;
                                sol.obj = obj;
                                sol.constraints = c;
                                
                                if isempty(feasibleSolutions)
                                    feasibleSolutions = sol;
                                else
                                    feasibleSolutions(end+1) = sol;
                                end
                                
                                if obj < bestObj
                                    bestObj = obj;
                                    bestParams = [t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinDia];
                                end
                            end
                        end
                    end
                end
            end
        end
        if mod(nEval, 100) == 0
            fprintf('Progress: %.1f%% complete, %d feasible solutions found\n', ...
                    100*nEval/totalCombinations, nFeasible);
        end
    end
    
    if ~isempty(bestParams)
        S.t = bestParams(1);
        S.nRows = round(bestParams(2));
        S.nPinsPerRow = round(bestParams(3));
        S.rowSpacing = bestParams(4);
        S.firstRowZ = bestParams(5);
        S.pinDia = bestParams(6);
        S.pinLen = 2 * S.t;
        S.edgeMarginEnd = S.t;
        
        success = true;
        fprintf('\n=== GRID SEARCH RESULTS ===\n');
        fprintf('Wall Thickness: %.4f in\n', S.t);
        fprintf('Axial Rows: %d\n', S.nRows);
        fprintf('Pins per Row: %d\n', S.nPinsPerRow);
        fprintf('Row Spacing: %.4f in\n', S.rowSpacing);
        fprintf('First Row Offset: %.4f in\n', S.firstRowZ);
        fprintf('Pin Diameter: %.4f in\n', S.pinDia);
        fprintf('Total Pins: %d\n', S.nRows * S.nPinsPerRow);
        fprintf('Objective Value: %.4f\n', bestObj);
        fprintf('Feasible solutions found: %d out of %d\n', nFeasible, nEval);
        
        % Show top 3 solutions if available
        if length(feasibleSolutions) > 1
            [~, sortIdx] = sort([feasibleSolutions.obj]);
            fprintf('\n--- Top 3 Solutions ---\n');
            for i = 1:min(3, length(sortIdx))
                sol = feasibleSolutions(sortIdx(i));
                fprintf('%d) t=%.3f, rows=%d, pins/row=%d, pinDia=%.3f, obj=%.4f\n', ...
                    i, sol.t, sol.nRows, sol.nPinsPerRow, sol.pinDia, sol.obj);
            end
        end
    else
        fprintf('No feasible solution found in grid search.\n');
        fprintf('Consider relaxing constraints or expanding bounds.\n');
        fprintf('\nDiagnostics:\n');
        fprintf('- Axial force: %.0f lbf (this is %.1fx higher than previous rocket)\n', ...
                F_axial, F_axial/18375);
        fprintf('- Min pins needed (est): %.0f\n', F_axial / (targets.pinShear_max * 1000 * pi/4 * 0.375^2));
    end
end
%% ---------- OBJECTIVE FUNCTION ----------
function obj = objectiveFunction(x, ID, MEOP_psi, DF, casingLength, density_CF, density_Al)
    % Extract variables
    t = x(1);
    nRows = round(x(2));
    nPinsPerRow = round(x(3));
    rowSpacing = x(4);
    firstRowZ = x(5);
    pinDia = x(6);
    
    % ===== TOTAL SYSTEM WEIGHT =====
    
    % 1. Carbon Fiber Casing Mass (80" long tube)
    OD = ID + 2*t;
    Volume_casing = (pi/4) * (OD^2 - ID^2) * casingLength;
    Mass_casing = Volume_casing * density_CF;
    
    % 2. Retention Ring Mass (aluminum, both ends)
    retention_ring_thickness = 0.25; % in (mandrel thickness)
    edgeMarginEnd = t;
    retention_length = edgeMarginEnd*2 + rowSpacing*(nRows-1);
    Volume_retention_ring = (pi/4) * (ID^2 - (ID - 2*retention_ring_thickness)^2) * retention_length;
    % Subtract pin holes from retention ring
    Volume_pins_in_ring = nPinsPerRow * nRows * pi * (pinDia/2)^2 * retention_ring_thickness;
    Volume_retention_net = Volume_retention_ring - Volume_pins_in_ring;
    Mass_retention_rings = Volume_retention_net * density_Al * 2; % x2 for both ends
    
    % 3. Pin Mass (aluminum cylinders)
    pinLen = 2 * t;  % Pins extend 2x wall thickness
    totalPins = nRows * nPinsPerRow * 2;  % x2 for both ends
    Volume_all_pins = totalPins * pi * (pinDia/2)^2 * pinLen;
    Mass_pins = Volume_all_pins * density_Al;
    
    % Total system mass
    totalMass = Mass_casing + Mass_retention_rings + Mass_pins;
    
    % Also penalize pin extension (last row position)
    lastRowX = firstRowZ + (nRows - 1) * rowSpacing;
    
    % Objective: minimize total mass with small penalty for extension
    % Normalize: typical casing ~15-25 lb, extension ~2-5 in
    w_mass = 1.0;
    w_extension = 0.1;  % Small weight since mass is primary concern
    
    obj = w_mass * (totalMass / 20.0) + w_extension * (lastRowX / 5.0);
end
%% ---------- NONLINEAR CONSTRAINTS ----------
function [c, ceq] = nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, minCircPitchFactor, minAxialPitchFactor, targets)
    % Extract variables
    t = x(1);
    nRows = round(x(2));
    nPinsPerRow = round(x(3));
    rowSpacing = x(4);
    firstRowZ = x(5);
    pinDia = x(6);
    
    % Derived geometry
    r_i = ID / 2;
    r_o = r_i + t;
    edgeMarginEnd = t;
    
    % Total pins
    totalPins = nRows * nPinsPerRow;
    
    % -------- GEOMETRIC CONSTRAINTS --------
    % Circumferential packing check
    r_center = r_o;
    circumference = 2 * pi * r_center;
    minCircPitch = minCircPitchFactor * pinDia;
    maxPinsCirc = floor(circumference / minCircPitch);
    c_packing = nPinsPerRow - maxPinsCirc;  % Must be <= 0
    
    % Axial packing check (row spacing >= min axial pitch)
    minAxialPitch = minAxialPitchFactor * pinDia;
    c_axialPitch = minAxialPitch - rowSpacing;  % Must be <= 0
    
    % Axial extent check
    lastRowX = firstRowZ + (nRows - 1) * rowSpacing;
    c_axial = (lastRowX + edgeMarginEnd) - L_casing;  % Must be <= 0
    
    % -------- STRESS CALCULATIONS --------
    n = 1:nRows;
    A_bore = pi * (ID/2)^2;
    A_ShearOut = sum((rowSpacing*(n-1) + firstRowZ) * t * 2 * nPinsPerRow);
    
    % Net tension area: cross-sectional area minus pin holes at ONE row
    A_wall_gross = pi * 0.25 * ((ID + 2*t)^2 - ID^2);
    A_pin_holes_per_row = pinDia * nPinsPerRow * t;
    A_tension = A_wall_gross - A_pin_holes_per_row;
    
    % Constraint: Net area must be positive (at least 30% of gross remaining)
    c_netArea = (0.3 * A_wall_gross) - A_tension;  % Must be <= 0
    
    % Safety: ensure positive area for stress calculations
    if A_tension <= 0, A_tension = 0.001; end
    if A_ShearOut <= 0, A_ShearOut = 0.001; end
    
    % Loads
    p_design = MEOP_psi * DF;
    F_axial = p_design * A_bore;
    F_perPin = F_axial / totalPins;
    
    % Stresses (KSI)
    Stress_ShearOut = F_axial / A_ShearOut / 1000;
    Net_Tension = (F_axial / A_tension) / 1000;
    Pin_Shear = (F_perPin / (pi * 0.25 * pinDia^2)) / 1000;
    Bearing = (F_perPin / (pinDia * t)) / 1000;
    Hoop = (p_design * (((ID + 2*t)^2 + ID^2) / ((ID + 2*t)^2 - ID^2))) / 1000;
    Axial_stress = Hoop / 2;
    
    % Factor of Safety
    Pin_Shear_Strength = 43.5; % KSI
    Pin_Shear_FOS = (Pin_Shear_Strength / Pin_Shear) - 2;
    
    % -------- STRESS CONSTRAINTS --------
    c_shearOut = Stress_ShearOut - targets.shearOut_max;
    c_netTension = Net_Tension - targets.netTension_max;
    c_pinShear = Pin_Shear - targets.pinShear_max;
    c_bearing = Bearing - targets.bearing_max;
    c_hoop = Hoop - targets.hoop_max;
    c_axialStress = Axial_stress - targets.axial_max;
    c_FOS = targets.pinShearFOS_min - Pin_Shear_FOS;
    
    % Combine all inequality constraints
    c = [c_packing; c_axialPitch; c_axial; c_netArea; ...
         c_shearOut; c_netTension; c_pinShear; c_bearing; c_hoop; c_axialStress; c_FOS];
    
    % No equality constraints
    ceq = [];
end
%% ---------- BUILD HUD / CONTROLS ----------
function S = buildHUDandControls(S)
    fig = S.fig;
    % Main container panel (top-left; extended downward, wider)
    S.hud = uipanel(fig, ...
        'Units', 'normalized', ...
        'Position', [0.01 0.02 0.52 0.96], ...
        'Title', 'Retention Ring Calculator', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'BorderType', 'line', ...
        'HighlightColor', [0 0 0]);
    % Inputs panel (interactive)
    S.pIn = uipanel(S.hud, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.52 0.47 0.46], ...
        'Title', 'Inputs', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);
    % Outputs container
    S.pOut = uipanel(S.hud, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.02 0.47 0.48], ...
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
        'FontSize', 9);
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
        'FontSize', 9);
    % ========== OPTIMIZATION PANEL ==========
    S.pOpt = uipanel(S.hud, ...
        'Units', 'normalized', ...
        'Position', [0.51 0.02 0.47 0.96], ...
        'Title', 'Optimization & Constraints', ...
        'ForegroundColor', [0 0 0.6], ...
        'BackgroundColor', [0.95 0.95 1]);
    
    % ----- Editable Target Constraints -----
    S.pOptTargets = uipanel(S.pOpt, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.62 0.96 0.36], ...
        'Title', 'Target Constraints', ...
        'ForegroundColor', [0 0.4 0], ...
        'BackgroundColor', [0.9 1 0.9]);
    
    % Layout for editable constraints (compact)
    yT = 0.85; dyT = 0.115;
    xLT = 0.02; wLT = 0.55;
    xCT = 0.58; wCT = 0.40;
    hT = 0.10;
    
    mkLabelT = @(txt, yy) uicontrol(S.pOptTargets, 'Style','text', 'Units','normalized', ...
        'Position',[xLT yy wLT hT], 'String',txt, 'HorizontalAlignment','left', ...
        'BackgroundColor',[0.9 1 0.9], 'ForegroundColor',[0 0 0], 'FontName','Consolas', 'FontSize',7);
    
    mkEditT = @(val, yy, fmt) uicontrol(S.pOptTargets, 'Style','edit', 'Units','normalized', ...
        'Position',[xCT yy wCT hT], 'String',num2str(val, fmt), ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], 'FontSize',8, ...
        'Callback',@onConstraintChanged);
    
    mkLabelT('Shear Out (KSI)', yT);       S.edTgtShearOut = mkEditT(S.targets.shearOut_max, yT, '%.4f');
    mkLabelT('Net Tension (KSI)', yT-dyT);   S.edTgtNetTens  = mkEditT(S.targets.netTension_max, yT-dyT, '%.4f');
    mkLabelT('Pin Shear (KSI)', yT-2*dyT);     S.edTgtPinShear = mkEditT(S.targets.pinShear_max, yT-2*dyT, '%.4f');
    mkLabelT('Bearing (KSI)', yT-3*dyT);       S.edTgtBearing  = mkEditT(S.targets.bearing_max, yT-3*dyT, '%.4f');
    mkLabelT('Hoop (KSI)', yT-4*dyT);          S.edTgtHoop     = mkEditT(S.targets.hoop_max, yT-4*dyT, '%.4f');
    mkLabelT('Axial (KSI)', yT-5*dyT);         S.edTgtAxial    = mkEditT(S.targets.axial_max, yT-5*dyT, '%.4f');
    mkLabelT('Pin FOS Min', yT-6*dyT);         S.edTgtFOS      = mkEditT(S.targets.pinShearFOS_min, yT-6*dyT, '%.4f');
    
    % ----- Editable Optimization Bounds -----
    S.pOptBounds = uipanel(S.pOpt, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.30 0.96 0.30], ...
        'Title', 'Optimization Bounds [min, max]', ...
        'ForegroundColor', [0.4 0 0.4], ...
        'BackgroundColor', [1 0.95 1]);
    
    % Layout for bounds (3 columns: label, min, max)
    yB = 0.82; dyB = 0.14;
    xLB = 0.02; wLB = 0.38;
    xMinB = 0.42; wMinB = 0.27;
    xMaxB = 0.71; wMaxB = 0.27;
    hB = 0.12;
    
    mkLabelB = @(txt, yy) uicontrol(S.pOptBounds, 'Style','text', 'Units','normalized', ...
        'Position',[xLB yy wLB hB], 'String',txt, 'HorizontalAlignment','left', ...
        'BackgroundColor',[1 0.95 1], 'ForegroundColor',[0 0 0], 'FontName','Consolas', 'FontSize',7);
    
    mkEditB = @(val, yy, xpos) uicontrol(S.pOptBounds, 'Style','edit', 'Units','normalized', ...
        'Position',[xpos yy wMinB hB], 'String',num2str(val, '%.3f'), ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], 'FontSize',8, ...
        'Callback',@onBoundsChanged);
    
    % Header row
    uicontrol(S.pOptBounds, 'Style','text', 'Units','normalized', ...
        'Position',[xMinB yB+dyB*0.3 wMinB hB*0.7], 'String','Min', 'HorizontalAlignment','center', ...
        'BackgroundColor',[1 0.95 1], 'ForegroundColor',[0.3 0 0.3], 'FontName','Consolas', 'FontSize',7, 'FontWeight','bold');
    uicontrol(S.pOptBounds, 'Style','text', 'Units','normalized', ...
        'Position',[xMaxB yB+dyB*0.3 wMaxB hB*0.7], 'String','Max', 'HorizontalAlignment','center', ...
        'BackgroundColor',[1 0.95 1], 'ForegroundColor',[0.3 0 0.3], 'FontName','Consolas', 'FontSize',7, 'FontWeight','bold');
    
    mkLabelB('Wall Thick (in)', yB);
    S.edBndTmin = mkEditB(S.optBounds.t(1), yB, xMinB);
    S.edBndTmax = mkEditB(S.optBounds.t(2), yB, xMaxB);
    
    mkLabelB('Axial Rows', yB-dyB);
    S.edBndRowsMin = mkEditB(S.optBounds.nRows(1), yB-dyB, xMinB);
    S.edBndRowsMax = mkEditB(S.optBounds.nRows(2), yB-dyB, xMaxB);
    
    mkLabelB('Pins/Row', yB-2*dyB);
    S.edBndPPRmin = mkEditB(S.optBounds.nPinsPerRow(1), yB-2*dyB, xMinB);
    S.edBndPPRmax = mkEditB(S.optBounds.nPinsPerRow(2), yB-2*dyB, xMaxB);
    
    mkLabelB('Row Space (in)', yB-3*dyB);
    S.edBndRowSpMin = mkEditB(S.optBounds.rowSpacing(1), yB-3*dyB, xMinB);
    S.edBndRowSpMax = mkEditB(S.optBounds.rowSpacing(2), yB-3*dyB, xMaxB);
    
    mkLabelB('1st Row Ofs (in)', yB-4*dyB);
    S.edBndFirstMin = mkEditB(S.optBounds.firstRowZ(1), yB-4*dyB, xMinB);
    S.edBndFirstMax = mkEditB(S.optBounds.firstRowZ(2), yB-4*dyB, xMaxB);
    
    mkLabelB('Pin Dia (in)', yB-5*dyB);
    S.edBndPinDiaMin = mkEditB(S.optBounds.pinDia(1), yB-5*dyB, xMinB);
    S.edBndPinDiaMax = mkEditB(S.optBounds.pinDia(2), yB-5*dyB, xMaxB);
    
    % Min pitch factor (single value, not a range)
    mkLabelB('Min Pitch (x dia)', yB-6*dyB);
    S.edMinCircPitch = uicontrol(S.pOptBounds, 'Style','edit', 'Units','normalized', ...
        'Position',[xMinB yB-6*dyB wMinB+wMaxB+0.02 hB], 'String',num2str(S.minCircPitchFactor, '%.1f'), ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], 'FontSize',8, ...
        'Callback',@onBoundsChanged);
    
    % Optimization status/results panel
    S.pOptResults = uipanel(S.pOpt, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.12 0.96 0.16], ...
        'Title', 'Constraint Status', ...
        'ForegroundColor', [0 0 0], ...
        'BackgroundColor', [1 1 1]);
    
    S.txtOptResults = uicontrol(S.pOptResults, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.02 0.96 0.96], ...
        'String', 'Adjust constraints/bounds, then Run Optimization', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0], ...
        'FontName', 'Consolas', ...
        'FontSize', 7);
    
    % Optimize button
    S.btnOptimize = uicontrol(S.pOpt, ...
        'Style', 'pushbutton', ...
        'Units', 'normalized', ...
        'Position', [0.1 0.02 0.8 0.08], ...
        'String', 'Run Optimization', ...
        'BackgroundColor', [0.2 0.6 0.2], ...
        'ForegroundColor', [1 1 1], ...
        'FontName', 'Consolas', ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'Callback', @onOptimizeClicked);
    
    S.txtOptStatus = uicontrol(S.pOpt, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.10 0.96 0.02], ...
        'String', '', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.95 0.95 1], ...
        'ForegroundColor', [0 0 0.6], ...
        'FontName', 'Consolas', ...
        'FontSize', 8);
    % ---- Controls layout inside Inputs panel ----
    % Normalized positions
    yTop = 0.90; dy = 0.075;
    xL  = 0.04; wL = 0.52;
    xC  = 0.58; wC = 0.38;
    h   = 0.07;
    % Helper to create label
    mkLabel = @(txt, yy) uicontrol(S.pIn, 'Style','text', 'Units','normalized', ...
        'Position',[xL yy wL h], 'String',txt, 'HorizontalAlignment','left', ...
        'BackgroundColor',[1 1 1], 'ForegroundColor',[0 0 0], 'FontName','Consolas', 'FontSize',8);
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
    mkLabel('First Row Offset (in)', yTop-8*dy);           S.edFirst   = mkEdit(S.firstRowZ,  yTop-8*dy,   'firstRowZ');
    mkLabel('Pin Diameter (in)', yTop-9*dy);               S.edPinDia  = mkEdit(S.pinDia,     yTop-9*dy,   'pinDia');
    mkLabel('Pin Engage Length (in)', yTop-10*dy);         S.edPinLen  = mkEdit(S.pinLen,     yTop-10*dy,  'pinLen');
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
        'FontSize', 8, ...
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
    
    % Also read constraint and bounds values from UI
    S = readConstraintsFromUI(S);
    S = readBoundsFromUI(S);
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
    % ===== TOTAL SYSTEM WEIGHT CALCULATION =====
    retention_ring_thickness = 0.25; % Set by mandrel thickness
    
    % 1. Carbon Fiber Casing Mass
    OD = S.ID + 2*S.t;
    Volume_casing = (pi/4) * (OD^2 - S.ID^2) * S.casingLength;
    Mass_casing = Volume_casing * S.density_CF;
    
    % 2. Retention Ring Mass (both ends)
    totalPins = S.nRows * S.nPinsPerRow;
    retention_length = S.edgeMarginEnd*2 + S.rowSpacing*(S.nRows-1);
    Volume_retention_ring = (pi/4) * (S.ID^2 - (S.ID - 2*retention_ring_thickness)^2) * retention_length;
    Volume_pins_in_ring = totalPins * pi * (S.pinDia/2)^2 * retention_ring_thickness;
    Mass_retention_rings = (Volume_retention_ring - Volume_pins_in_ring) * S.density_Al * 2;
    
    % 3. Pin Mass (both ends)
    pinLen = 2 * S.t;
    Volume_all_pins = totalPins * 2 * pi * (S.pinDia/2)^2 * pinLen;
    Mass_pins = Volume_all_pins * S.density_Al;
    
    % Total system mass
    Total_Mass = Mass_casing + Mass_retention_rings + Mass_pins;
    geomOutText = sprintf([ ...
        'Total Pins: %d\n' ...
        'Last Row X: %.2f in\n' ...
        'Max Pins/Row: %d\n' ...
        'Packing OK: %s\n' ...
        'Axial OK: %s\n' ...
        '-------------------\n' ...
        'Casing Mass: %.2f lb\n' ...
        'Retention Mass: %.2f lb\n' ...
        'Pin Mass: %.2f lb\n' ...
        'TOTAL MASS: %.2f lb\n' ], ...
        totalPins, lastRowX, maxPinsCirc, tfStr(packingOK), tfStr(axialOK), ...
        Mass_casing, Mass_retention_rings, Mass_pins, Total_Mass);
    set(S.txtGeom, 'String', geomOutText);
    %%%%%%%%%%%%%%%% -------- AREA -------- %%%%%%%%%%%%%%%%
    n = 1:S.nRows;
    A_bore   = pi * (S.ID/2)^2;       % in^2
    A_ShearOut = sum((S.rowSpacing*(n-1) + S.firstRowZ) * S.t * 2 * S.nPinsPerRow); % Sum all shearout areas
    
    % Net tension area: cross-sectional area minus pin holes at ONE row (critical section)
    % Each row is at a different axial location, so we only subtract one row's pins
    A_wall_gross = pi*0.25*((S.ID + 2*S.t)^2 - S.ID^2);  % Gross wall cross-section
    A_pin_holes_per_row = S.pinDia * S.nPinsPerRow * S.t; % Area removed by pins in one row
    A_tension = A_wall_gross - A_pin_holes_per_row;
    
    % Safety check: net area must be positive
    if A_tension <= 0
        A_tension = 0.001;  % Prevent division by zero, will show very high stress
    end
    %%%%%%%%%%%%%%%% ----- Load calcs ----- %%%%%%%%%%%%%%%%
    p_design = S.MEOP_psi * S.DF;     % psi
    F_axial  = p_design * A_bore;     % lbf
    F_perPin = F_axial / totalPins;   % lbf average
    Stress_ShearOut = F_axial/A_ShearOut/1000; % KSI
    Net_Tension = (F_axial/A_tension)/1000;    % KSI
    Pin_Shear = (F_perPin/(pi*0.25*(S.pinDia)^2))/1000; % KSI
    Bearing = (F_perPin/(S.pinDia*S.t))/1000; % KSI
    Hoop = (p_design*(((S.ID+2*S.t)^2 + (S.ID)^2)/((S.ID+2*S.t)^2 - (S.ID)^2)))/1000; % KSI
    Axial = Hoop/2; % KSI
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
        p_design, F_axial, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);
    set(S.txtStress, 'String', stressOutText);
    
    % Update constraint status panel
    S = updateConstraintStatus(S, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);
    drawnow limitrate;
end
%% ---------- UPDATE CONSTRAINT STATUS ----------
function S = updateConstraintStatus(S, shearOut, netTension, pinShear, bearing, hoop, axial, FOS)
    % Check each constraint
    status = {};
    allPass = true;
    
    % Shear Out
    if shearOut <= S.targets.shearOut_max
        status{end+1} = sprintf('Shear Out: %.2f <= %.2f OK', shearOut, S.targets.shearOut_max);
    else
        status{end+1} = sprintf('Shear Out: %.2f > %.2f FAIL', shearOut, S.targets.shearOut_max);
        allPass = false;
    end
    
    % Net Tension
    if netTension <= S.targets.netTension_max
        status{end+1} = sprintf('Net Tens:  %.2f <= %.2f OK', netTension, S.targets.netTension_max);
    else
        status{end+1} = sprintf('Net Tens:  %.2f > %.2f FAIL', netTension, S.targets.netTension_max);
        allPass = false;
    end
    
    % Pin Shear
    if pinShear <= S.targets.pinShear_max
        status{end+1} = sprintf('Pin Shear: %.2f <= %.2f OK', pinShear, S.targets.pinShear_max);
    else
        status{end+1} = sprintf('Pin Shear: %.2f > %.2f FAIL', pinShear, S.targets.pinShear_max);
        allPass = false;
    end
    
    % Bearing
    if bearing <= S.targets.bearing_max
        status{end+1} = sprintf('Bearing:   %.2f <= %.2f OK', bearing, S.targets.bearing_max);
    else
        status{end+1} = sprintf('Bearing:   %.2f > %.2f FAIL', bearing, S.targets.bearing_max);
        allPass = false;
    end
    
    % Hoop
    if hoop <= S.targets.hoop_max
        status{end+1} = sprintf('Hoop:      %.2f <= %.2f OK', hoop, S.targets.hoop_max);
    else
        status{end+1} = sprintf('Hoop:      %.2f > %.2f FAIL', hoop, S.targets.hoop_max);
        allPass = false;
    end
    
    % Axial
    if axial <= S.targets.axial_max
        status{end+1} = sprintf('Axial:     %.2f <= %.2f OK', axial, S.targets.axial_max);
    else
        status{end+1} = sprintf('Axial:     %.2f > %.2f FAIL', axial, S.targets.axial_max);
        allPass = false;
    end
    
    % FOS
    if FOS >= S.targets.pinShearFOS_min
        status{end+1} = sprintf('Pin FOS:   %.2f >= %.2f OK', FOS, S.targets.pinShearFOS_min);
    else
        status{end+1} = sprintf('Pin FOS:   %.2f < %.2f FAIL', FOS, S.targets.pinShearFOS_min);
        allPass = false;
    end
    
    % Add summary
    if allPass
        status{end+1} = '--- ALL CONSTRAINTS MET ---';
        set(S.pOptResults, 'BackgroundColor', [0.9 1 0.9]);
    else
        status{end+1} = '--- CONSTRAINTS VIOLATED ---';
        set(S.pOptResults, 'BackgroundColor', [1 0.9 0.9]);
    end
    
    set(S.txtOptResults, 'String', strjoin(status, '\n'));
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