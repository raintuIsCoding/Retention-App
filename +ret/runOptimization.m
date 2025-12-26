function [S, success] = runOptimization(S)
success = false;

% Pull constants from S (optimizer treats these as fixed)
ID = S.ID;
MEOP_psi = S.MEOP_psi;
DF = S.DF;

minCircPitchFactor = S.minCircPitchFactor;
minAxialPitchFactor = S.minAxialPitchFactor;
L_casing = S.L_casing;

casingLength = S.casingLength;
density_CF   = S.density_CF;
density_Al   = S.density_Al;
density_pin  = S.density_pin;

% Retention ring thickness is fixed during optimization (not a decision var)
retRingThk = S.retRingThk;

targets = S.targets;
bounds  = S.optBounds;

% --- Pin diameter is a discrete choice: optimize its INDEX ---
if ~isfield(S,'allowedPinDias') || isempty(S.allowedPinDias)
    error('S.allowedPinDias must be defined (e.g., [0.25 0.3125 0.375 0.5]).');
end
[~, pinIdx0] = min(abs(S.allowedPinDias - S.pinDia));  % robust

x0 = [S.t, S.nRows, S.nPinsPerRow, S.rowSpacing, S.firstRowZ, pinIdx0]; %#ok<NASGU>

lb = [bounds.t(1), bounds.nRows(1), bounds.nPinsPerRow(1), ...
      bounds.rowSpacing(1), bounds.firstRowZ(1), 1];

ub = [bounds.t(2), bounds.nRows(2), bounds.nPinsPerRow(2), ...
      bounds.rowSpacing(2), bounds.firstRowZ(2), numel(S.allowedPinDias)];

% Integer decision variables: nRows, nPinsPerRow, pinIdx
intcon = [2, 3, 6];

% Updated objective + constraints (match updateAll physics)
allowedPinDias = S.allowedPinDias;

objFun = @(x) ret.objectiveFunction(x, ID, casingLength, density_CF, density_Al, density_pin, retRingThk, allowedPinDias);

nonlconFun = @(x) ret.nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, ...
                                          minCircPitchFactor, minAxialPitchFactor, ...
                                          retRingThk, targets, allowedPinDias);

options = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 200, ...
    'FunctionTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6, ...
    'UseParallel', false, ...
    'PlotFcn', []);

try
    [xOpt, fval, exitflag] = ga(objFun, 6, [], [], [], [], lb, ub, nonlconFun, intcon, options);

    if exitflag > 0
        % Write optimized design back into S
        S.t           = xOpt(1);
        S.nRows       = round(xOpt(2));
        S.nPinsPerRow = round(xOpt(3));
        S.rowSpacing  = xOpt(4);
        S.firstRowZ   = xOpt(5);

        pinIdx = round(xOpt(6));
        pinIdx = max(1, min(pinIdx, numel(S.allowedPinDias)));
        S.pinDia = S.allowedPinDias(pinIdx);
        S.pinDiaIdx = pinIdx;

        % Derived quantities (keep consistent with updateAll)
        S.pinLen = S.t + S.retRingThk;

        % Optional: remove legacy field if it exists (avoid confusion)
        if isfield(S,'edgeMarginEnd')
            S = rmfield(S,'edgeMarginEnd');
        end

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
        [S, success] = ret.gridSearchOptimization(S);
    end

catch ME
    fprintf('Optimization error: %s\n', ME.message);
    fprintf('Attempting grid search fallback...\n');
    [S, success] = ret.gridSearchOptimization(S);
end

end