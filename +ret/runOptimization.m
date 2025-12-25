function [S, success] = runOptimization(S)
success = false;

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

x0 = [S.t, S.nRows, S.nPinsPerRow, S.rowSpacing, S.firstRowZ, S.pinDia]; %#ok<NASGU>

lb = [bounds.t(1), bounds.nRows(1), bounds.nPinsPerRow(1), ...
      bounds.rowSpacing(1), bounds.firstRowZ(1), bounds.pinDia(1)];
ub = [bounds.t(2), bounds.nRows(2), bounds.nPinsPerRow(2), ...
      bounds.rowSpacing(2), bounds.firstRowZ(2), bounds.pinDia(2)];

intcon = [2, 3];

objFun = @(x) ret.objectiveFunction(x, ID, MEOP_psi, DF, casingLength, density_CF, density_Al);
nonlconFun = @(x) ret.nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, ...
                                          minCircPitchFactor, minAxialPitchFactor, targets);

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
        S.t = xOpt(1);
        S.nRows = round(xOpt(2));
        S.nPinsPerRow = round(xOpt(3));
        S.rowSpacing = xOpt(4);
        S.firstRowZ = xOpt(5);
        S.pinDia = xOpt(6);
        S.pinLen = S.t + S.retRingThk;
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
        [S, success] = ret.gridSearchOptimization(S);
    end
catch ME
    fprintf('Optimization error: %s\n', ME.message);
    fprintf('Attempting grid search fallback...\n');
    [S, success] = ret.gridSearchOptimization(S);
end

end
