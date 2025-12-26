function [S, success] = gridSearchOptimization(S)
success = false;

fprintf('\n=== RUNNING GRID SEARCH OPTIMIZATION ===\n');

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

p_design = MEOP_psi * DF;
A_bore = pi * (ID/2)^2;
F_axial = p_design * A_bore;
fprintf('Design axial force: %.0f lbf\n', F_axial);

if bounds.t(2) > bounds.t(1), t_vals = linspace(bounds.t(1), bounds.t(2), 8); else, t_vals = bounds.t(1); end
nRows_vals = bounds.nRows(1):bounds.nRows(2);
if bounds.nPinsPerRow(2) > bounds.nPinsPerRow(1), nPinsPerRow_vals = bounds.nPinsPerRow(1):2:bounds.nPinsPerRow(2); else, nPinsPerRow_vals = bounds.nPinsPerRow(1); end
if bounds.rowSpacing(2) > bounds.rowSpacing(1), rowSpacing_vals = linspace(bounds.rowSpacing(1), bounds.rowSpacing(2), 6); else, rowSpacing_vals = bounds.rowSpacing(1); end
if bounds.firstRowZ(2) > bounds.firstRowZ(1), firstRowZ_vals = linspace(bounds.firstRowZ(1), bounds.firstRowZ(2), 5); else, firstRowZ_vals = bounds.firstRowZ(1); end
% Option A: pin decision variable is an index into allowedPinDias
pinIdx_vals = 1:numel(S.allowedPinDias);

bestObj = inf;
bestParams = [];

totalCombinations = numel(t_vals)*numel(nRows_vals)*numel(nPinsPerRow_vals)*numel(rowSpacing_vals)*numel(firstRowZ_vals)*numel(pinIdx_vals);
fprintf('Total combinations to evaluate: %d\n', totalCombinations);

nEval = 0;
nFeasible = 0;

for t = t_vals
    for nRows = nRows_vals
        for nPinsPerRow = nPinsPerRow_vals
            for rowSpacing = rowSpacing_vals
                for firstRowZ = firstRowZ_vals
                    for pinIdx = pinIdx_vals
                        x = [t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinIdx];

                        nEval = nEval + 1;

                        [c, ~] = ret.nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, ...
                             minCircPitchFactor, minAxialPitchFactor, ...
                             S.retRingThk, targets, S.allowedPinDias);

                        if all(c <= 0)
                            nFeasible = nFeasible + 1;
                            obj = ret.objectiveFunction(x, ID, casingLength, density_CF, density_Al, S.density_pin, S.retRingThk, S.allowedPinDias);

                            if obj < bestObj
                                bestObj = obj;
                                bestParams = x;
                            end
                        end
                    end
                end
            end
        end
    end
end

if ~isempty(bestParams)
    S.t = bestParams(1);
    S.nRows = round(bestParams(2));
    S.nPinsPerRow = round(bestParams(3));
    S.rowSpacing = bestParams(4);
    S.firstRowZ = bestParams(5);
    S.pinDiaIdx = round(bestParams(6));
    S.pinDiaIdx = max(1, min(S.pinDiaIdx, numel(S.allowedPinDias)));
    S.pinDia = S.allowedPinDias(S.pinDiaIdx);
    S.pinLen = S.t + S.retRingThk;
    S.edgeMarginEnd = S.t;

    success = true;
    fprintf('\n=== GRID SEARCH RESULTS ===\n');
    fprintf('Wall Thickness: %.4f in\n', S.t);
    fprintf('Axial Rows: %d\n', S.nRows);
    fprintf('Pins per Row: %d\n', S.nPinsPerRow);
    fprintf('Row Spacing: %.4f in\n', S.rowSpacing);
    fprintf('First Row Offset: %.4f in\n', S.firstRowZ);
    fprintf('Pin Diameter: %.4f in\n', S.pinDia);
    fprintf('Objective Value: %.4f\n', bestObj);
    fprintf('Feasible solutions found: %d out of %d\n', nFeasible, nEval);
else
    fprintf('No feasible solution found in grid search.\n');
end

end
