function S = updateConstraintStatus(S, pressureVessel, netTension, shearOut, bearing, pinShear)

% Expect 8 lines: 7 constraints + 1 banner
if ~isfield(S,'txtOptLines') || numel(S.txtOptLines) < 8 || ~all(isgraphics(S.txtOptLines))
    return; % UI not built yet
end

targets = S.targets;

% Soft PASS margin (5% by default)
marginFrac = S.marginFrac;   % 0.5 = allow 5% violation to still show PASS

% Formatting helpers
fmtLine = @(name, val, lim, pass, sense) sprintf('%-16s: %6.2f %2s %6.2f  %s', ...
    name, val, sense, lim, tern(pass,'PASS','FAIL'));

colPass = [0.85 1.00 0.85];
colFail = [1.00 0.85 0.85];

passAll = true;

% 1 Pressure Vessel (<=)
lim = targets.pressure_vessel_max;
softLim = lim * (1 + marginFrac);
pass = (pressureVessel <= softLim); passAll = passAll && pass;
setLine(1, fmtLine('Pressure Vessel', pressureVessel, lim, pass, '<='), pass);

% 2 Net Tension (<=)
lim = targets.netTension_max;
softLim = lim * (1 + marginFrac);
pass = (netTension <= softLim); passAll = passAll && pass;
setLine(2, fmtLine('Net Tension', netTension, lim, pass, '<='), pass);

% 3 Shear Out (<=)
lim = targets.shearOut_max;
softLim = lim * (1 + marginFrac);
pass = (shearOut <= softLim); passAll = passAll && pass;
setLine(3, fmtLine('Shear Out', shearOut, lim, pass, '<='), pass);

% 4 Bearing (<=)
lim = targets.bearing_max;
softLim = lim * (1 + marginFrac);
pass = (bearing <= softLim); passAll = passAll && pass;
setLine(4, fmtLine('Bearing', bearing, lim, pass, '<='), pass);

% 5 Pin Shear (<=)
lim = targets.pinShear_max;
softLim = lim * (1 + marginFrac);
pass = (pinShear <= softLim); passAll = passAll && pass;
setLine(5, fmtLine('Pin Shear', pinShear, lim, pass, '<='), pass);

% 7 Banner + panel background
if passAll
    set(S.pOptResults, 'BackgroundColor', [0.90 1.00 0.90]);
    setLine(6, '--- ALL CONSTRAINTS MET ---', true);
else
    set(S.pOptResults, 'BackgroundColor', [1.00 0.90 0.90]);
    setLine(6, '--- CONSTRAINTS VIOLATED ---', false);
end

    function setLine(i, str, isPass)
        set(S.txtOptLines(i), 'String', str);
        set(S.txtOptLines(i), 'BackgroundColor', tern(isPass, colPass, colFail));
        set(S.txtOptLines(i), 'ForegroundColor', [0 0 0]);
    end

end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end
