function S = updateConstraintStatus(S, shearOut, netTension, pinShear, bearing, hoop, axial, FOS)

status = {};
allPass = true;

if shearOut <= S.targets.shearOut_max
    status{end+1} = sprintf('Shear Out: %.2f <= %.2f OK', shearOut, S.targets.shearOut_max);
else
    status{end+1} = sprintf('Shear Out: %.2f > %.2f FAIL', shearOut, S.targets.shearOut_max);
    allPass = false;
end

if netTension <= S.targets.netTension_max
    status{end+1} = sprintf('Net Tens:  %.2f <= %.2f OK', netTension, S.targets.netTension_max);
else
    status{end+1} = sprintf('Net Tens:  %.2f > %.2f FAIL', netTension, S.targets.netTension_max);
    allPass = false;
end

if pinShear <= S.targets.pinShear_max
    status{end+1} = sprintf('Pin Shear: %.2f <= %.2f OK', pinShear, S.targets.pinShear_max);
else
    status{end+1} = sprintf('Pin Shear: %.2f > %.2f FAIL', pinShear, S.targets.pinShear_max);
    allPass = false;
end

if bearing <= S.targets.bearing_max
    status{end+1} = sprintf('Bearing:   %.2f <= %.2f OK', bearing, S.targets.bearing_max);
else
    status{end+1} = sprintf('Bearing:   %.2f > %.2f FAIL', bearing, S.targets.bearing_max);
    allPass = false;
end

if hoop <= S.targets.hoop_max
    status{end+1} = sprintf('Hoop:      %.2f <= %.2f OK', hoop, S.targets.hoop_max);
else
    status{end+1} = sprintf('Hoop:      %.2f > %.2f FAIL', hoop, S.targets.hoop_max);
    allPass = false;
end

if axial <= S.targets.axial_max
    status{end+1} = sprintf('Axial:     %.2f <= %.2f OK', axial, S.targets.axial_max);
else
    status{end+1} = sprintf('Axial:     %.2f > %.2f FAIL', axial, S.targets.axial_max);
    allPass = false;
end

if FOS >= S.targets.pinShearFOS_min
    status{end+1} = sprintf('Pin FOS:   %.2f >= %.2f OK', FOS, S.targets.pinShearFOS_min);
else
    status{end+1} = sprintf('Pin FOS:   %.2f < %.2f FAIL', FOS, S.targets.pinShearFOS_min);
    allPass = false;
end

if allPass
    status{end+1} = '--- ALL CONSTRAINTS MET ---';
    set(S.pOptResults, 'BackgroundColor', [0.9 1 0.9]);
else
    status{end+1} = '--- CONSTRAINTS VIOLATED ---';
    set(S.pOptResults, 'BackgroundColor', [1 0.9 0.9]);
end

set(S.txtOptResults, 'String', strjoin(status, '\n'));

end
