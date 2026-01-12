function S = readConstraintsFromUI(S)
S.targets.shearOut_max    = str2double(get(S.edTgtShearOut, 'String'));
S.targets.netTension_max  = str2double(get(S.edTgtNetTens, 'String'));
S.targets.pinShear_max    = str2double(get(S.edTgtPinShear, 'String'));
S.targets.bearing_max     = str2double(get(S.edTgtBearing, 'String'));
S.targets.pressure_vessel_max        = str2double(get(S.edTgtPressureVessel, 'String'));
end
