function S = readConstraintsFromUI(S)
S.targets.shearOut_max    = str2double(get(S.edTgtShearOut, 'String'));
S.targets.netTension_max  = str2double(get(S.edTgtNetTens, 'String'));
S.targets.pinShear_max    = str2double(get(S.edTgtPinShear, 'String'));
S.targets.bearing_max     = str2double(get(S.edTgtBearing, 'String'));
S.targets.hoop_max        = str2double(get(S.edTgtHoop, 'String'));
S.targets.axial_max       = str2double(get(S.edTgtAxial, 'String'));
S.targets.pinShearFOS_min = str2double(get(S.edTgtFOS, 'String'));
end
