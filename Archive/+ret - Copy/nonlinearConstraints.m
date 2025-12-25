function [c, ceq] = nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, minCircPitchFactor, minAxialPitchFactor, targets)
t = x(1);
nRows = round(x(2));
nPinsPerRow = round(x(3));
rowSpacing = x(4);
firstRowZ = x(5);
pinDia = x(6);

r_i = ID / 2;
r_o = r_i + t;
edgeMarginEnd = t;

totalPins = nRows * nPinsPerRow;

% Circumferential packing
r_center = r_o;
circumference = 2 * pi * r_center;
minCircPitch = minCircPitchFactor * pinDia;
maxPinsCirc = floor(circumference / minCircPitch);
c_packing = nPinsPerRow - maxPinsCirc;

% Axial pitch
minAxialPitch = minAxialPitchFactor * pinDia;
c_axialPitch = minAxialPitch - rowSpacing;

% Axial extent
lastRowX = firstRowZ + (nRows - 1) * rowSpacing;
c_axial = (lastRowX + edgeMarginEnd) - L_casing;

% Areas
n = 1:nRows;
A_bore = pi * (ID/2)^2;
A_ShearOut = sum((rowSpacing*(n-1) + firstRowZ) * t * 2 * nPinsPerRow);

A_wall_gross = pi * 0.25 * ((ID + 2*t)^2 - ID^2);
A_pin_holes_per_row = pinDia * nPinsPerRow * t;
A_tension = A_wall_gross - A_pin_holes_per_row;

c_netArea = (0.3 * A_wall_gross) - A_tension;

if A_tension <= 0, A_tension = 0.001; end
if A_ShearOut <= 0, A_ShearOut = 0.001; end

p_design = MEOP_psi * DF;
F_axial = p_design * A_bore;
F_perPin = F_axial / max(1,totalPins);

Stress_ShearOut = F_axial / A_ShearOut / 1000;
Net_Tension     = (F_axial / A_tension) / 1000;
Pin_Shear       = (F_perPin / (pi * 0.25 * pinDia^2)) / 1000;
Bearing         = (F_perPin / (pinDia * t)) / 1000;
Hoop            = (p_design * (((ID + 2*t)^2 + ID^2) / ((ID + 2*t)^2 - ID^2))) / 1000;
Axial_stress    = Hoop / 2;

Pin_Shear_Strength = 43.5;
Pin_Shear_FOS = (Pin_Shear_Strength / Pin_Shear) - 2;

c_shearOut    = Stress_ShearOut - targets.shearOut_max;
c_netTension  = Net_Tension - targets.netTension_max;
c_pinShear    = Pin_Shear - targets.pinShear_max;
c_bearing     = Bearing - targets.bearing_max;
c_hoop        = Hoop - targets.hoop_max;
c_axialStress = Axial_stress - targets.axial_max;
c_FOS         = targets.pinShearFOS_min - Pin_Shear_FOS;

c = [c_packing; c_axialPitch; c_axial; c_netArea; ...
     c_shearOut; c_netTension; c_pinShear; c_bearing; c_hoop; c_axialStress; c_FOS];

ceq = [];

end
