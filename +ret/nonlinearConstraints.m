function [c, ceq] = nonlinearConstraints(x, ID, MEOP_psi, DF, L_casing, minCircPitchFactor, minAxialPitchFactor, retRingThk, targets, allowedPinDias)

% ret.nonlinearConstraints
% Nonlinear inequality constraints c(x) <= 0
%
% x = [t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinDia]

ceq = [];

% Unpack
t           = x(1);
nRows       = round(x(2));
nPinsPerRow = round(x(3));
rowSpacing  = x(4);
firstRowZ   = x(5);
pinIdx = round(x(6));
pinIdx = max(1, min(pinIdx, numel(allowedPinDias)));
pinDia = allowedPinDias(pinIdx);

% Guard against nonsense (GA sometimes samples negatives)
nRows       = max(1, nRows);
nPinsPerRow = max(1, nPinsPerRow);
t           = max(1e-6, t);
rowSpacing  = max(1e-6, rowSpacing);
firstRowZ   = max(0, firstRowZ);
pinDia      = max(1e-6, pinDia);

% Derived radii
r_i = ID/2;
r_o = r_i + t;

% -----------------------
% Geometry constraints
% -----------------------

% (1) Circumferential packing
circumference = 2*pi*r_o;
minPitch = minCircPitchFactor * pinDia;
maxPinsCirc = floor(circumference / minPitch);
c_pack = nPinsPerRow - maxPinsCirc;   % <= 0

% (2) Axial pitch minimum
c_axPitch = (minAxialPitchFactor*pinDia) - rowSpacing;  % <= 0  (i.e., rowSpacing >= minAxialPitchFactor*pinDia)

% (3) Axial extent must fit in modeled region (symmetric margins based on firstRowZ)
lastRowX = firstRowZ + (nRows - 1)*rowSpacing;
% Need enough room for "same margin" after last row too:
c_axExtent = (lastRowX + firstRowZ) - L_casing;         % <= 0

% -----------------------
% Stress constraints
% -----------------------
% Loads
A_bore   = pi * r_i^2;
p_design = MEOP_psi * DF;
F_axial  = p_design * A_bore;  % lbf

pins_oneEnd = nRows * nPinsPerRow;
F_perPin = F_axial / max(1, pins_oneEnd);

% Areas and stresses (match updateAll logic)
n = 1:nRows;

A_ShearOut = sum((rowSpacing*(n-1) + firstRowZ) * t * 2 * nPinsPerRow);

A_wall_gross    = (pi/4) * ((ID + 2*t)^2 - ID^2);
A_holes_per_row = pinDia * nPinsPerRow * t;

kInteract = 3.0;
if rowSpacing <= kInteract * pinDia
    N_eff = nRows;
else
    N_eff = 1;
end
A_tension = A_wall_gross - N_eff * A_holes_per_row;

% clamps
if ~isfinite(A_tension)  || A_tension  <= 0, A_tension  = 1e-6; end
if ~isfinite(A_ShearOut) || A_ShearOut <= 0, A_ShearOut = 1e-6; end

% KSI
Stress_ShearOut = (F_axial / A_ShearOut) / 1000;
Net_Tension     = (F_axial / A_tension)  / 1000;

A_pin = (pi/4) * pinDia^2;
if ~isfinite(A_pin) || A_pin <= 0, A_pin = 1e-6; end
Pin_Shear = (F_perPin / A_pin) / 1000;

Bearing = (F_perPin / (pinDia * t)) / 1000;

Hoop = (p_design * (((ID+2*t)^2 + ID^2) / ((ID+2*t)^2 - (ID)^2))) / 1000;
Axial = Hoop/2;

% Pin shear FOS (keep your same model/offset)
Pin_Shear_Strength = 43.5; % KSI
if Pin_Shear <= 0
    Pin_Shear_FOS = inf;
else
    Pin_Shear_FOS = (Pin_Shear_Strength/Pin_Shear) - 2;
end

% Inequalities: stress <= allowable  -> (stress - allowable) <= 0
c_shearOut = Stress_ShearOut - targets.shearOut_max;
c_netTen   = Net_Tension     - targets.netTension_max;
c_pinShear = Pin_Shear       - targets.pinShear_max;
c_bearing  = Bearing         - targets.bearing_max;
c_hoop     = Hoop            - targets.hoop_max;
c_axial    = Axial           - targets.axial_max;

% FOS constraint: FOS >= min -> (min - FOS) <= 0
c_fos = targets.pinShearFOS_min - Pin_Shear_FOS;

% Bundle constraints
c = [
    c_pack
    c_axPitch
    c_axExtent
    c_shearOut
    c_netTen
    c_pinShear
    c_bearing
    c_hoop
    c_axial
    c_fos
];

end