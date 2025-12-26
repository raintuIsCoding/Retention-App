function S = updateAll(S)

% Derived geometry
r_i = S.ID/2;
r_o = r_i + S.t;

% -----------------------
% Pin counts (define once)
% -----------------------
totalPins_oneEnd = S.nRows * S.nPinsPerRow;   % pins engaged at ONE end (for axial load split)
totalPins_total  = totalPins_oneEnd * 2;      % full hardware mass assumes BOTH ends

% Packing check
r_center = r_o;
circumference = 2*pi*r_center;
minPitch = S.minCircPitchFactor * S.pinDia;
maxPinsCirc = floor(circumference / minPitch);
packingOK = (S.nPinsPerRow <= maxPinsCirc);

% Axial check
lastRowX = S.firstRowZ + (S.nRows - 1)*S.rowSpacing;
axialOK = (lastRowX + S.firstRowZ <= S.L_casing);

% Retention ring axial length (symmetric margins based on firstRowZ)
retLen = 2*S.firstRowZ + (S.nRows - 1)*S.rowSpacing;

% Update casing + pins
S = ret.updateCasingSurfaces(S, r_i, r_o);
S = ret.updatePins(S, r_o);

%% ===== MASS CALCS =====
% Mass model assumes FULL hardware (both ends) regardless of rendered length.
nEnds_mass = 2;

% --- Retention ring thickness (in) ---
if ~isfield(S,'retRingThk') || ~isfinite(S.retRingThk) || S.retRingThk <= 0
    S.retRingThk = 0.25;  % default guess; can be overridden
end
retRingThk = S.retRingThk;

% -----------------------
% 1) Casing mass (CF)
%    Option: subtract pin-hole volume through wall thickness
% -----------------------
OD = S.ID + 2*S.t;
A_annulus      = (pi/4) * (OD^2 - S.ID^2);
Volume_casing  = A_annulus * S.casingLength;

% Subtract pin-hole volume in the casing wall (simple cylindrical hole model)
Volume_holes_casing = totalPins_total * pi * (S.pinDia/2)^2 * S.t;
Volume_casing_net   = Volume_casing - Volume_holes_casing;

Mass_casing = Volume_casing_net * S.density_CF;

% -----------------------
% 2) Retention ring mass (Al)
% -----------------------
retLen = 2*S.firstRowZ + (S.nRows-1)*S.rowSpacing;  % axial ring length
A_ring = (pi/4) * (S.ID^2 - (S.ID - 2*retRingThk)^2);

Volume_ring_oneEnd      = A_ring * retLen;
Volume_pinHoles_oneEnd  = totalPins_oneEnd * pi * (S.pinDia/2)^2 * retRingThk;

Mass_retention_rings = (Volume_ring_oneEnd - Volume_pinHoles_oneEnd) ...
                        * S.density_Al * nEnds_mass;

% -----------------------
% 3) Pin mass (steel)
% -----------------------
pinLen = S.t + retRingThk;
Volume_pins_total = totalPins_total * pi * (S.pinDia/2)^2 * pinLen;
Mass_pins = Volume_pins_total * S.density_pin;

% -----------------------
% Total
% -----------------------
Total_Mass = Mass_casing + Mass_retention_rings + Mass_pins;

% -----------------------
% UI readout
% -----------------------
geomOutText = sprintf([ ...
    'Total Pins (one end): %d\n' ...
    'Total Pins (both ends): %d\n' ...
    'Last Row X: %.2f in\n' ...
    'Retention Ring Length: %.2f in\n' ...
    'Max Pins/Row: %d\n' ...
    'Packing OK: %s\n' ...
    'Axial OK: %s\n' ...
    '-------------------\n' ...
    'Casing Mass: %.2f lb\n' ...
    'Retention Mass: %.2f lb\n' ...
    'Pin Mass: %.2f lb\n' ...
    '\n' ...
    'TOTAL MASS: %.2f lb\n' ], ...
    totalPins_oneEnd, totalPins_total, lastRowX, retLen, maxPinsCirc, ...
    ret.tfStr(packingOK), ret.tfStr(axialOK), ...
    Mass_casing, Mass_retention_rings, Mass_pins, Total_Mass);

set(S.txtGeom, 'String', geomOutText);

%% ===== STRESS CALCS =====
n = 1:S.nRows;

% -------- Areas --------
A_bore = pi * (S.ID/2)^2;  % in^2

% Shear-out area
A_ShearOut = sum((S.rowSpacing*(n-1) + S.firstRowZ) * S.t * 2 * S.nPinsPerRow);

% Net tension area
A_wall_gross    = (pi/4) * ((S.ID + 2*S.t)^2 - S.ID^2);
A_holes_per_row = S.pinDia * S.nPinsPerRow * S.t;

% Hole interaction model
kInteract = 0;
if S.rowSpacing <= kInteract * S.pinDia
    N_eff = S.nRows;
else
    N_eff = 1;
end
A_tension = A_wall_gross - N_eff * A_holes_per_row;

% Safety clamps
if ~isfinite(A_tension)  || A_tension  <= 0, A_tension  = 1e-6; end
if ~isfinite(A_ShearOut) || A_ShearOut <= 0, A_ShearOut = 1e-6; end

% -------- Loads --------
p_design_casing = S.MEOP_psi * S.DF_casing;  % psi
p_design_pin    = S.MEOP_psi * S.DF_pin;     % psi

F_axial_casing = p_design_casing * A_bore;   % lbf
F_axial_pin    = p_design_pin    * A_bore;   % lbf

% Two per-pin forces:
F_perPin_casing = F_axial_casing / max(1, totalPins_oneEnd); % for member-bearing
F_perPin_pin    = F_axial_pin    / max(1, totalPins_oneEnd); % for pin shear

% -------- Stresses (KSI) --------
Stress_ShearOut = (F_axial_casing / A_ShearOut) / 1000;
Net_Tension     = (F_axial_casing / A_tension)  / 1000;

A_pin = (pi/4) * S.pinDia^2;
if ~isfinite(A_pin) || A_pin <= 0, A_pin = 1e-6; end

Pin_Shear = (F_perPin_pin / A_pin) / 1000;

% Bearing should be casing DF (member failure):
Bearing  = (F_perPin_casing / (S.pinDia * S.t)) / 1000;

% Lame thick-wall hoop stress at inner surface (conservative); p_external = 0
Hoop  = (p_design_casing * (((S.ID+2*S.t)^2 + S.ID^2) / ((S.ID+2*S.t)^2 - (S.ID)^2))) / 1000;
Axial = Hoop/2;

if isfield(S,'allowedPinDias')
    [~,i] = min(abs(S.allowedPinDias - S.pinDia));
    S.pinDia = S.allowedPinDias(i);
end

% -------- FOS --------
Pin_Shear_Strength = 75; % KSI
if ~isfinite(Pin_Shear) || Pin_Shear <= 0
    Pin_Shear_FOS = inf;
else
    Pin_Shear_FOS = (Pin_Shear_Strength/Pin_Shear) - 2;
end

stressOutText = sprintf([ ...
    'Design Pressure (Casing): %.1f psi\n' ...
    'Design Pressure (Pins):   %.1f psi\n' ...
    'Axial End Force (Casing): %.0f lbf\n' ...
    'Axial End Force (Pins):   %.0f lbf\n' ...
    'Avg Force per Pin (Casing DF): %.1f lbf\n' ...
    'Avg Force per Pin (Pin DF):    %.1f lbf\n' ...
    'Shear Out Stress: %.4f KSI\n' ...
    'Net Tension: %.4f KSI\n' ...
    'Pin Shear: %.4f KSI\n' ...
    'Bearing: %.4f KSI\n' ...
    'Hoop: %.4f KSI\n' ...
    'Axial: %.4f KSI\n' ...
    '\n' ...
    'Pin Shear FOS: %.4f \n' ], ...
    p_design_casing, p_design_pin, F_axial_casing, F_axial_pin, ...
    F_perPin_casing, F_perPin_pin, ...
    Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);

set(S.txtStress, 'String', stressOutText);

S = ret.updateConstraintStatus(S, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);

drawnow limitrate;

end