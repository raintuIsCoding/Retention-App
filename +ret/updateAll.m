function S = updateAll(S)

% =========================
% Derived geometry
% =========================
r_i = S.ID/2;
r_o = r_i + S.t;

% -----------------------
% Pin counts (define once)
% -----------------------
totalPins_oneEnd = S.nRows * S.nPinsPerRow;   % pins engaged at ONE end (for axial load split)
totalPins_total  = totalPins_oneEnd * 2;      % full hardware mass assumes BOTH ends

% =========================
% Retention ring length logic
% =========================
% Row-driven required ring length per end (this SHOULD grow with nRows)
retLen_rows = 2*S.firstRowZ + (S.nRows - 1)*S.rowSpacing;

% Optional: user-controlled minimum retention ring length (per end)
if ~isfield(S,'retRingLength') || ~isfinite(S.retRingLength) || S.retRingLength < 0
    S.retRingLength = 1.0;     % default minimum per end
end
retLen_manual = S.retRingLength;

% Actual ring length used everywhere (never shorter than what rows require)
retLen_end = max(retLen_manual, retLen_rows);

% =========================
% Casing length logic (render + mass)
% =========================
if ~isfield(S,'casingBaseLength') || ~isfinite(S.casingBaseLength) || S.casingBaseLength <= 0
    S.casingBaseLength = 90;   % fixed minimum tube length
end
L_base = S.casingBaseLength;

S.casingLength_eff = L_base + 2*retLen_end;
S.L_casing = S.casingLength_eff;   % renderer uses this

% =========================
% Packing check
% =========================
r_center = r_o;
circumference = 2*pi*r_center;
minPitch = S.minCircPitchFactor * S.pinDia;
maxPinsCirc = floor(circumference / minPitch);
packingOK = (S.nPinsPerRow <= maxPinsCirc);

% =========================
% Axial geometry check (now uses correct length)
% =========================
lastRowX = S.firstRowZ + (S.nRows - 1)*S.rowSpacing;
axialOK = (lastRowX + S.firstRowZ <= S.L_casing);

% =========================
% Update casing + pins (rendering)
% =========================
S = ret.updateCasingSurfaces(S, r_i, r_o);
S = ret.updatePins(S, r_o);

%% ===== MASS CALCS =====
nEnds_mass = 2;

% --- Retention ring thickness (in) ---
if ~isfield(S,'retRingThk') || ~isfinite(S.retRingThk) || S.retRingThk <= 0
    S.retRingThk = 0.25;  % default guess; can be overridden
end
retRingThk = S.retRingThk;

% -----------------------
% 1) Casing mass (CF)
% -----------------------
OD = S.ID + 2*S.t;
A_annulus      = (pi/4) * (OD^2 - S.ID^2);
Volume_casing  = A_annulus * S.casingLength_eff;

% Subtract pin-hole volume in the casing wall (simple cylindrical hole model)
Volume_holes_casing = totalPins_total * pi * (S.pinDia/2)^2 * S.t;
Volume_casing_net   = Volume_casing - Volume_holes_casing;

Mass_casing = Volume_casing_net * S.density_CF;

% -----------------------
% 2) Retention ring mass (Al)
% -----------------------
A_ring = (pi/4) * (S.ID^2 - (S.ID - 2*retRingThk)^2);

Volume_ring_oneEnd      = A_ring * retLen_end;
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
    'Ret Ring Len (rows req, per end): %.2f in\n' ...
    'Ret Ring Len (manual min, per end): %.2f in\n' ...
    'Ret Ring Len (USED, per end): %.2f in\n' ...
    'Base Casing Length: %.2f in\n' ...
    'Total Casing Length: %.2f in\n' ...
    'Max Pins/Row: %d\n' ...
    'Packing OK: %s\n' ...
    'Axial OK: %s\n' ...
    '-------------------\n' ...
    'Casing Mass: %.2f lb\n' ...
    'Retention Mass: %.2f lb\n' ...
    'Pin Mass: %.2f lb\n' ...
    '\n' ...
    'TOTAL MASS: %.2f lb\n' ], ...
    totalPins_oneEnd, totalPins_total, lastRowX, ...
    retLen_rows, retLen_manual, retLen_end, ...
    L_base, S.casingLength_eff, maxPinsCirc, ...
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

F_perPin_casing = F_axial_casing / max(1, totalPins_oneEnd);
F_perPin_pin    = F_axial_pin    / max(1, totalPins_oneEnd);

% -------- Stresses (KSI) --------
Stress_ShearOut = (F_axial_casing / A_ShearOut) / 1000;
Net_Tension     = (F_axial_casing / A_tension)  / 1000;

A_pin = (pi/4) * S.pinDia^2;
if ~isfinite(A_pin) || A_pin <= 0, A_pin = 1e-6; end

Pin_Shear = (F_perPin_pin / A_pin) / 1000;
Bearing  = (F_perPin_casing / (S.pinDia * S.t)) / 1000;

% Thick-wall pressure vessel hoop at inner surface (KSI)
Pressure_Vessel  = (p_design_casing * (((S.ID+2*S.t)^2 + S.ID^2) / ((S.ID+2*S.t)^2 - (S.ID)^2))) / 1000;

if isfield(S,'allowedPinDias')
    [~,i] = min(abs(S.allowedPinDias - S.pinDia));
    S.pinDia = S.allowedPinDias(i);
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
    'Pressure Vessel: %.4f KSI\n' ], ...
    p_design_casing, p_design_pin, F_axial_casing, F_axial_pin, ...
    F_perPin_casing, F_perPin_pin, ...
    Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Pressure_Vessel);

set(S.txtStress, 'String', stressOutText);

S = ret.updateConstraintStatus(S, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Pressure_Vessel);

drawnow limitrate;

end