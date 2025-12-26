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
%    Assumption: ring has through-holes for pins;
%    we subtract the cylindrical hole volume using pinDia and retRingThk.
% -----------------------
retLen = 2*S.firstRowZ + (S.nRows-1)*S.rowSpacing;  % axial ring length
A_ring = (pi/4) * (S.ID^2 - (S.ID - 2*retRingThk)^2);

Volume_ring_oneEnd      = A_ring * retLen;
Volume_pinHoles_oneEnd  = totalPins_oneEnd * pi * (S.pinDia/2)^2 * retRingThk;

Mass_retention_rings = (Volume_ring_oneEnd - Volume_pinHoles_oneEnd) ...
                        * S.density_Al * nEnds_mass;

% -----------------------
% 3) Pin mass (steel)
%    Engagement length derived from casing thickness + ring thickness.
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
% Notes:
% - Hoop stress uses Lame thick-wall hoop stress at the INNER surface (conservative); assumes p_external = 0.
% - Net-tension hole interaction: assume all rows contribute when rows are "close" (conservative).
%   A future improvement would be to reduce the effective rows counted when row spacing is large.

n = 1:S.nRows;

% -------- Areas --------
A_bore = pi * (S.ID/2)^2;  % in^2

% Shear-out area (your current model)
A_ShearOut = sum((S.rowSpacing*(n-1) + S.firstRowZ) * S.t * 2 * S.nPinsPerRow);

% Net tension area
A_wall_gross    = (pi/4) * ((S.ID + 2*S.t)^2 - S.ID^2);   % gross wall cross-section
A_holes_per_row = S.pinDia * S.nPinsPerRow * S.t;         % "slot" approx: d*t per pin per row

% Conservative hole interaction model:
% If rows are within ~kInteract diameters, treat rows as one interacting damaged zone.
% If rows are spaced farther than that, this "all rows interact" assumption becomes redundant.
kInteract = 1.0;
if S.rowSpacing <= kInteract * S.pinDia
    N_eff = S.nRows;   % close proximity -> count all rows (conservative)
else
    N_eff = 1;         % spaced out -> single-row critical section (less conservative)
end

A_tension = A_wall_gross - N_eff * A_holes_per_row;

% Safety clamps to avoid divide-by-zero / NaNs
if ~isfinite(A_tension)  || A_tension  <= 0, A_tension  = 1e-6; end
if ~isfinite(A_ShearOut) || A_ShearOut <= 0, A_ShearOut = 1e-6; end

% -------- Loads --------
p_design = S.MEOP_psi * S.DF;     % psi
F_axial  = p_design * A_bore;     % lbf
F_perPin = F_axial / max(1, totalPins_oneEnd);   % lbf average per pin (one end)

% -------- Stresses (KSI) --------
Stress_ShearOut = (F_axial / A_ShearOut) / 1000;
% KSI - Axial end load divided by your modeled net shear-out area
% (2 * "distance to edge" * wall thickness summed across rows).
% This is a useful trend metric; it can become optimistic/nonphysical if pin size grows large
% relative to the remaining ligament geometry.

Net_Tension = (F_axial / A_tension) / 1000;
% KSI - Conservative net-section model: counts multiple rows as interacting when row spacing is small.
% If rows are far apart (spacing >> pin diameter), only a single-row critical section should dominate.
% That "row decoupling" case is not fully modeled yet; we approximate it via N_eff.

A_pin = (pi/4) * S.pinDia^2;    % in^2
if ~isfinite(A_pin) || A_pin <= 0, A_pin = 1e-6; end
Pin_Shear = (F_perPin / A_pin) / 1000;           % KSI - single-plane shear per pin (average)

Bearing = (F_perPin / (S.pinDia * S.t)) / 1000;  % KSI - bearing per pin (average)

% Lame thick-wall hoop stress at inner surface (conservative); p_external = 0
Hoop = (p_design * (((S.ID+2*S.t)^2 + S.ID^2) / ((S.ID+2*S.t)^2 - (S.ID)^2))) / 1000;
% KSI - Thick-wall (Lame) hoop stress at the inner surface.
% This will generally be higher than a thin-wall estimate used in some simpler calculators.

Axial = Hoop/2; % KSI - closed-end cylinder axial stress approximation

if isfield(S,'allowedPinDias')
    [~,i] = min(abs(S.allowedPinDias - S.pinDia));
    S.pinDia = S.allowedPinDias(i);
end

% -------- FOS --------
Pin_Shear_Strength = 43.5; % KSI
if ~isfinite(Pin_Shear) || Pin_Shear <= 0
    Pin_Shear_FOS = inf;
else
    Pin_Shear_FOS = (Pin_Shear_Strength/Pin_Shear) - 2; % keep your calibration offset
end

stressOutText = sprintf([ ...
    'Design Pressure: %.1f psi\n' ...
    'Axial End Force: %.0f lbf\n' ...
    'Avg Force per Pin: %.1f lbf\n' ...
    'Shear Out Stress: %.4f KSI\n' ...
    'Net Tension: %.4f KSI\n' ...
    'Pin Shear: %.4f KSI\n' ...
    'Bearing: %.4f KSI\n' ...
    'Hoop: %.4f KSI\n' ...
    'Axial: %.4f KSI\n' ...
    '\n' ...
    'Pin Shear FOS: %.4f \n' ], ...
    p_design, F_axial, F_perPin, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);

set(S.txtStress, 'String', stressOutText);

S = ret.updateConstraintStatus(S, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Hoop, Axial, Pin_Shear_FOS);

drawnow limitrate;

end