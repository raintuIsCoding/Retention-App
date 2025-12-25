function S = updateAll(S)

% Derived geometry
r_i = S.ID/2;
r_o = r_i + S.t;

% Packing check
r_center = r_o;
circumference = 2*pi*r_center;
minPitch = S.minCircPitchFactor * S.pinDia;
maxPinsCirc = floor(circumference / minPitch);
packingOK = (S.nPinsPerRow <= maxPinsCirc);

% Axial check
lastRowX = S.firstRowZ + (S.nRows - 1)*S.rowSpacing;
axialOK = (lastRowX + S.edgeMarginEnd <= S.L_casing);

% Update casing + pins
S = ret.updateCasingSurfaces(S, r_i, r_o);
S = ret.updatePins(S, r_o);

%% ===== MASS CALCS =====
if ~isfield(S,'retRingThk') || ~isfinite(S.retRingThk) || S.retRingThk<=0
    S.retRingThk = 0.25;
end
retention_ring_thickness = S.retRingThk;

OD = S.ID + 2*S.t;
Volume_casing = (pi/4) * (OD^2 - S.ID^2) * S.casingLength;
Mass_casing = Volume_casing * S.density_CF;

totalPins_oneEnd = S.nRows * S.nPinsPerRow;
retention_length = S.edgeMarginEnd*2 + S.rowSpacing*(S.nRows-1);
Volume_retention_ring = (pi/4) * (S.ID^2 - (S.ID - 2*retention_ring_thickness)^2) * retention_length;
Volume_pins_in_ring = totalPins_oneEnd * pi * (S.pinDia/2)^2 * retention_ring_thickness;
Mass_retention_rings = (Volume_retention_ring - Volume_pins_in_ring) * S.density_Al * 2;

pinLen = S.t + retention_ring_thickness;  % derived engagement length
Volume_all_pins = totalPins_oneEnd * 2 * pi * (S.pinDia/2)^2 * pinLen;
Mass_pins = Volume_all_pins * S.density_Al;

Total_Mass = Mass_casing + Mass_retention_rings + Mass_pins;

geomOutText = sprintf([ ...
    'Total Pins (one end): %d\n' ...
    'Last Row X: %.2f in\n' ...
    'Max Pins/Row: %d\n' ...
    'Packing OK: %s\n' ...
    'Axial OK: %s\n' ...
    '-------------------\n' ...
    'Casing Mass: %.2f lb\n' ...
    'Retention Mass: %.2f lb\n' ...
    'Pin Mass: %.2f lb\n' ...
    'TOTAL MASS: %.2f lb\n' ], ...
    totalPins_oneEnd, lastRowX, maxPinsCirc, ret.tfStr(packingOK), ret.tfStr(axialOK), ...
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
kInteract = 3.0;
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