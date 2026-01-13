function S = updateAll(S)

% =========================
% Derived geometry
% =========================
r_i = S.ID/2;
r_o = r_i + S.t;

% =========================
% Snap pin diameter early (if you use discrete sizes)
% =========================
if isfield(S,'allowedPinDias') && ~isempty(S.allowedPinDias) && isfinite(S.pinDia)
    [~,i] = min(abs(S.allowedPinDias - S.pinDia));
    S.pinDia = S.allowedPinDias(i);
end

% -----------------------
% Pin counts (define once)
% -----------------------
totalPins_oneEnd = S.nRows * S.nPinsPerRow;   % pins engaged at ONE end (for axial load split)
totalPins_total  = totalPins_oneEnd * 2;      % full hardware mass assumes BOTH ends

% =========================
% Retention ring length logic (per end)
% =========================
% Driven by: edge offset (firstRowZ), row spacing, number of rows
retLen_rows = 2*S.firstRowZ + (S.nRows - 1)*S.rowSpacing;

% Optional minimum (keeps it from going absurdly small)
if ~isfield(S,'retRingLengthMin') || ~isfinite(S.retRingLengthMin) || S.retRingLengthMin < 0
    S.retRingLengthMin = 1.0; % in, per end
end

S.retLen_end = max(S.retRingLengthMin, retLen_rows);

% =========================
% Physical casing length logic (true length)
% =========================
if ~isfield(S,'casingBaseLength') || ~isfinite(S.casingBaseLength) || S.casingBaseLength <= 0
    S.casingBaseLength = 90.0; % fixed minimum tube length
end

S.L_phys = S.casingBaseLength + 2*S.retLen_end;
S.casingLength_eff = S.L_phys; % mass uses physical length

% =========================
% Packing check
% =========================
r_center = r_o;
circumference = 2*pi*r_center;
minPitch = S.minCircPitchFactor * S.pinDia;
maxPinsCirc = floor(circumference / minPitch);
packingOK = (S.nPinsPerRow <= maxPinsCirc);

% =========================
% Axial geometry check (must use PHYSICAL length)
% =========================
lastRowX = S.firstRowZ + (S.nRows - 1)*S.rowSpacing;
axialOK = (lastRowX + S.firstRowZ <= S.L_phys);

% =========================
% Apply render mode (full vs config) + force visual refresh
% (This sets S.L_casing for drawing only)
% =========================
S = ret.updateView(S, r_o);

% =========================
% Update casing + pins (rendering)
% =========================
S = ret.updateCasingSurfaces(S, r_i, r_o);
S = ret.updatePins(S, r_o);
S = ret.updateRetRings(S, r_i);

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
if ~isfinite(Volume_casing_net) || Volume_casing_net < 0
    Volume_casing_net = max(0, Volume_casing); % safety clamp
end

Mass_casing = Volume_casing_net * S.density_CF;

% -----------------------
% 2) Retention ring mass (Al)
% -----------------------
A_ring = (pi/4) * (S.ID^2 - (S.ID - 2*retRingThk)^2);

Volume_ring_oneEnd      = A_ring * S.retLen_end;
Volume_pinHoles_oneEnd  = totalPins_oneEnd * pi * (S.pinDia/2)^2 * retRingThk;

Volume_ring_net_oneEnd  = Volume_ring_oneEnd - Volume_pinHoles_oneEnd;
if ~isfinite(Volume_ring_net_oneEnd) || Volume_ring_net_oneEnd < 0
    Volume_ring_net_oneEnd = max(0, Volume_ring_oneEnd);
end

Mass_retention_rings = Volume_ring_net_oneEnd * S.density_Al * nEnds_mass;

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
    'Ret Ring Length: %.2f in\n' ...
    'Base Casing Length: %.2f in\n' ...
    'Full Casing Length: %.2f in\n' ...
    '-------------------\n' ...
    'Casing Mass: %.2f lb\n' ...
    'Retention Mass: %.2f lb\n' ...
    'Pin Mass: %.2f lb\n' ...
    '\n' ...
    'TOTAL MASS: %.2f lb\n' ], ...
    totalPins_oneEnd, totalPins_total, ...
    S.retLen_end, ...
    S.casingBaseLength, S.L_phys, ...
    Mass_casing, Mass_retention_rings, Mass_pins, Total_Mass);

if isfield(S,'txtGeom') && isgraphics(S.txtGeom)
    set(S.txtGeom, 'String', geomOutText);
end

%% ===== STRESS CALCS =====
n = 1:S.nRows;

% -------- Areas --------
A_bore = pi * (S.ID/2)^2;  % in^2

% Shear-out area (your existing style; kept)
A_ShearOut = sum((S.rowSpacing*(n-1) + S.firstRowZ) * S.t * 2 * S.nPinsPerRow);

% Net tension area
A_wall_gross    = (pi/4) * ((S.ID + 2*S.t)^2 - S.ID^2);
A_holes_per_row = S.pinDia * S.nPinsPerRow * S.t;

% Hole interaction model (make kInteract a tunable state var)
if ~isfield(S,'kInteract') || ~isfinite(S.kInteract) || S.kInteract < 0
    S.kInteract = 1.0; % default: interaction if spacing <= 1*D
end

if S.rowSpacing <= S.kInteract * S.pinDia
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
Bearing   = (F_perPin_casing / (S.pinDia * S.t)) / 1000;

% Thick-wall pressure vessel hoop at inner surface (KSI)
OD = S.ID + 2*S.t;
Pressure_Vessel  = (p_design_casing * ((OD^2 + S.ID^2) / (OD^2 - S.ID^2))) / 1000;

stressOutText = sprintf([ ...
    'Design Pressure (Casing): %.1f psi\n' ...
    'Design Pressure (Pins):   %.1f psi\n' ...
    'Axial End Force (Casing): %.0f lbf\n' ...
    'Axial End Force (Pins):   %.0f lbf\n' ...
    'Avg Force per Pin (Casing DF): %.1f lbf\n' ...
    'Avg Force per Pin (Pin DF):    %.1f lbf\n' ...
    '-------------------\n' ...
    '\n' ...
    'Shear Out Stress: %.4f KSI\n' ...
    'Net Tension: %.4f KSI\n' ...
    'Pin Shear: %.4f KSI\n' ...
    'Bearing: %.4f KSI\n' ...
    'Pressure Vessel: %.4f KSI\n' ], ...
    p_design_casing, p_design_pin, F_axial_casing, F_axial_pin, ...
    F_perPin_casing, F_perPin_pin, ...
    Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Pressure_Vessel);

if isfield(S,'txtStress') && isgraphics(S.txtStress)
    set(S.txtStress, 'String', stressOutText);
end

% ===== MARGIN STATUS (MEOP load, DF as FOS) =====
p_meop = S.MEOP_psi;             % psi
F_axial_meop = p_meop * A_bore;  % lbf

F_perPin_meop = F_axial_meop / max(1, totalPins_oneEnd);

% Loads at MEOP (KSI)
Load_ShearOut_MEOP = (F_axial_meop / A_ShearOut) / 1000;
Load_NetTens_MEOP  = (F_axial_meop / A_tension)  / 1000;
Load_Bearing_MEOP  = (F_perPin_meop / (S.pinDia * S.t)) / 1000;

A_pin = (pi/4)*S.pinDia^2;
Load_PinShear_MEOP = (F_perPin_meop / max(1e-6, A_pin)) / 1000;

OD = S.ID + 2*S.t;
Load_PV_MEOP = (p_meop * ((OD^2 + S.ID^2) / (OD^2 - S.ID^2))) / 1000;

cases = struct([]);

cases(1).name  = 'Shear Out';
cases(1).load  = Load_ShearOut_MEOP;
cases(1).allow = S.allow.yield.shearOut;
cases(1).fos   = S.DF_casing;

cases(2).name  = 'Net Tension';
cases(2).load  = Load_NetTens_MEOP;
cases(2).allow = S.allow.yield.netTension;
cases(2).fos   = S.DF_casing;

cases(3).name  = 'Pin Shear';
cases(3).load  = Load_PinShear_MEOP;
cases(3).allow = S.allow.yield.pinShear;
cases(3).fos   = S.DF_pin;

cases(4).name  = 'Bearing';
cases(4).load  = Load_Bearing_MEOP;
cases(4).allow = S.allow.yield.bearing;
cases(4).fos   = S.DF_casing;

cases(5).name  = 'P Vessel';
cases(5).load  = Load_PV_MEOP;
cases(5).allow = S.allow.yield.pressureVessel;
cases(5).fos   = S.DF_casing;

S = ret.updateMarginStatus(S, cases);

S = ret.updateConstraintStatus(S, Stress_ShearOut, Net_Tension, Pin_Shear, Bearing, Pressure_Vessel);

% ===== Margin Table (Yield + Ultimate) =====
% Load values should be at MEOP (NOT multiplied by DF/FOS)
p_meop = S.MEOP_psi;
F_axial_meop = p_meop * A_bore;
F_perPin_meop = F_axial_meop / max(1, totalPins_oneEnd);

Load_ShearOut = (F_axial_meop / A_ShearOut) / 1000;
Load_NetTens  = (F_axial_meop / A_tension)  / 1000;
Load_Bearing  = (F_perPin_meop / (S.pinDia * S.t)) / 1000;

A_pin = (pi/4)*S.pinDia^2;
Load_PinShear = (F_perPin_meop / max(1e-6, A_pin)) / 1000;

OD = S.ID + 2*S.t;
Load_PV = (p_meop * ((OD^2 + S.ID^2) / (OD^2 - S.ID^2))) / 1000;

marg = @(allow, load, fos) (allow./(load.*fos)) - 1;

rows = struct([]);

% ---- CASING ----
rows(1).component = "Casing";
rows(1).caseName  = "Shear-Out";
rows(1).marginY   = marg(S.allow.yield.shearOut, Load_ShearOut, S.FOS.casing);
rows(1).marginU   = marg(S.allow.ult.shearOut,   Load_ShearOut, S.FOS.casing);

rows(2).component = "Casing";
rows(2).caseName  = "Net Tension";
rows(2).marginY   = marg(S.allow.yield.netTension, Load_NetTens, S.FOS.casing);
rows(2).marginU   = marg(S.allow.ult.netTension,   Load_NetTens, S.FOS.casing);

rows(3).component = "Casing";
rows(3).caseName  = "Bearing";
rows(3).marginY   = marg(S.allow.yield.bearing, Load_Bearing, S.FOS.casing);
rows(3).marginU   = marg(S.allow.ult.bearing,   Load_Bearing, S.FOS.casing);

rows(4).component = "Casing";
rows(4).caseName  = "Pressure Vessel";
rows(4).marginY   = marg(S.allow.yield.pressureVessel, Load_PV, S.FOS.casing);
rows(4).marginU   = marg(S.allow.ult.pressureVessel,   Load_PV, S.FOS.casing);

% ---- PINS ----
rows(5).component = "Pins";
rows(5).caseName  = "Pin Shear";
rows(5).marginY   = marg(S.allow.yield.pinShear, Load_PinShear, S.FOS.comp_y);
rows(5).marginU   = marg(S.allow.ult.pinShear,   Load_PinShear, S.FOS.comp_u);

% ---- RET RING ---- (defaults to N/A if allowables are NaN)
rows(6).component = "Ret Ring";
rows(6).caseName  = "Shear-Out";
rows(6).marginY   = marg(S.allow.yield.ret_shearOut, Load_ShearOut, S.FOS.comp_y);
rows(6).marginU   = marg(S.allow.ult.ret_shearOut,   Load_ShearOut, S.FOS.comp_u);

rows(7).component = "Ret Ring";
rows(7).caseName  = "Net Tension";
rows(7).marginY   = marg(S.allow.yield.ret_netTension, Load_NetTens, S.FOS.comp_y);
rows(7).marginU   = marg(S.allow.ult.ret_netTension,   Load_NetTens, S.FOS.comp_u);

rows(8).component = "Ret Ring";
rows(8).caseName  = "Bearing";
rows(8).marginY   = marg(S.allow.yield.ret_bearing, Load_Bearing, S.FOS.comp_y);
rows(8).marginU   = marg(S.allow.ult.ret_bearing,   Load_Bearing, S.FOS.comp_u);

S = ret.updateMarginTable(S, rows);

drawnow limitrate;

end