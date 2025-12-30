function obj = objectiveFunction(x, ID, casingLength, density_CF, density_Al, density_pin, retRingThk, allowedPinDias)
% ret.objectiveFunction
% Objective: minimize total mass of full hardware (both ends), consistent with updateAll.
%
% x = [t, nRows, nPinsPerRow, rowSpacing, firstRowZ, pinDia]
% ID          = casing inner diameter (in)
% casingLength= full casing length for mass (in) (ex: 96)
% density_*   = lb/in^3
% retRingThk  = retention ring thickness (in)

% Unpack design variables
t           = x(1);
nRows       = round(x(2));
nPinsPerRow = round(x(3));
rowSpacing  = x(4);
firstRowZ   = x(5);
pinIdx = round(x(6));
pinIdx = max(1, min(pinIdx, numel(allowedPinDias)));
pinDia = allowedPinDias(pinIdx);

% Full hardware assumption
nEnds_mass = 2;

% Pin counts
pins_oneEnd  = max(0, nRows) * max(0, nPinsPerRow);
pins_total   = pins_oneEnd * nEnds_mass;

% -----------------------
% 1) Casing mass (CF)
% -----------------------
OD = ID + 2*t;
A_annulus = (pi/4) * (OD^2 - ID^2);
V_casing  = A_annulus * casingLength;

% Subtract casing hole volume (simple cylinder through wall thickness)
V_holes_casing = pins_total * pi * (pinDia/2)^2 * t;
V_casing_net   = V_casing - V_holes_casing;

Mass_casing = V_casing_net * density_CF;

% -----------------------
% 2) Retention ring mass (Al) (two rings)
% -----------------------
% Ring OD = casing ID (per your design choice)
ringOD = ID;
ringID = ringOD - 2*retRingThk;
A_ring = (pi/4) * (ringOD^2 - ringID^2);

% Ring axial length based on symmetric margins = firstRowZ on both ends
retLen = 2*firstRowZ + (max(0,nRows)-1)*rowSpacing;

V_ring_oneEnd = A_ring * retLen;

% Subtract pin hole volume through ring thickness (hole depth = retRingThk)
V_ring_holes_oneEnd = pins_oneEnd * pi * (pinDia/2)^2 * retRingThk;

Mass_rings = (V_ring_oneEnd - V_ring_holes_oneEnd) * density_Al * nEnds_mass;

% -----------------------
% 3) Pins mass (steel)
% -----------------------
pinLen = t + retRingThk;
V_pins = pins_total * pi * (pinDia/2)^2 * pinLen;
Mass_pins = V_pins * density_pin;

% Total mass objective
Mass_total = Mass_casing + Mass_rings + Mass_pins;

% Optional tiny preference for compact pinfield (helps GA choose "cleaner" designs)
lastRowX = firstRowZ + (max(0,nRows)-1)*rowSpacing;

w_mass      = 1.0;
w_extension = 0.05;  % small nudge only

obj = w_mass*Mass_total + w_extension*lastRowX;
end
