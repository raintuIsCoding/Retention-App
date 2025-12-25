function obj = objectiveFunction(x, ID, MEOP_psi, DF, casingLength, density_CF, density_Al)
t = x(1);
nRows = round(x(2));
nPinsPerRow = round(x(3));
rowSpacing = x(4);
firstRowZ = x(5);
pinDia = x(6);

OD = ID + 2*t;
Volume_casing = (pi/4) * (OD^2 - ID^2) * casingLength;
Mass_casing = Volume_casing * density_CF;

retention_ring_thickness = 0.25;
edgeMarginEnd = t;
retention_length = edgeMarginEnd*2 + rowSpacing*(nRows-1);
Volume_retention_ring = (pi/4) * (ID^2 - (ID - 2*retention_ring_thickness)^2) * retention_length;

Volume_pins_in_ring = nPinsPerRow * nRows * pi * (pinDia/2)^2 * retention_ring_thickness;
Volume_retention_net = Volume_retention_ring - Volume_pins_in_ring;
Mass_retention_rings = Volume_retention_net * density_Al * 2;

retRingThk = 0.25;
pinLen = t + retRingThk;
totalPins = nRows * nPinsPerRow * 2;
Volume_all_pins = totalPins * pi * (pinDia/2)^2 * pinLen;
Mass_pins = Volume_all_pins * density_Al;

totalMass = Mass_casing + Mass_retention_rings + Mass_pins;

lastRowX = firstRowZ + (nRows - 1) * rowSpacing;

w_mass = 1.0;
w_extension = 0.1;

obj = w_mass * (totalMass / 20.0) + w_extension * (lastRowX / 5.0);

end
