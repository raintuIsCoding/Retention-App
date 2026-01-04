function S = defaults()
%RET.DEFAULTS Initialize default state struct

%% ---------- DEFAULT INPUTS ----------
S.MIII = (8/6);
S.marginFrac = 0.01;
S.ID          = 6.0;      % inner diameter of casing (in)
S.t           = 0.20;     % wall thickness (in)
S.L_casing    = 10.0;     % length of modeled region (in)
S.lengthMode  = "aftOnly"; % "aftOnly" (10 in) or "full" (96 in)
S.mirrorPins  = false;     % when true, mirror pins to forward end for display
S.lengthMode  = "aftOnly"; % "aftOnly" (10 in) or "full" (96 in)
S.mirrorPins  = false;     % when true, mirror pins to forward end for display
S.nRows       = 3;        % axial rows of pins
S.nPinsPerRow = 10;       % pins around circumference per row
S.rowSpacing  = 0.5;      % axial spacing between rows (in)
S.firstRowZ   = 0.75;     % axial position of first row from x=0 (in)
S.pinDia      = 0.25;    % pin diameter (in)
S.retRingThk  = 0.25;     % retention ring thickness (in)
S.pinLen      = S.t + S.retRingThk;  % derived pin engagement length (in)

% Loads inputs
S.MEOP_psi   = 650;   % psi
S.DF_casing = 2;
S.DF_pin    = 1.75;

S.DF_casing_from_CLT = 1; % 3.4;

S.theta_scale     = deg2rad(50);
S.theta_ref = atan(sqrt(2));    % 54.7356 deg

S.hoop_scaling_from_theta  = 1; % (sin(S.theta_scale).^2) ./ (sin(S.theta_ref).^2);
S.axial_scaling_from_theta = 1; %(cos(S.theta_scale).^2) ./ (cos(S.theta_ref).^2);

% Pattern toggle
S.pinPatternMode = "progressive";   % "progressive" or "alternating"
S.altStartPhase  = 0;               % 0 or 1

S.allowedPinDias = [0.25, 0.3125, 0.375, 0.5];  % in

%% ---------- TARGET STRESS LIMITS ---------- (All set by MIII max pressure motor pressure)
S.targets.shearOut_max    = 2.45;    % KSI
S.targets.netTension_max  = 10.82;    % KSI
S.targets.pinShear_max    = 21.84;          % KSI
S.targets.bearing_max     = 24.5;    % KSI
S.targets.hoop_max        = 20.17; % * (S.DF_casing_from_CLT/S.DF_casing) * S.hoop_scaling_from_theta;        % KSI
S.targets.axial_max       = 10.09; % * (S.DF_casing_from_CLT/S.DF_casing) * S.axial_scaling_from_theta;        % KSI
S.targets.pinShearFOS_min = 1.43;    % Minimum FOS

%% ---------- CASING PARAMETERS ----------
S.casingLength = 91;              % Total casing length (in) - for weight calc
S.density_CF = 0.0535;            % Carbon fiber density (lb/in^3)
S.density_Al = 0.0975;            % 6061 Aluminum density (lb/in^3)
S.density_pin = 0.289;            % lb/in^3 for 316 SS (typical)

%% ---------- GEOMETRIC CONSTRAINTS ----------
S.minCircPitchFactor  = 3.0;       % min center spacing = factor * pinDia
S.minAxialPitchFactor = 2.5;       % min axial spacing = factor * pinDia
S.edgeMarginEnd       = S.t;       % safety margin from last row to casing end

%% ---------- OPTIMIZATION BOUNDS ----------
S.optBounds.t           = [0.2, 0.75];
S.optBounds.nRows       = [1, 5];
S.optBounds.nPinsPerRow = [6, 30];
S.optBounds.rowSpacing  = [0.5, 2.0];
S.optBounds.firstRowZ   = [0.75, 2.5];
S.optBounds.pinDia      = [0.25, 0.5];

%% ---------- RENDER CONSTANTS ----------
S.nCirc    = 150;  % casing mesh resolution
S.nPinCirc = 30;   % pin mesh resolution
S.maxPins  = 400;  % pool size

end
