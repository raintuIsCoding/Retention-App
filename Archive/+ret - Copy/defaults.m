function S = defaults()
%RET.DEFAULTS Initialize default state struct

%% ---------- DEFAULT INPUTS ----------
S.ID          = 8.0;      % inner diameter of casing (in)
S.t           = 0.25;     % wall thickness (in)
S.L_casing    = 10.0;     % length of modeled region (in)
S.nRows       = 3;        % axial rows of pins
S.nPinsPerRow = 12;       % pins around circumference per row
S.rowSpacing  = 0.5;      % axial spacing between rows (in)
S.firstRowZ   = 0.75;     % axial position of first row from x=0 (in)
S.pinDia      = 0.375;    % pin diameter (in)
S.pinLen      = 2 * S.t;  % how far pins stick out radially (in)

% Loads inputs
S.MEOP_psi   = 850;   % psi
S.DF         = 1.5;   % design factor (MEOP multiplier)

% Pattern toggle
S.pinPatternMode = "progressive";   % "progressive" or "alternating"
S.altStartPhase  = 0;               % 0 or 1

%% ---------- TARGET STRESS LIMITS ----------
S.targets.shearOut_max    = 1.837831702;    % KSI
S.targets.netTension_max  = 11.50766592;    % KSI
S.targets.pinShear_max    = 18.72;          % KSI
S.targets.bearing_max     = 18.37831702;    % KSI
S.targets.hoop_max        = 15.1125;        % KSI
S.targets.axial_max       = 7.55625;        % KSI
S.targets.pinShearFOS_min = 0.323717949;    % Minimum FOS (must be >= this)

%% ---------- CASING PARAMETERS ----------
S.casingLength = 96;              % Total casing length (in) - for weight calc
S.density_CF = 0.0535;            % Carbon fiber density (lb/in^3)
S.density_Al = 0.0975;            % 6061 Aluminum density (lb/in^3)

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
