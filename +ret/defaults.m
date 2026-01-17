function S = defaults()
%RET.DEFAULTS Initialize default state struct

%% ---------- DEFAULT INPUTS ----------
S.MIII = (8/6);
S.marginFrac = 0.01;

% Main geometry defaults (updated)
S.ID          = 8.0;      % inner diameter of casing (in)
S.t           = 0.28;     % wall thickness (in)
S.L_casing    = 10.0;     % length of modeled region (in)
S.lengthMode  = "aftOnly"; % "aftOnly" (config view) or "full" (full casing)
S.mirrorPins  = false;     % when true, mirror pins/ret rings to forward end for display

S.nRows       = 3;        % axial rows of pins
S.nPinsPerRow = 10;       % pins around circumference per row
S.rowSpacing  = 0.75;     % axial spacing between rows (in)
S.firstRowZ   = 1.0;      % axial position of first row from x=0 (in)
S.pinDia      = 0.375;    % pin diameter (in)

S.retRingThk  = 0.25;     % retention ring thickness (in)
S.pinLen      = S.t + S.retRingThk;  % derived pin engagement length (in)

S.casingBaseLength   = 18.0; % fixed minimum physical casing length
S.retRingLengthMin   = 1.0;  % minimum retention length per end
S.configViewLength   = 10.0; % rendered length in config view

% Loads inputs
S.MEOP_psi   = 850;   % psi
S.DF_casing = 2;
S.DF_pin    = 1.75;

S.theta_scale     = deg2rad(50);
S.theta_ref = atan(sqrt(2));    % 54.7356 deg

% Pattern toggle
S.pinPatternMode = "spiral";   % "spiral" or "alternating"
S.altStartPhase  = 0;           % 0 or 1

S.allowedPinDias = [0.25, 0.3125, 0.375, 0.5];  % in

% Ret ring visibility (default on)
S.showRetRings = true;

%% ---------- DEFAULT ALLOWABLES (KSI) ----------
% Casing (single input; applied to yield + ultimate internally)
S.allow = struct();
S.allow.yield = struct();
S.allow.ult   = struct();

S.allow.yield.shearOut       = 36.27;
S.allow.yield.netTension     = 13.93;
S.allow.yield.bearing        = 31.42;
S.allow.yield.pressureVessel = 30.10;

% Mirror casing allowables to ultimate by default (single-box UI behavior)
S.allow.ult.shearOut         = S.allow.yield.shearOut;
S.allow.ult.netTension       = S.allow.yield.netTension;
S.allow.ult.bearing          = S.allow.yield.bearing;
S.allow.ult.pressureVessel   = S.allow.yield.pressureVessel;

% Pins
S.allow.yield.pinShear       = 58.0;
S.allow.ult.pinShear         = 104.0;

% Ret ring
S.allow.yield.ret_shearOut   = 20.0;
S.allow.ult.ret_shearOut     = 27.0;
S.allow.yield.ret_netTension = 35.0;
S.allow.ult.ret_netTension   = 42.0;
S.allow.yield.ret_bearing    = 58.0;
S.allow.ult.ret_bearing      = 88.0;

%% ---------- CASING PARAMETERS ----------
S.casingLength = 90;              % Total casing length (in) - for weight calc
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

if ~isfield(S,'minCircPitchFactor') || ~isfinite(S.minCircPitchFactor) || S.minCircPitchFactor <= 0
    S.minCircPitchFactor = 3.0;
end

end
