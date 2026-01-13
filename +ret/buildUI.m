function S = buildUI(S)
%RET.BUILDUI Create figure/axes + graphics pool + HUD/UI

%% ---------- FIGURE / AXES ----------
S.fig = figure('Renderer','opengl', 'Position', [100, 100, 1400, 800]);

S.ax = axes( ...
    'Parent', S.fig, ...
    'Units', 'normalized', ...
    'Position', [0.55 0.10 0.40 0.80]);
hold(S.ax, 'on');
axis(S.ax, 'equal');
axis(S.ax, 'vis3d');
axis(S.ax, 'off');
set(S.ax, 'Visible','off');
lighting(S.ax, 'gouraud');
camlight(S.ax, 'headlight');
view(S.ax, 180, 0);
rotate3d(S.fig, 'on');
pan(S.fig, 'off');
zoom(S.fig, 'off');

%% ---------- DRAW CASING (handles once) ----------
theta = linspace(0, 2*pi, S.nCirc);
x2    = [0, 1]; % placeholder; will be set in updateCasingSurfaces()
[Theta, X] = meshgrid(theta, x2); %#ok<ASGLU>

S.hCasingOuter = surf(S.ax, X, X, X, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','none');
S.hCasingInner = surf(S.ax, X, X, X, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none');
S.hCap0        = surf(S.ax, X, X, X, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none');
S.hCapL        = surf(S.ax, X, X, X, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none');

%% ---------- CREATE PIN GRAPHICS POOL ----------
S.pinSide = gobjects(S.maxPins,1);
S.pinCap  = gobjects(S.maxPins,1);
for k = 1:S.maxPins
    S.pinSide(k) = surf(S.ax, nan(2,S.nPinCirc), nan(2,S.nPinCirc), nan(2,S.nPinCirc), ...
        'FaceColor',[0.8 0.2 0.2], 'EdgeColor','none', 'Visible','off');
    S.pinCap(k) = patch(S.ax, nan(1,S.nPinCirc), nan(1,S.nPinCirc), nan(1,S.nPinCirc), ...
        [0.8 0.2 0.2], 'EdgeColor','none', 'Visible','off');
end

%% ---------- RETENTION RING GRAPHICS (2 ends) ----------
% Rendered as a simple thin-walled aluminum sleeve that slides into the casing ID.
% Geometry is set in ret.updateRetRings() (OD = casing ID; length = S.retLen_end).
x2 = [0, 1];
[Theta2, X2] = meshgrid(theta, x2); %#ok<ASGLU>

S.hRetRingOuter = gobjects(2,1);
S.hRetRingInner = gobjects(2,1);
S.hRetRingFace  = gobjects(2,1); % inner face (at the inboard end)
for e = 1:2
    S.hRetRingOuter(e) = surf(S.ax, X2, X2, X2, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor','none', 'Visible','off');
    S.hRetRingInner(e) = surf(S.ax, X2, X2, X2, 'FaceColor', [0.78 0.78 0.78], 'EdgeColor','none', 'Visible','off');
    S.hRetRingFace(e)  = surf(S.ax, X2, X2, X2, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor','none', 'Visible','off');
    S.hRetRingCap(e)   = surf(S.ax, X2, X2, X2, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor','none', 'Visible','off');
end

%% ---------- HUD + CONTROLS ----------
S = ret.buildHUDandControls(S);

% Store the (figure) handle for guidata robustness
guidata(S.fig, S);

end
