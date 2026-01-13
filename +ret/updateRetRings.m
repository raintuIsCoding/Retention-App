function S = updateRetRings(S, r_i)
%RET.UPDATERETRINGS Render retention rings at the casing ends.
%
% Requirements:
% - Ring OD matches casing ID (outer radius = r_i)
% - Render one end in detail view; render both ends in full view
% - Rings align to the casing ends and extend inward by S.retLen_end
% - Ring wall thickness uses S.retRingThk (same as mass model)

% Safety defaults
if ~isfield(S,'retRingThk') || ~isfinite(S.retRingThk) || S.retRingThk <= 0
    S.retRingThk = 0.25; % in
end
if ~isfield(S,'retLen_end') || ~isfinite(S.retLen_end) || S.retLen_end <= 0
    % Fall back to geometric logic used in updateAll
    S.retLen_end = max(1.0, 2*S.firstRowZ + (S.nRows - 1)*S.rowSpacing);
end
if ~isfield(S,'L_casing') || ~isfinite(S.L_casing) || S.L_casing <= 0
    return;
end

% Radii
r_o_ring = r_i;                         % OD matches casing ID
r_i_ring = max(0, r_o_ring - S.retRingThk);

theta = linspace(0, 2*pi, S.nCirc);

% Build parametric cylinder surfaces at x in [0, L]
% Outer and inner sleeves
Xspan = [0, S.retLen_end];
[Th, X] = meshgrid(theta, Xspan);

Y_outer = r_o_ring*cos(Th);
Z_outer = r_o_ring*sin(Th);

Y_inner = r_i_ring*cos(Th);
Z_inner = r_i_ring*sin(Th);

% Inner face (disc) at the inboard end of each ring
% (at x = S.retLen_end for aft ring; x = L_casing - S.retLen_end for fore ring)
R = linspace(r_i_ring, r_o_ring, 2);
[ThF, RF] = meshgrid(theta, R);
Y_face = RF.*cos(ThF);
Z_face = RF.*sin(ThF);

% End 1 (aft): ring spans x=[0, retLen]
set(S.hRetRingOuter(1), 'XData', X, 'YData', Y_outer, 'ZData', Z_outer);
set(S.hRetRingInner(1), 'XData', X, 'YData', Y_inner, 'ZData', Z_inner);

X_face1 = ones(size(Y_face))*S.retLen_end;
set(S.hRetRingFace(1), 'XData', X_face1, 'YData', Y_face, 'ZData', Z_face);

% Outboard cap at the casing end (x = 0)
X_cap1 = zeros(size(Y_face));
set(S.hRetRingCap(1), 'XData', X_cap1, 'YData', Y_face, 'ZData', Z_face);

% End 2 (fore): mirror along midplane by translating to the far end
X2 = X + (S.L_casing - S.retLen_end);
set(S.hRetRingOuter(2), 'XData', X2, 'YData', Y_outer, 'ZData', Z_outer);
set(S.hRetRingInner(2), 'XData', X2, 'YData', Y_inner, 'ZData', Z_inner);

X_face2 = ones(size(Y_face))*(S.L_casing - S.retLen_end);
set(S.hRetRingFace(2), 'XData', X_face2, 'YData', Y_face, 'ZData', Z_face);

% Outboard cap at the casing end (x = L_casing)
X_cap2 = ones(size(Y_face))*S.L_casing;
set(S.hRetRingCap(2), 'XData', X_cap2, 'YData', Y_face, 'ZData', Z_face);

% Visibility toggle
if ~isfield(S,'showRetRings') || ~islogical(S.showRetRings)
    S.showRetRings = true;
end

if ~S.showRetRings
    set([S.hRetRingOuter(1) S.hRetRingInner(1) S.hRetRingFace(1) S.hRetRingCap(1) ...
         S.hRetRingOuter(2) S.hRetRingInner(2) S.hRetRingFace(2) S.hRetRingCap(2)], 'Visible','off');
    return;
end

% Visibility by mode
if isfield(S,'mirrorPins') && S.mirrorPins
    set([S.hRetRingOuter(1) S.hRetRingInner(1) S.hRetRingFace(1) S.hRetRingCap(1) ...
         S.hRetRingOuter(2) S.hRetRingInner(2) S.hRetRingFace(2) S.hRetRingCap(2)], 'Visible','on');
else
    set([S.hRetRingOuter(1) S.hRetRingInner(1) S.hRetRingFace(1) S.hRetRingCap(1)], 'Visible','on');
    set([S.hRetRingOuter(2) S.hRetRingInner(2) S.hRetRingFace(2) S.hRetRingCap(2)], 'Visible','off');
end

end
