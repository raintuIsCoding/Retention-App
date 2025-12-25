function S = updateCasingSurfaces(S, r_i, r_o)

theta = linspace(0, 2*pi, S.nCirc);
x     = [0, S.L_casing];
[Theta, X] = meshgrid(theta, x);

Yo = r_o*cos(Theta);
Zo = r_o*sin(Theta);
Yi = r_i*cos(Theta);
Zi = r_i*sin(Theta);

set(S.hCasingOuter, 'XData', X, 'YData', Yo, 'ZData', Zo);
set(S.hCasingInner, 'XData', X, 'YData', Yi, 'ZData', Zi);

Xcap0 = [zeros(1,S.nCirc); zeros(1,S.nCirc)];
Ycap0 = [r_i*cos(theta);   r_o*cos(theta)];
Zcap0 = [r_i*sin(theta);   r_o*sin(theta)];
set(S.hCap0, 'XData', Xcap0, 'YData', Ycap0, 'ZData', Zcap0);

XcapL = [S.L_casing*ones(1,S.nCirc); S.L_casing*ones(1,S.nCirc)];
YcapL = [r_i*cos(theta);            r_o*cos(theta)];
ZcapL = [r_i*sin(theta);            r_o*sin(theta)];
set(S.hCapL, 'XData', XcapL, 'YData', YcapL, 'ZData', ZcapL);

end
