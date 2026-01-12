function S = updateView(S, r_o)
%RET.UPDATEVIEW Sets render length + UI label + axes limits (no physics here)

% Defaults (safe if older states exist)
if ~isfield(S,'configViewLength') || ~isfinite(S.configViewLength) || S.configViewLength <= 0
    S.configViewLength = 10.0; % inches
end
if ~isfield(S,'lengthMode') || ~(S.lengthMode=="full" || S.lengthMode=="aftOnly")
    S.lengthMode = "aftOnly";
end

% Render length (what the 3D casing surface uses)
if S.lengthMode == "full"
    S.L_casing = S.L_phys;
    S.mirrorPins = true;
else
    % Show enough to contain the aft retention region + a small buffer
    buf = 0.5; % inches
    L_need = S.retLen_end + buf;
    S.L_casing = max(S.configViewLength, L_need);
    S.mirrorPins = false;
end

% Update the toggle button text (show both view + physical so it's obvious)
if isfield(S,'btnLengthMode') && isgraphics(S.btnLengthMode)
    if S.lengthMode == "full"
        set(S.btnLengthMode,'Value',1);
        set(S.btnLengthMode,'String','Full Casing');
    else
        set(S.btnLengthMode,'Value',0);
        set(S.btnLengthMode,'String','Retention Detail View');
    end
end

% Axes limits (this is the "visual not updating" part most people miss)
if isfield(S,'ax') && isgraphics(S.ax)
    padR = max(0.25, 0.25*max(1,r_o));
    pinExtra = 0;
    if isfield(S,'pinLen') && isfinite(S.pinLen), pinExtra = S.pinLen; end

    xlim(S.ax, [0, S.L_casing]);
    limR = (r_o + pinExtra + padR);
    ylim(S.ax, [-limR, limR]);
    zlim(S.ax, [-limR, limR]);
end

end