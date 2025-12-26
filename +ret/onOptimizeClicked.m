function onOptimizeClicked(src, ~)
S = guidata(src);

S = ret.readControlsToState(S);

% ensure pinDia matches dropdown if present
if isfield(S,'ddPinDia') && isgraphics(S.ddPinDia)
    idx = get(S.ddPinDia,'Value');
    idx = max(1, min(idx, numel(S.allowedPinDias)));
    S.pinDia = S.allowedPinDias(idx);
end

S = ret.readConstraintsFromUI(S);
S = ret.readBoundsFromUI(S);

set(S.txtOptStatus, 'String', 'Optimizing... Please wait.');
drawnow;

[S, success] = ret.runOptimization(S);

if success
    set(S.edT, 'String', num2str(S.t, '%.4f'));
    set(S.edRows, 'String', num2str(S.nRows));
    set(S.edPPR, 'String', num2str(S.nPinsPerRow));
    set(S.edRowSp, 'String', num2str(S.rowSpacing, '%.4f'));
    set(S.edFirst, 'String', num2str(S.firstRowZ, '%.4f'));
    set(S.edPinDia, 'String', num2str(S.pinDia, '%.4f'));
    set(S.txtOptStatus, 'String', 'Optimization Complete!');
else
    set(S.txtOptStatus, 'String', 'Optimization Failed - Try adjusting bounds');
end

S = ret.updateAll(S);
guidata(S.fig, S);

end
