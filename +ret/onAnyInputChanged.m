function onAnyInputChanged(src, ~)
S = guidata(src);
S = ret.readControlsToState(S);
S = ret.updateAll(S);
guidata(S.fig, S);
end
