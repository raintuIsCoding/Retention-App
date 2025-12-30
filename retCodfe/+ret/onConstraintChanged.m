function onConstraintChanged(src, ~)
S = guidata(src);
S = ret.readConstraintsFromUI(S);
S = ret.updateAll(S);
guidata(S.fig, S);
end
