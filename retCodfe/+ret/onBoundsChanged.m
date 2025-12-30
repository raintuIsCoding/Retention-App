function onBoundsChanged(src, ~)
S = guidata(src);
S = ret.readBoundsFromUI(S);
S = ret.readConstraintsFromUI(S);
guidata(S.fig, S);
end
