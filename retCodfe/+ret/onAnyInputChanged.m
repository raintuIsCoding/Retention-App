function onAnyInputChanged(src, ~)

% Pull state FIRST
S = guidata(src);

% Detect control type
ctrlStyle = get(src, 'Style');

switch ctrlStyle
    case 'popupmenu'
        % Discrete pin diameter selection
        if isfield(S,'allowedPinDias')
            idx = get(src,'Value');
            idx = max(1, min(idx, numel(S.allowedPinDias)));
            S.pinDia = S.allowedPinDias(idx);
        end

    otherwise
        % For edit boxes, toggles, etc.
        S = ret.readControlsToState(S);
end

% Update physics + UI
S = ret.updateAll(S);

% Push state back
guidata(S.fig, S);

end