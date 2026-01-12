function S = updateMarginTable(S, rows)
% rows is struct array with fields:
%  .component (string/char)
%  .caseName  (string/char)
%  .marginY   (numeric or NaN)   % Yield margin
%  .marginU   (numeric or NaN)   % Ultimate margin

if ~isfield(S,'marginCells') || isempty(S.marginCells) || ~all(isgraphics(S.marginCells(:)))
    return;
end

green = [0.85 1.00 0.85];
red   = [1.00 0.85 0.85];
white = [1 1 1];

n = min(numel(rows), size(S.marginCells,1));

for i = 1:n
    set(S.marginCells(i,1), 'String', char(string(rows(i).component)));
    set(S.marginCells(i,2), 'String', char(string(rows(i).caseName)));

    % Yield margin cell
    my = rows(i).marginY;
    if ~isfinite(my)
        set(S.marginCells(i,3), 'String', 'N/A', 'BackgroundColor', white, 'ForegroundColor', [0 0 0]);
    else
        set(S.marginCells(i,3), 'String', sprintf('%+.2f', my));
        if my >= 0
            set(S.marginCells(i,3), 'BackgroundColor', green, 'ForegroundColor', [0 0.35 0]);
        else
            set(S.marginCells(i,3), 'BackgroundColor', red, 'ForegroundColor', [0.65 0 0]);
        end
    end

    % Ultimate margin cell
    mu = rows(i).marginU;
    if ~isfinite(mu)
        set(S.marginCells(i,4), 'String', 'N/A', 'BackgroundColor', white, 'ForegroundColor', [0 0 0]);
    else
        set(S.marginCells(i,4), 'String', sprintf('%+.2f', mu));
        if mu >= 0
            set(S.marginCells(i,4), 'BackgroundColor', green, 'ForegroundColor', [0 0.35 0]);
        else
            set(S.marginCells(i,4), 'BackgroundColor', red, 'ForegroundColor', [0.65 0 0]);
        end
    end
end

% Clear remaining rows
for i = (n+1):size(S.marginCells,1)
    for c = 1:4
        set(S.marginCells(i,c), 'String', '', 'BackgroundColor', white, 'ForegroundColor', [0 0 0]);
    end
end

end