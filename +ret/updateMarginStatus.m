function S = updateMarginStatus(S, cases)
% cases: struct array with fields:
%   .name (char/string)
%   .load (numeric)
%   .allow (numeric)
%   .fos (numeric)

% Header is line 1; rows start at line 2
if ~isfield(S,'txtOptLines') || isempty(S.txtOptLines) || ~all(isgraphics(S.txtOptLines))
    return;
end

% Safety
nShow = min(numel(cases), 6); % keep UI compact; adjust if you add more
passAll = true;

for i = 1:nShow
    nm = string(cases(i).name);
    L  = cases(i).load;
    A  = cases(i).allow;
    F  = cases(i).fos;

    % Robust clamps
    if ~isfinite(L) || L <= 0, L = NaN; end
    if ~isfinite(A) || A <= 0, A = NaN; end
    if ~isfinite(F) || F <= 0, F = NaN; end

    margin = (A/(L*F)) - 1;
    isPass = isfinite(margin) && (margin >= 0);
    passAll = passAll && isPass;

    status = "PASS";
    if ~isPass, status = "FAIL"; end

    % Pretty fixed-width columns
    % Case(14) Load(8) Allow(8) FOS(5) Margin(8) Status
    line = sprintf('%-14s %8.2f %8.2f %5.2f %8.2f   %s', ...
        char(nm), L, A, F, margin, char(status));

    set(S.txtOptLines(i+1), 'String', line);

    % Color per row
    if isPass
        set(S.txtOptLines(i+1), 'ForegroundColor', [0 0.45 0]);
    else
        set(S.txtOptLines(i+1), 'ForegroundColor', [0.75 0 0]);
    end
end

% Clear unused lines (up to line 7)
for k = (nShow+2):7
    set(S.txtOptLines(k), 'String', '', 'ForegroundColor', [0 0 0]);
end

% Banner line 8
if passAll
    set(S.txtOptLines(8), 'String', '--- ALL MARGINS PASS ---', 'ForegroundColor', [0 0.45 0]);
else
    set(S.txtOptLines(8), 'String', '--- MARGINS VIOLATED ---', 'ForegroundColor', [0.75 0 0]);
end

end