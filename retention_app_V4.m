function retention_app_V4
clc; close all;

S = ret.defaults();
S = ret.buildUI(S);        % creates figure + controls + stores handles in S

guidata(S.fig, S);         % store state ASAP (so callbacks always find S)

S = ret.updateAll(S);      % computes + renders once

guidata(S.fig, S);         % store state again after updateAll has populated outputs

end
