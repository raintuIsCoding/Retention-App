function retention_app_V4
clc; close all;

S = ret.defaults();
S = ret.buildUI(S);        % creates figure + controls + stores handles in S
S = ret.updateAll(S);      % computes + renders once
guidata(S.fig, S);         % store state

end