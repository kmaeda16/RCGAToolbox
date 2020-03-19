function [f, g] = wrapper_threestep_con_mex(x)

x = 10 .^ x;

[f, g] = threestep_con_mex(x);

% if max(isnan([f g])) > 0
%     f = 1e+22;
%     g = 1e+22 * ones(1,24);
% end
