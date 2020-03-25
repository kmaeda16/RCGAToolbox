function [f, g] = wrapper_hiv_con_mex(x)

% x = 10 .^ x;

% x = [ 3.76e+2 6.57e+0 5.65e+5 2.24e+1 6.28e+0 2.46e+1 2.35e+1 2.72e+1 1.82e+1 1.73e+1 5.09e-3 5.17e-3 6.00e-3 4.50e-3 4.15e-3 -3.93e-3 -6.31e-3 -1.63e-2 -9.32e-3 1.02e-4 ];

[f, g] = hiv_con_mex(x);

% f = sum(x.^2);%1e+22;              % For debug
% g = 1e+22 * ones(1,12); % For debug
% if rand < 0.5
%     g = zeros(1,12);
% end

% if max(isnan([f g])) > 0
%     f = 1e+22;
%     g = 1e+22 * ones(1,12);
% end
