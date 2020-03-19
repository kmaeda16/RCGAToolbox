function [c, ceq] = cst_wrapper_nodecoding(fitnessfun, x)
% cst_wrapper_nodecoding returns g of fitnessfun.
% 
% [SYNTAX]
% [c, ceq] = cst_wrapper(fitnessfun, decodingfun, gene)
% 
% [INPUT]
% x   : Decision variables
% 
% [OUTPUT]
% c   :  Inequality constraint function values
% ced :  Equality constraint function values

% x = decodingfun(gene);
[~, c] = fitnessfun(x);
ceq = [];
