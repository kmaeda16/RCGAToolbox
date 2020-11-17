function [c, ceq] = cst_wrapper_nodecoding(fitnessfun, x)
% cst_wrapper_nodecoding returns the constraint function value vector g of 
% fitnessfun.
% 
% [SYNTAX]
% [c, ceq] = cst_wrapper_nodecoding(fitnessfun, x)
% 
% [INPUT]
% fitnessfun :  Function handle for a fitness function.
% x          :  Decision variables.
% 
% [OUTPUT]
% c          :  Inequality constraint function values (g).
% ceq        :  Equality constraint function values (empty).


[~, c] = fitnessfun(x);
ceq = [];
