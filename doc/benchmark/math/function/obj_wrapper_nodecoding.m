function f = obj_wrapper_nodecoding(fitnessfun, x)
% obj_wrapper_nodecoding returns the objective function value f of 
% fitnessfun.
% 
% [SYNTAX]
% f = obj_wrapper_nodecoding(fitnessfun, x)
% 
% [INPUT]
% fitnessfun :  Function handle for a fitness function.
% x          :  Decision variables.
% 
% [OUTPUT]
% f          :  Objective function value.


f = fitnessfun(x);
