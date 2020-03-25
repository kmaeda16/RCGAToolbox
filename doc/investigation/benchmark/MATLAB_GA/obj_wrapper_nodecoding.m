function f = obj_wrapper_nodecoding(fitnessfun, x)
% obj_wrapper_nodecoding returns f of fitnessfun.
% 
% [SYNTAX]
% f = obj_wrapper(fitnessfun, decodingfun, gene)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f :  Objective function value

% x = decodingfun(gene);
f = fitnessfun(x);
