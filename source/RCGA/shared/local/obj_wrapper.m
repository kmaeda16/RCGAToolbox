function f = obj_wrapper(fitnessfun, decodingfun, gene)
% obj_wrapper returns the objective function value f of fitnessfun.
% 
% [SYNTAX]
% f = obj_wrapper(fitnessfun, decodingfun, gene)
% 
% [INPUT]
% fitnessfun  :  Function handle for a fitness function.
% decodingfun :  Function handle for a decoding function.
% gene        :  Genes (vector of encoded decision variables).
% 
% [OUTPUT]
% f           :  Objective function value.

x = decodingfun(gene);
f = fitnessfun(x);
