function [c, ceq] = cst_wrapper(fitnessfun, decodingfun, gene)
% cst_wrapper returns the constraint function value vector g of fitnessfun.
% 
% [SYNTAX]
% [c, ceq] = cst_wrapper(fitnessfun, decodingfun, gene)
% 
% [INPUT]
% fitnessfun  :  Function handle for a fitness function.
% decodingfun :  Function handle for a decoding function.
% gene        :  Genes (vector of encoded decision variables).
% 
% [OUTPUT]
% c           :  Inequality constraint function values (g).
% ceq         :  Equality constraint function values (empty).

x = decodingfun(gene);
[~, c] = fitnessfun(x);
ceq = [];
