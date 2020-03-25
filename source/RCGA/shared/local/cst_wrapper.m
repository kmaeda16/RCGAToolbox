function [c, ceq] = cst_wrapper(fitnessfun, decodingfun, gene)
% cst_wrapper returns g of fitnessfun.
% 
% [SYNTAX]
% [c, ceq] = cst_wrapper(fitnessfun, decodingfun, gene)
% 
% [INPUT]
% fitnessfun  :  Fitness function
% decodingfun :  Decoding function
% gene        :  Genes (encoded decision variables)
% 
% [OUTPUT]
% c           :  Inequality constraint function values
% ced         :  Equality constraint function values

x = decodingfun(gene);
[~, c] = fitnessfun(x);
ceq = [];
