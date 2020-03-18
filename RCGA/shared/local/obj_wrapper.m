function f = obj_wrapper(fitnessfun, decodingfun, gene)
% obj_wrapper returns f of fitnessfun.
% 
% [SYNTAX]
% f = obj_wrapper(fitnessfun, decodingfun, gene)
% 
% [INPUT]
% fitnessfun  :  Fitness function
% decodingfun :  Decoding function
% gene        :  Genes (encoded decision variables)
% 
% [OUTPUT]
% f           :  Objective function value

x = decodingfun(gene);
f = fitnessfun(x);
