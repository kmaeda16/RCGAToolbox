function x = MMbenchmark_decode(gene)
% Decoding function for MMbenchmark.
% 
% [SYNTAX]
% x = MMbenchmark_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -1;
ub =  1;

x = gene * ( ub - lb ) + lb;
