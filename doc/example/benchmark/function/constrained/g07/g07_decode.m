function x = g07_decode(gene)
% Decoding function for g07.
% 
% [SYNTAX]
% x = g07_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (10 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (10 dimensional)


lb = -10;
ub =  10;

x = gene * ( ub - lb ) + lb;
