function x = g11_decode(gene)
% Decoding function for g11.
% 
% [SYNTAX]
% x = g11_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (2 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (2 dimensional)


lb = -1;
ub =  1;

x = gene * ( ub - lb ) + lb;
