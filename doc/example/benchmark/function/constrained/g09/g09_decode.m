function x = g09_decode(gene)
% Decoding function for g09.
% 
% [SYNTAX]
% x = g09_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (7 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (7 dimensional)


lb = -10;
ub =  10;

x = gene * ( ub - lb ) + lb;
