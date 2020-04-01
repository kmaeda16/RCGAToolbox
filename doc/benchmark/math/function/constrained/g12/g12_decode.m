function x = g12_decode(gene)
% Decoding function for g12.
% 
% [SYNTAX]
% x = g12_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (3 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (3 dimensional)


lb = 0;
ub = 10;

x = gene * ( ub - lb ) + lb;
