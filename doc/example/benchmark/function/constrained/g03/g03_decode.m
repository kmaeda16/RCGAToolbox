function x = g03_decode(gene)
% Decoding function for g03.
% 
% [SYNTAX]
% x = g03_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (10 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (10 dimensional)


lb = 0;
ub = 1;

x = gene * ( ub - lb ) + lb;
