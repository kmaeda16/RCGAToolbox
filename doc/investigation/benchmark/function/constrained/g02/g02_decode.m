function x = g02_decode(gene)
% Decoding function for g02.
% 
% [SYNTAX]
% x = g02_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (20 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (20 dimensional)


lb = 0;
ub = 10;

x = gene * ( ub - lb ) + lb;
