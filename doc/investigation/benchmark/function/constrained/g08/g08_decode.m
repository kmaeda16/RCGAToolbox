function x = g08_decode(gene)
% Decoding function for g08.
% 
% [SYNTAX]
% x = g08_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (2 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (2 dimensional)


lb = 0;
ub = 10;

x = gene * ( ub - lb ) + lb;
