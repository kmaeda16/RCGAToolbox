function x = Schwefel_decode(gene)
% Decoding function for Schwefel.
% 
% [SYNTAX]
% x = Schwefel_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -512;
ub =  512;

x = gene * ( ub - lb ) + lb;
