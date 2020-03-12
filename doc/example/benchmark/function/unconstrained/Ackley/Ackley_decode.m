function x = Ackley_decode(gene)
% Decoding function for Ackley.
% 
% [SYNTAX]
% x = Ackley_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -32.768;
ub =  32.768;

x = gene * ( ub - lb ) + lb;
