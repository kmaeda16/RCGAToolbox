function x = Rosenbrock_star_decode(gene)
% Decoding function for Rosenbrock_star.
% 
% [SYNTAX]
% x = Rosenbrock_star_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -2.048;
ub =  2.048;

x = gene * ( ub - lb ) + lb;
