function x = Rosenbrock_chain_decode(gene)
% Decoding function for Rosenbrock_chain.
% 
% [SYNTAX]
% x = Rosenbrock_chain_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -2.048;
ub =  2.048;

x = gene * ( ub - lb ) + lb;
