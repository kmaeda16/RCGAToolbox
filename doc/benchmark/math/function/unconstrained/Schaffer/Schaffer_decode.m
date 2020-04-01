function x = Schaffer_decode(gene)
% Decoding function for Schaffer.
% 
% [SYNTAX]
% x = Schaffer_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -100;
ub =  100;

x = gene * ( ub - lb ) + lb;
