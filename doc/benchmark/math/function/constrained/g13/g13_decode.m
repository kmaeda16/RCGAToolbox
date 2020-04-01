function x = g13_decode(gene)
% Decoding function for g13.
% 
% [SYNTAX]
% x = g13_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (5 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (5 dimensional)


lb = zeros(1,5);
ub = zeros(1,5);

lb(1:2) = -2.3;
ub(1:2) =  2.3;

lb(3:5) = -3.2;
ub(3:5) =  3.2;

x = gene .* ( ub - lb ) + lb;
