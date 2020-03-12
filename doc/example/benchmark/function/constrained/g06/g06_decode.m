function x = g06_decode(gene)
% Decoding function for g06.
% 
% [SYNTAX]
% x = g06_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (2 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (2 dimensional)


lb = zeros(1,2);
ub = zeros(1,2);

lb(1) = 13; 
ub(1) = 100;

lb(2) = 0;
ub(2) = 100;

x = gene .* ( ub - lb ) + lb;
