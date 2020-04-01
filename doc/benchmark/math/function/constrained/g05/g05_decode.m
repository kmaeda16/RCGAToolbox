function x = g05_decode(gene)
% Decoding function for g05.
% 
% [SYNTAX]
% x = g05_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (4 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (4 dimensional)


lb = zeros(1,4);
ub = zeros(1,4);

lb(1:2) = 0; 
ub(1:2) = 1200;

lb(3:4) = -0.55;
ub(3:4) =  0.55;

x = gene .* ( ub - lb ) + lb;
