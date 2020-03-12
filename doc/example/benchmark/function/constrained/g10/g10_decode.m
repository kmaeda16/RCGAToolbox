function x = g10_decode(gene)
% Decoding function for g10.
% 
% [SYNTAX]
% x = g10_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (8 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (8 dimensional)


lb = zeros(1,8);
ub = zeros(1,8);

lb(1) = 2; 
ub(1) = 4;

lb(2:3) = 3;
ub(2:3) = 4;

lb(4:8) = 1;
ub(4:8) = 3;

x = 10 .^ ( gene .* ( ub - lb ) + lb );
