function x = g01_decode(gene)
% Decoding function for g01.
% 
% [SYNTAX]
% x = g01_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables (13 dimensional)
% 
% [OUTPUT]
% x    : Decoded decision variables (13 dimensional)


lb = zeros(1,13);
ub = zeros(1,13);

lb(1:9) = 0; 
ub(1:9) = 1;

lb(10:12) = 0;
ub(10:12) = 100;

lb(13) = 0;
ub(13) = 1;

x = gene .* ( ub - lb ) + lb;
