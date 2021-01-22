function x = RotationMatrix_decode(gene)
% Decoding function for RotationMatrix.
% 
% [SYNTAX]
% x = RotationMatrix_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -1;
ub =  1;

x = gene * ( ub - lb ) + lb;
