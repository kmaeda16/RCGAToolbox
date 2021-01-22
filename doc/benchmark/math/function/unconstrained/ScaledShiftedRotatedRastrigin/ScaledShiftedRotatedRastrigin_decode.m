function x = ScaledShiftedRotatedRastrigin_decode(gene)
% Decoding function for ScaledShiftedRotatedRastrigin_decode.
% 
% [SYNTAX]
% x = ScaledShiftedRotatedRastrigin_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -5.12;
ub =  5.12;

x = gene * ( ub - lb ) + lb;
