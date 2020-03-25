function x = Cigar_decode(gene)
% Decoding function for Cigar.
% 
% [SYNTAX]
% x = Cigar_decode(gene)
% 
% [INPUT]
% gene : Coded decision variables
% 
% [OUTPUT]
% x    : Decoded decision variables


lb = -5.12;
ub =  5.12;

x = gene * ( ub - lb ) + lb;
