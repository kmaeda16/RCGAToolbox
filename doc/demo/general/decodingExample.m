function x = decodingExample(gene)
% decodingExample is an example of decoding function. "gene" takes values
% from 0 to 1. The purpose of decoding functions is to change the value
% range, i.e. to decode "gene" and return it as x.
% 
% [SYNTAX]
% x = decodingExample(gene)
% 
% [INPUT]
% gene : Encoded decision variables.
% 
% [OUTPUT]
% x    : Decoded decision variables.
% 
% --------------------- Example Problem ---------------------
% Minimize:
%   f = x(1)^2 + x(2)^2 + ... + x(10)^2
% 	
% Subject to:
%   g(1) = x(1) * x(2) + 1 <= 0
%   g(2) = x(1) + x(2) + 1 <= 0
% 	-5.12 <= x(i) <= 5.12 for all i
% -----------------------------------------------------------


lb = -5.12;
ub =  5.12;

x = gene * ( ub - lb ) + lb;
