function [a, b] = RCGAswap(a, b)
% swap receives two input variables, swap them, and return them.
% 
% [SYNTAX]
% [a, b] = swap(a, b)
% 
% [INPUT]
% a :  The first input
% b :  The second input
% 
% [OUTPUT]
% a :  This equals the second input
% b :  This equals the first input

c = a;
a = b;
b = c;
