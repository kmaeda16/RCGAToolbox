function [a, b] = RCGAswap(a, b)
% RCGAswap receives two input variables, swap them, and return them. The
% inputs can be anything, e.g., vectors, structures, and objects.
% 
% [SYNTAX]
% [a, b] = RCGAswap(a, b)
% 
% [INPUT]
% a :  The first input.
% b :  The second input.
% 
% [OUTPUT]
% a :  This equals the second input.
% b :  This equals the first input.

c = a;
a = b;
b = c;
