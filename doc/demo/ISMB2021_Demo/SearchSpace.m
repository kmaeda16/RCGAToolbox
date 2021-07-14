function x = SearchSpace(gene)
% SearchSpace is an example of decoding function for
% "Model_Example_odefun.m". "gene" takes values from 0 to 1. The purpose of
% decoding functions is to change the value range, i.e. to decode "gene"
% and return it as x.
% 
% [SYNTAX]
% x = SearchSpace(gene)
% 
% [INPUT]
% gene : Encoded decision variables.
% 
% [OUTPUT]
% x    : Decoded decision variables.
% 
% === x and corresponding parameters ===
%   x(1) : EmptySet
%   x(2) : kappa
%   x(3) : k6
%   x(4) : k4
%   x(5) : k4prime
%   x(6) : cell


%% Linear search space
x = [1, 0.015, 1, 180, 0.018, 1]; % Default parameter values
lb = 0.33;
ub = 3.0;
x = x .* ( gene * ( ub - lb ) + lb );
