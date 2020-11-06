function x = Decoding_Example(gene)
% Decoding_Example is an example of decoding function for
% "Model_Example_odefun.m". "gene" takes values from 0 to 1. The purpose of
% decoding functions is to change the value range, i.e. to decode "gene"
% and return it as x.
% 
% [SYNTAX]
% x = Decoding_Example(gene)
% 
% [INPUT]
% gene : Encoded decision variables.
% 
% [OUTPUT]
% x    : Decoded decision variables.
% 
% === x and corresponding parameters ===
%   x(1) : X0
%   x(2) : k1
%   x(3) : k2
%   x(4) : k3
%   x(5) : K2
%   x(6) : K3
%   x(7) : rootCompartment


%% Linear search space
lb = 0;
ub = 2;
x = gene * ( ub - lb ) + lb;

%% Logarithmic search space
% lb = -2;
% ub = 2;
% x = 10 .^ ( gene * ( ub - lb ) + lb );

%% You can fix some of parameter values
% x(1) = 0.1;
% x(7) = 1;

%% This is the parameter values to create "measurementExample.xls"
% x = [ 0.1, 1, 1, 1, 1, 1, 1 ];
