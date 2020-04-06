function x = decodingExample(gene)
% decodingExample is an example of decoding function for
% "modelExample.xml". "gene" takes values from 0 to 1. The purpose of
% decoding functions is to change the value range, i.e. to decode "gene"
% and return it as x.
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
% === x and corresponding parameters ===
%   x(1) : S4
%   x(2) : S0
%   x(3) : J1_Vmax
%   x(4) : J1_n
%   x(5) : J1_K
%   x(6) : J2_J2_k
%   x(7) : J3_J3_k
%   x(8) : J0_J0_k
%   x(9) : compart


%% Linear search space
lb = 0;
ub = 10;
x = gene * ( ub - lb ) + lb;

%% Logarithmic search space
% lb = -2;
% ub = 1;
% x = 10 .^ ( gene * ( ub - lb ) + lb );

%% You can fix some of parameter values
% x(1) = 0;
% x(end) = 1;

%% This is the parameter values to create "measurementExample.xls"
% x = [0, 5, 5.5, 4, 0.5, 0.1, 0.1, 0.01, 1];
