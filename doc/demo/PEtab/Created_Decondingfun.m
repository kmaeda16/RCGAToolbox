function param = Created_Decondingfun(gene)
% This function was created by RCGAToolbox RCGAcreateDecodingfun
% based on PEtab_Parameter.tsv.
% Created: 19-Apr-2021

% X0
param(1) = 1.000000e-01;

% k1
lb = 0.000000e+00;
ub = 1.000000e+01;
param(2) = lb + ( ub - lb ) * gene(2);

% k2
lb = 0.000000e+00;
ub = 1.000000e+01;
param(3) = lb + ( ub - lb ) * gene(3);

% k3
lb = log10(1.000000e-01);
ub = log10(1.000000e+01);
param(4) = 10 ^ ( lb + ( ub - lb ) * gene(4) );

% K2
lb = log10(1.000000e-01);
ub = log10(1.000000e+01);
param(5) = 10 ^ ( lb + ( ub - lb ) * gene(5) );

% K3
lb = log10(1.000000e-01);
ub = log10(1.000000e+01);
param(6) = 10 ^ ( lb + ( ub - lb ) * gene(6) );

% rootCompartment
param(7) = 1.000000e+00;

