function dydt = Model_Example_conciseOdefun(t, y)
% Model_Example_conciseOdefun is an example of the concise ODE file
% (RCGAToolbox format) ready for conversion into an IQMmodel object by
% RCGAreadConciseODEfile. In this file, the only sections sandwiched
% between "BEGIN" and "END" are used for the conversion.
% 
% [SYNTAX]
% dydt = Model_Example(t, y)
% 
% [INPUT]
% t    :  Time.
% y    :  X1 and X2.
% 
% [OUTPUT]
% dydt :  dX1/dt and dX2/dt.


X1 = y(1);
X2 = y(2);

% === BEGIN NAME ===
% Model_Example
% === END NAME ===

% === BEGIN NOTES ===
% Simple metabolic pathway with two Michaelis-Menten rate equations.
% === END NOTES ===

% === BEGIN INITIAL CONDITION ===
X1_0 = 0;
X2_0 = 0;
% === END INITIAL CONDITION ===

% === BEGIN PARAMETERS ===
X0 = 0.1;
k1 = 1;
k2 = 1;
k3 = 1;
K2 = 1;
K3 = 1;
% === END PARAMETERS ===

% === BEGIN VARIABLES ===
X12 = X1 + X2;
% === END VARIABLES ===

% === BEGIN REACTIONS ===
v1 = k1 * X0;
v2 = k2 * X1 / ( K2 + X1 );
v3 = k3 * X2 / ( K3 + X2 );
% === END REACTIONS ===

% === BEGIN BALANCE ===
X1_dot = v1 - v2;
X2_dot = v2 - v3;
% === END BALANCE ===

dydt = zeros(2,1);
dydt(1) = X1_dot;
dydt(2) = X2_dot;
