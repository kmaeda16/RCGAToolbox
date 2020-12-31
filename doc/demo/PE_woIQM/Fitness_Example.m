function ssr = Fitness_Example(x)
% Fitness_Example is an example of fitness function for parameter 
% estimation.
% 
% [SYNTAX]
% ssr = Fitness_Example(x)
% 
% [INPUT]
% x   : Decision variables.
% 
% [OUTPUT]
% ssr : Objective function to be minimized. In this example, it is sum of 
%       squared resudials (SSR).
% 
% ------------------------ Example Kinetic Model ------------------------
% - INITIAL CONDITION
% X1 = 0
% X2 = 0
% 
% - PARAMETERS
% X0 = 0.1
% k1 = 1
% k2 = 1
% k3 = 1
% K2 = 1
% K3 = 1
% rootCompartment = 1
% 
% - VARIABLES
% X12 = X1 + X2
% 
% - REACTIONS
% v1 = k1 * X0
% v2 = k2 * (X1/rootCompartment) / (K2 + (X1/rootCompartment))
% v3 = k3 * (X2/rootCompartment) / (K3 + (X2/rootCompartment))
% 
% - BALANCE
% X1_dot = v1 - v2;
% X2_dot = v2 - v3;
% -----------------------------------------------------------------------


global experimentaldata;


%% Checking time errors
tspan = experimentaldata(:,1);
y0 = [0 0];
param = x;


%% Running simulation
[ ~, Y ] = RCGAsimulateODE(@Model_Example_odefun, tspan, y0, param);


%% If simulation faild, return ssr = 1e+10
[n_row, ~] = size(Y);
if max(max(isnan(Y))) || length(tspan) ~= n_row
    warning('Simulation failed!');
    ssr = 1e+10;
    return;
end


%% Preparing Y_sim and Y_exp
Y_sim = real(Y);
[n_row, n_col] = size(Y_sim);
Y_exp = zeros([n_row n_col]);
for i = 1 : n_col
    Y_exp(:,i) = experimentaldata(:,1+i);
end


%% Calculating SSR
ssr = 0;
for i = 1 : n_row
    for j = 1 : n_col
        if ~isnan(Y_exp(i,j))
            % ssr = ssr + abs( Y_sim(i,j) - Y_exp(i,j) );
            % ssr = ssr + abs( ( Y_sim(i,j) - Y_exp(i,j) ) / Y_exp(i,j) );
            ssr = ssr + ( Y_sim(i,j) - Y_exp(i,j) ) ^ 2;
            % ssr = ssr + ( ( Y_sim(i,j) - Y_exp(i,j) ) / Y_exp(i,j) ) ^ 2;
        end
    end
end
