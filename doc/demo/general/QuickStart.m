% This script demonstrates how to run a real-coded genetic algorithm to
% solve an example constrained optimization problem.
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


clearvars;

% ========= Problem Settings ========= %
problem.n_gene = 10; % # Decision Variables
problem.n_constraint = 2; % # Constraints
problem.fitnessfun = @fitnessExample; % Fitness Function
problem.decodingfun = @decodingExample; % Decoding Function

% ========== Executing RCGA ========== %
% Results = RCGA_UNDXMGG(problem); % UNDX/MGG
Results = RCGA_REXstarJGG(problem); % REXstar/JGG
