% This script demonstrates how to run a user-defined custom real-coded
% genetic algorithm to solve an example constrained optimization problem.
% 
% --------------------- Example Problem ---------------------
% Minimize:
%   f = x(1)^2 + x(2)^2 + ... + x(10)^2
% 
% Subject to:
%   g(1) = x(1) * x(2) + 1 <= 0
%   g(2) = x(1) + x(2) + 1 <= 0
% 	-5.12 <= x(i) <= 5.12 for all i
% 
% 
% Global minimum is f = 3, g(1) = 0, g(2) = 0 at x = (-1.618, 0.6180, 0, 0,
% 0, 0, 0, 0, 0, 0) or at x = (0.6180, -1.618, 0, 0, 0, 0, 0, 0, 0, 0)
% -----------------------------------------------------------


clearvars;

% ========= Problem Settings ========= %
problem.n_gene = 10; % # Decision Variables
problem.n_constraint = 2; % # Constraints
problem.fitnessfun = @Fitness_Example; % Fitness Function
problem.decodingfun = @Decoding_Example; % Decoding Function

% ========== Executing RCGA ========== %
% GenerationAlternation_Example is a user-defined custom generation
% alternation function
Results = RCGA_CustomRCGA(problem, @GenerationAlternation_Example);
