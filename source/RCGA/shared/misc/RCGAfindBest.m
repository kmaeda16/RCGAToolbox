function index = RCGAfindBest(Population)
% RCGAfindBest returns index of the best individual in Population
% 
% [SYNTAX]
% index = RCGAfindBest(Population)
% 
% [INPUT]
% Population :  Array of individuals
% 
% [OUTPUT]
% index      :  Index of the best individual in Population


%% Preparation
n_population = length(Population);
f = Inf;
phi = Inf;
index = 0;


%% Finding the best individual
for i = 1 : n_population
    if Population(i).phi < phi ...
            || ( Population(i).phi == phi && Population(i).f < f )
        f = Population(i).f;
        phi = Population(i).phi;
        index = i;
    end
end


%% If f or phi are inf or nan, it means fitness calculation went wrong
for i = 1 : n_population
    if isnan(Population(i).f) || isnan(Population(i).phi) ...
            || isinf(Population(i).f) || isinf(Population(i).phi)
        error('Population(%d) has inf or nan for f or phi!',i);
    end
end
