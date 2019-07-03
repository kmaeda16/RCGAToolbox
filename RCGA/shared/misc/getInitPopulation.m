function Population = getInitPopulation(problem,opts)
% getInitPopulation returns randomly generated initial population.
% 
% [SYNTAX]
% Population = getInitPopulation(problem,opts)
% 
% [INPUT]
% problem :  Problem structure
% opts    :  RCGA options. See XXXXXXXXXXX for options.
% 
% [OUTPUT]
% Population :  Randomly generated initial population


%% Shortening variable names
n_gene = problem.n_gene;
n_population = opts.n_population;
n_constraint = problem.n_constraint;
Pf = opts.Pf;


%% Preparation
Population(1,1:n_population) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',0,'phi',0);


%% Generating population
for i = 1 : n_population
    for j = 1 : n_gene
        Population(i).gene(j) = rand();
    end
end

Population = requestFitnessCalc(problem,opts,Population);


%% Sorting population
Population = SRsort(Population,Pf);
