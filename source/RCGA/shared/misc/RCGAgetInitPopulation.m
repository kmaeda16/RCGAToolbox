function Population = RCGAgetInitPopulation(problem, opts)
% RCGAgetInitPopulation returns randomly generated initial population.
% 
% [SYNTAX]
% Population = RCGAgetInitPopulation(problem, opts)
% 
% [INPUT]
% problem    :  Problem structure:
%               - problem.n_gene: Number of decision variables.
%               - problem.n_constraint: Number of constraint functions. 
%                  For unconstained problems, this must be zero.
%               - problem.fitnessfun: Function handle for a fitness 
%                  function.
%               - problem.decodingfun: Function handle for a decoding 
%                  function.
% opts       :  Option structure:
%               - opts.n_population: Population size.
%               - opts.n_children: Number of children.
%               - opts.n_parent: Number of parents.
%               - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%               - opts.selection_type: Selection type for REXstar/JGG 
%                  (0 or 1).
%               - opts.Pf: Probability that only the objective function f 
%                  is used in comparisons of individuals in the stochastic
%                  ranking.
%               - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                  local optimizer is used.
%               - opts.localopts: Options for the local optimizer.
%               - opts.maxgen: Maximum number of generations.
%               - opts.maxtime: Maximum time (sec).
%               - opts.maxeval: Maximum number of fitnessfun evaluations.
%               - opts.maxstallgen: Maximum number of stall generations for
%                  early stopping.
%               - opts.vtr: Value to be reached.
%               - opts.n_par: Number of workers in parallel computation.
%               - opts.initial_population: n x n_gene matrix in which each
%                  row represents an individual. Note that each gene 
%                  should be 0 ~ 1. The first n_population individuals of 
%                  the designated initial population are used, and others 
%                  are ignored. If n < n_population, 
%                  n_population - n individuals are randomly generated. 
%               - opts.output_intvl: Interval generation for updating the 
%                  transition file and the report file.
%               - opts.out_transition: Name of an output file called the 
%                  transition file.
%               - opts.out_best: Name of an output file called the best 
%                  individual file.
%               - opts.out_population: Name of an output file called the 
%                  final population file.
%               - opts.out_report: Name of an output file called the 
%                  report file.
%               - opts.interimreportfun: Function handle for the interim 
%                  report function.
%               - opts.finalreportfun: Function handle for the final 
%                  report function.
% 
% [OUTPUT]
% Population :  Randomly generated initial population (Array of
%               individuals).


%% Shortening variable names
n_gene = problem.n_gene;
n_population = opts.n_population;
n_constraint = problem.n_constraint;
Pf = opts.Pf;


%% Preparation
Population(1,1:n_population) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',0,'phi',0);


%% Generating population
for i = 1 : n_population
    Population(i).gene = rand(1,n_gene);
end


%% Setting designated initial population
if isfield(opts,'initial_population') && ~isempty(opts.initial_population)
    [ n_designated, ~] = size(opts.initial_population);
    for i = 1 : min(n_population,n_designated)
        Population(i).gene = opts.initial_population(i,1:n_gene);
        Population(i).gene(Population(i).gene>1) = 1;
        Population(i).gene(Population(i).gene<0) = 0;
    end
end


%% Fitness calculation
Population = RCGArequestFitnessCalc(problem,opts,Population);


%% Sorting population
Population = RCGAsrsort(Population,Pf);
