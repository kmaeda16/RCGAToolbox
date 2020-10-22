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
%               - opts.n_generation: Number of maximum generations.
%               - opts.maxtime: Maximum time (sec).
%               - opts.maxeval: Maximum number of fitnessfun evaluations.
%               - opts.vtr: Value to be reached.
%               - opts.n_par: Number of workers in parallel computation.
%               - opts.output_intvl: Interval generation for updating the 
%                  transition file and the report file.
%               - opts.out_transition: Name of an output file called the 
%                  transition file.
%               - opts.out_best: Name of an output file called the best 
%                  individual file.
%               - opts.out_population: Name of an output file called the 
%                  final population file.
%               - opts.out_report: Name of an output file called the report
%                  file.
%               - opts.interimreportfun: Function handle for the interim 
%                  report function.
%               - opts.finalreportfun: Function handle for the final report 
%                  function.
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
    for j = 1 : n_gene
        Population(i).gene(j) = rand();
    end
end

Population = RCGArequestFitnessCalc(problem,opts,Population);


%% Sorting population
Population = RCGAsrsort(Population,Pf);
