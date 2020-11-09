function Population = GenerationAlternation_Example(problem, opts, Population)
% GenerationAlternation_Example updates population by a simple generation
% alternation algorithm. This function can be used as a template for custom
% generation alternation functions.
% 
% [SYNTAX]
% Population = GenerationAlternation_Example(problem, opts, Population)
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
% Population :  Population (Array of individuals).
% 
% [OUTPUT]
% Population :  Updated population (Array of individuals).


%% Shortening variable names
n_population = opts.n_population;
n_children = n_population - 1; % Note that opts.n_children is not used.
n_constraint = problem.n_constraint;
n_gene = problem.n_gene;
Pf = opts.Pf;


%% Generating new children
ip = randperm(n_population,2);

c(1,1:n_children) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',zeros,'phi',zeros);

for i = 1 : n_children
    
    % Crossover
    c(i).gene = 0.5 * ( Population(ip(1)).gene + Population(ip(2)).gene );
    
    % Mutation
    mutation_rate = 0.1;
    for j = 1 : n_gene
        if rand < mutation_rate
            c(i).gene(j) = rand;
        end
    end

end
c = RCGArequestFitnessCalc(problem,opts,c);


%% Making a new population
Population = [ Population(1) c ];
Population = RCGAsrsort(Population,Pf);
