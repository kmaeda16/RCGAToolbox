function Population = RCGA_JGG(problem, opts, Population)
% RCGA_JGG updates population by using Just Generation Gap (JGG).
% 
% [SYNTAX]
% Population = RCGA_JGG(problem, opts, Population)
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
% Population :  Population (Array of individuals).
% 
% [OUTPUT]
% Population :  Updated population (Array of individuals).
% 
% 
% See "Kimura S, Sato M, Okada-Hatakeyama M: An Effective Method for the
% Inference of Reduced S-system Models of Genetic Networks. Information and
% Media Technologies 2015, 10(1):166-174.".


%% Shortening variable names
n_population = opts.n_population;
n_parent = opts.n_parent;
Pf = opts.Pf;
selection_type = opts.selection_type;


%% Pick up parents from main population
ip = randperm(n_population,n_parent);


%% Generating children
p = Population(ip);
c = RCGA_REXstar(problem,opts,p);


%% Updating population
switch selection_type
    case 0
        % Chosen from children (Kobayashi, 2009)
        Population(ip) = c(1:n_parent);
    case 1
        % Chosen from family (Kimura et al., 2015)
        f = [p c];
        f = RCGAsrsort(f,Pf);
        Population(ip) = f(1:n_parent);
    otherwise
        error('Unexpected selection_type!');
end


%% Stochastic ranking sort
Population = RCGAsrsort(Population,Pf);
