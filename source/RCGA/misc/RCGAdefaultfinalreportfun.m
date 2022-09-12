function RCGAdefaultfinalreportfun(elapsedTime, generation, problem, opts, Population, best)
% RCGAdefaultfinalreportfun is called at the end of RCGA, making two output
% files, the best individual file and the final population file.
% 
% [SYNTAX]
% RCGAdefaultfinalreportfun(elapsedTime, generation, problem, opts, Population, best)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec).
% generation  :  Generation.
% problem     :  Problem structure:
%                - problem.n_gene: Number of decision variables.
%                - problem.n_constraint: Number of constraint functions. 
%                   For unconstained problems, this must be zero.
%                - problem.fitnessfun: Function handle for a fitness 
%                   function.
%                - problem.decodingfun: Function handle for a decoding 
%                   function.
% opts        :  Option structure:
%                - opts.n_population: Population size.
%                - opts.n_children: Number of children.
%                - opts.n_parent: Number of parents.
%                - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%                - opts.selection_type: Selection type for REXstar/JGG 
%                   (0 or 1).
%                - opts.Pf: Probability that only the objective function f 
%                   is used in comparisons of individuals in the stochastic
%                   ranking.
%                - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                   local optimizer is used.
%                - opts.localopts: Options for the local optimizer.
%                - opts.maxgen: Maximum number of generations.
%                - opts.maxtime: Maximum time (sec).
%                - opts.maxeval: Maximum number of fitnessfun evaluations.
%                - opts.maxstallgen: Maximum number of stall generations 
%                   for early stopping.
%                - opts.vtr: Value to be reached.
%                - opts.n_par: Number of workers in parallel computation.
%                - opts.initial_population: n x n_gene matrix in which each
%                   row represents an individual. Note that each gene 
%                   should be 0 ~ 1. The first n_population individuals of 
%                   the designated initial population are used, and others 
%                   are ignored. If n < n_population, 
%                   n_population - n individuals are randomly generated. 
%                - opts.output_intvl: Interval generation for updating the 
%                   transition file and the report file.
%                - opts.out_transition: Name of an output file called the 
%                   transition file.
%                - opts.out_best: Name of an output file called the best 
%                   individual file.
%                - opts.out_population: Name of an output file called the 
%                   final population file.
%                - opts.out_report: Name of an output file called the 
%                   report file.
%                - opts.interimreportfun: Function handle for the interim 
%                   report function.
%                - opts.finalreportfun: Function handle for the final 
%                   report function.
% Population  :  Population (Array of individuals).
% best        :  Structure of the the best individual.
%                - best.gene: Decision variable vector.
%                - best.g: Constraint function value vector.
%                - best.f: Fitness function value.
%                - best.phi: Penalty function value.


RCGAwriteBest(elapsedTime,generation,problem,opts,best);
RCGAwritePopulation(problem,opts,Population);
