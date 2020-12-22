function [improvedChrom, localneval] = RCGAlocalOptimize(problem, opts, chrom)
% RCGAlocalOptimize improves an individual (chrom) by the local optimizer 
% fmincon. Optiimzation Toolbox is required to use fmincon.
% 
% [SYNTAX]
% [improvedChrom, localneval] = RCGAlocalOptimize(problem, opts, chrom)
% 
% [INPUT]
% problem       :  Problem structure:
%                  - problem.n_gene: Number of decision variables.
%                  - problem.n_constraint: Number of constraint functions. 
%                     For unconstained problems, this must be zero.
%                  - problem.fitnessfun: Function handle for a fitness 
%                     function.
%                  - problem.decodingfun: Function handle for a decoding 
%                     function.
% opts          :  Option structure:
%                  - opts.n_population: Population size.
%                  - opts.n_children: Number of children.
%                  - opts.n_parent: Number of parents.
%                  - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%                  - opts.selection_type: Selection type for REXstar/JGG 
%                     (0 or 1).
%                  - opts.Pf: Probability that only the objective function 
%                     f is used in comparisons of individuals in the 
%                     stochastic ranking.
%                  - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                     local optimizer is used.
%                  - opts.localopts: Options for the local optimizer.
%                  - opts.maxgen: Maximum number of generations.
%                  - opts.maxtime: Maximum time (sec).
%                  - opts.maxeval: Maximum number of fitnessfun 
%                     evaluations.
%                  - opts.vtr: Value to be reached.
%                  - opts.n_par: Number of workers in parallel computation.
%                  - opts.output_intvl: Interval generation for updating 
%                     the transition file and the report file.
%                  - opts.out_transition: Name of an output file called the 
%                     transition file.
%                  - opts.out_best: Name of an output file called the best 
%                     individual file.
%                  - opts.out_population: Name of an output file called the 
%                     final population file.
%                  - opts.out_report: Name of an output file called the 
%                     report file.
%                  - opts.interimreportfun: Function handle for the interim 
%                     report function.
%                  - opts.finalreportfun: Function handle for the final 
%                     report function.
% chrom         :  Individual (chrom structure).
% 
% [OUTPUT]
% improvedChrom :  Improved individual (chrom structure).


localopts = opts.localopts;

LB = zeros(1,problem.n_gene); % Lower bound
UB = ones(1,problem.n_gene);  % Upper bound
    
ObjectiveFunction  = @(gene) obj_wrapper(problem.fitnessfun, problem.decodingfun, gene);
ConstraintFunction = @(gene) cst_wrapper(problem.fitnessfun, problem.decodingfun, gene);

try
    if problem.n_constraint == 0
        [improvedChrom.gene, ~, ~, output] = fmincon(ObjectiveFunction,chrom.gene,[],[],[],[],LB,UB,[],localopts);
    else
        [improvedChrom.gene, ~, ~, output] = fmincon(ObjectiveFunction,chrom.gene,[],[],[],[],LB,UB,ConstraintFunction,localopts);
    end
catch ME
    warning('Local optimization failed in RCGAlocalOptimize: %s',ME.message);
    improvedChrom = chrom;
    localneval = 1;
    return;
end

[ improvedChrom.f, improvedChrom.g, improvedChrom.phi ] = RCGAgetFitness(problem,improvedChrom);

% If improvedChrom is worse than chrom, then chrom is returned.
if ~( improvedChrom.phi < chrom.phi || ( improvedChrom.phi == chrom.phi && improvedChrom.f < chrom.f ) )
    improvedChrom = chrom;
end

localneval = output.funcCount + 1;
