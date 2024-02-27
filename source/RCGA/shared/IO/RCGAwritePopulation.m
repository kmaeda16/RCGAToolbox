function RCGAwritePopulation(problem, opts, Population)
% RCGAwritePopulation makes an output file, the final population file.
% 
% [SYNTAX]
% RCGAwritePopulation(problem, opts, Population)
% 
% [INPUT]
% problem    :  Problem structure:
%               - problem.n_gene: Number of decision variables.
%               - problem.n_constraint: Number of constraint functions. For 
%                  unconstained problems, this must be zero.
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
% Population :  Array of individual (chrom structure).


%% Shortening variable names
decodingfun = problem.decodingfun;
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
n_population = opts.n_population;
out_population = opts.out_population;


%% If out_population is 'none', nothing done.
if isempty(out_population) || strcmpi('none',out_population)
    return;
end


%% Opening out_population
[ dirname, ~, ~ ] = fileparts(out_population);
if ~exist(dirname,"dir") && ~isempty(dirname)
    mkdir(dirname);
end

out = fopen(out_population,'w');
if out == -1
    warning('cannot open %s!\n',out_population);
    return;
end

fprintf(out,'f\t');
if n_constraint > 0
    fprintf(out,'phi\t');
end
for j = 1 : n_gene
    fprintf(out,'x(%d)\t',j);
end
for j = 1 : n_constraint
    fprintf(out,'g(%d)\t',j);
end
fprintf(out,'\n');


%% Making outputs
for i = 1 : n_population
    x = decodingfun(Population(i).gene);
    fprintf(out,'%e\t',Population(i).f);
    if n_constraint > 0
        fprintf(out,'%e\t',Population(i).phi);
    end
    for j = 1 : n_gene
        fprintf(out,'%e\t',x(j));
    end
    for j = 1 : n_constraint
        fprintf(out,'%e\t',Population(i).g(j));
    end
    fprintf(out,'\n');
end

fclose(out);
