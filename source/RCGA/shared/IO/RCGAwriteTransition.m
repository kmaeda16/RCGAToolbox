function RCGAwriteTransition(elapsedTime, generation, problem, opts, best)
% RCGAwriteTransition makes an output file, the transition file.
% 
% [SYNTAX]
% RCGAwriteTransition(elapsedTime, generation, problem, opts, best)
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
%                - opts.n_generation: Number of maximum generations.
%                - opts.maxtime: Maximum time (sec).
%                - opts.maxeval: Maximum number of fitnessfun evaluations.
%                - opts.vtr: Value to be reached.
%                - opts.n_par: Number of workers in parallel computation.
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
% best        :  Structure of the the best individual.
%                - best.gene: Decision variable vector.
%                - best.g: Constraint function value vector.
%                - best.f: Fitness function value.
%                - best.phi: Penalty function value.


%% Shortening variable names
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
decodingfun = problem.decodingfun;
out_transition = opts.out_transition;
n_population = opts.n_population;
n_children = opts.n_children;


%% If out_transition is 'none', nothing done.
if  isempty(out_transition) || strcmpi('none',out_transition)
    return;
end


%% Opening out_transition
if generation == 1
    out = fopen(out_transition,'w');
    if out == -1
        warning('cannot open %s!\n',out_transition);
        return;
    end
    fprintf(out,'Time\tNEval\tGeneration\tf\t');
    if n_constraint > 0
        fprintf(out,'phi\t');
    end
    for i = 1 : n_gene
        fprintf(out,'x(%d)\t',i);
    end
    for i = 1 : n_constraint
        fprintf(out,'g(%d)\t',i);
    end
    fprintf(out,'\n');
else
    out = fopen(out_transition,'a');
    if out == -1
        warning('cannot open %s!\n',out_transition);
        return;
    end
end


%% Decoding genes to x
x = decodingfun(best.gene);
length_x = length(x);

if length_x ~= n_gene
    error('decodingfun should return a vector with %d elements but it returned a vector with %d elements.',n_gene,length_x);
end


%% Calculating the number of fitness function evaluations
global RCGA_LOCALNEVAL;
neval = n_population + ( generation - 1 ) * n_children + RCGA_LOCALNEVAL;


%% Making outputs
fprintf(out,'%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,best.f);
if n_constraint > 0
    fprintf(out,'%e\t',best.phi);
end
for i = 1 : n_gene
    fprintf(out,'%e\t',x(i));
end
for i = 1 : n_constraint
    fprintf(out,'%e\t',best.g(i));
end
fprintf(out,'\n');

fclose(out);
