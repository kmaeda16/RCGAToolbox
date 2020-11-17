function writeBestAlt(elapsedTime, generation, problem, opts, x, neval)
% writeBestAlt makes an output file, the best individual file.
% 
% [SYNTAX]
% writeBestAlt(elapsedTime, generation, problem, opts, x, neval)
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
% x           :  Decision variables of the the best individual.
% neval       :  Number of evaluation of fitnessfun.


n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
fitnessfun = problem.fitnessfun;
out_best = opts.out_best;


if strcmpi('none',out_best)
    return;
end

out = fopen(out_best,'w');
if out == -1
    warning('cannot open %s!\n',out_best);
    return;
end


%% Get f, g, and phi
if n_constraint == 0
    if nargout(fitnessfun) > 1
        warning('n_constraint was set to %d, but %s returns g.',n_constraint,func2str(fitnessfun));
    end
    f = feval(fitnessfun,x);
    g = 0;
    phi = 0;
else
    if nargout(fitnessfun) == 1
        error('n_constraint was set to %d, but %s does not returns g.',n_constraint,func2str(fitnessfun));
    end
    [f, g] = feval(fitnessfun,x);
    phi = sum( max(0,g) .^ 2 );
    length_g = length(g);
    if n_constraint ~= length_g
        error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
    end
end


%% Write
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

fprintf(out,'%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,f);

if n_constraint > 0
    fprintf(out,'%e\t',phi);
end

for i = 1 : n_gene
    fprintf(out,'%e\t',x(i));
end

for i = 1 : n_constraint
    fprintf(out,'%e\t',g(i));
end

fprintf(out,'\n');

fclose(out);
