function writeBest(elapsedTime, generation, problem, opts, x)
% RCGAwriteBest make a elapsed time, generation, fitness, and phi of the
% best individual.
% 
% [SYNTAX]
% RCGAwriteBest(elapsedTime, generation, problem, opts, best)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec)
% generation  :  Generation
% problem     :  Problem structure
% best        :  Structure of the the best individual

n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
fitnessfun = problem.fitnessfun;
out_best = opts.out_best;
% n_population = opts.n_population;
% n_children = opts.n_children;

if strcmpi('none',out_best)
    return;
end

out = fopen(out_best,'w');
if out == -1
    warning('cannot open %s!\n',out_best);
    return;
end

%% Getting f, g, and phi
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

%%
% fprintf(out,'Time\tNEval\tGeneration\tf\t');
fprintf(out,'Time\tGeneration\tf\t');
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

% x = decodingfun(best.gene);
% neval = n_population + ( generation - 1 ) * n_children;

% fprintf(out,'%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,best.f);
fprintf(out,'%e\t%e\t%e\t',elapsedTime,generation,f);
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



% function writeBest_wrapper(elapsedTime, n_gene, n_constraint, fitnessfun, decodingfun, out_best, gene, generation)
% 
% problem.n_gene = n_gene;
% problem.n_constraint = n_constraint;
% problem.decodingfun = decodingfun;
% opts.out_best = out_best;
% opts.n_population = nan;
% opts.n_children = nan;
% 
% best.gene = gene;
% ObjectiveFunction  = @(gene) obj_wrapper(fitnessfun, decodingfun, gene);
% best.f = ObjectiveFunction(gene);
% if n_constraint > 0
%     ConstraintFunction = @(gene) cst_wrapper(fitnessfun, decodingfun, gene);
%     best.g = ConstraintFunction(gene);
% else
%     best.g = 0;
% end
% best.phi = sum( max(0,best.g) .^2 );
% 
% RCGAwriteBest(elapsedTime, generation, problem, opts, best)
