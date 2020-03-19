function writeBest(elapsedTime, generation, problem, opts, x, neval)
% writeBest makes an output file which includes elapsed time, generation,
% fitness, and phi of the best individual.
% 
% [SYNTAX]
% writeBest(elapsedTime, generation, problem, opts, x, neval)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec)
% generation  :  Generation
% problem     :  Problem structure
% opts        :  RCGA options. See XXXXXXXXXXX for options.
% x           :  Decision variables
% neval       :  Number of evaluation of fitnessfun

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

% x = decodingfun(best.gene);
% neval = n_population + ( generation - 1 ) * n_children;

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
