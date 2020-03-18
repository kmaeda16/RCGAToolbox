function RCGAwriteBest(elapsedTime, generation, problem, opts, best)
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
% opts        :  RCGA options. See XXXXXXXXXXX for options.
% best        :  Structure of the the best individual


n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
decodingfun = problem.decodingfun;
out_best = opts.out_best;
n_population = opts.n_population;
n_children = opts.n_children;

if strcmpi('none',out_best)
    return;
end

out = fopen(out_best,'w');
if out == -1
    warning('cannot open %s!\n',out_best);
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

x = decodingfun(best.gene);
global RCGA_LOCALNEVAL;
neval = n_population + ( generation - 1 ) * n_children + RCGA_LOCALNEVAL;

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
