function writePopulation(problem, opts, Population)

decodingfun = problem.decodingfun;
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
n_population = opts.n_population;
out_population = opts.out_population;

if  strcmpi('none',out_population)
    return;
end

out = fopen(out_population,'w');
if out == -1
    warning('cannot open %s!\n',out_transition);
    return;
end

fprintf(out,'No\tf\tphi\t');
for j = 1 : n_gene
    fprintf(out,'x[%d]\t',j);
end
for j = 1 : n_constraint
    fprintf(out,'g[%d]\t',j);
end
fprintf(out,'\n');

for i = 1 : n_population
    x = decodingfun(Population(i).gene);
    fprintf(out,'%d\t%e\t%e\t',i,Population(i).f,Population(i).phi);
    for j = 1 : n_gene
        fprintf(out,'%e\t',x(j));
    end
    for j = 1 : n_constraint
        fprintf(out,'%e\t',Population(i).g(j));
    end
    fprintf(out,'\n');
end

fclose(out);
