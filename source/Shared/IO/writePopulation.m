function writePopulation(Param, Population)

out_population = Param.out_population;
decodingfun = Param.decodingfun;
n_population = Param.n_population;
n_gene = Param.n_gene;
n_constraint = Param.n_constraint;

if  strcmp('NONE',out_population) == 1 ...
        || strcmp('None',out_population) == 1 ...
        || strcmp('none',out_population) == 1
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
