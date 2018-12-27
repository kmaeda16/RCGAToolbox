function writeSolution(elapsedTime, problem, opts, chrom)

n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
decodingfun = problem.decodingfun;
out_solution = opts.out_solution;

if  strcmpi('none',out_solution)
    return;
end

out = fopen(out_solution,'w');
if out == -1
    warning('cannot open %s!\n',out_solution);
    return;
end

fprintf(out,'Time\tf\tphi\t');
for i = 1 : n_gene
    fprintf(out,'x[%d]\t',i);
end
for i = 1 : n_constraint
    fprintf(out,'g[%d]\t',i);
end
fprintf(out,'\n');

x = decodingfun(chrom.gene);

fprintf(out,'%e\t%e\t%e\t',elapsedTime,chrom.f,chrom.phi);
for i = 1 : n_gene
    fprintf(out,'%e\t',x(i));
end
for i = 1 : n_constraint
    fprintf(out,'%e\t',chrom.g(i));
end
fprintf(out,'\n');

fclose(out);
