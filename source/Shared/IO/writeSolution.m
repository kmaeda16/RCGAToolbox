function writeSolution(elapsedTime, Param, chrom)

n_gene = Param.n_gene;
n_constraint = Param.n_constraint;
out_solution = Param.out_solution;
decodingfun = Param.decodingfun;

if  strcmp('NONE',out_solution) == 1 ...
        || strcmp('None',out_solution) == 1 ...
        || strcmp('none',out_solution) == 1
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
