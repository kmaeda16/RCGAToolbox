function writeTransition(elapsedTime, generation, Param, chrom)

n_gene = Param.n_gene;
n_constraint = Param.n_constraint;
out_transition = Param.out_transition;
decodingfun = Param.decodingfun;

if  strcmp('NONE',out_transition) == 1 ...
        || strcmp('None',out_transition) == 1 ...
        || strcmp('none',out_transition) == 1
    return;
end
    
if generation == 1
    out = fopen(out_transition,'w');
    if out == -1
        warning('cannot open %s!\n',out_transition);
        return;
    end
    fprintf(out,'Time\tGeneration\tf\tphi\t');
    for i = 1 : n_gene
        fprintf(out,'x[%d]\t',i);
    end
    for i = 1 : n_constraint
        fprintf(out,'g[%d]\t',i);
    end
    fprintf(out,'\n');
else
    out = fopen(out_transition,'a');
    if out == -1
        warning('cannot open %s!\n',out_transition);
        return;
    end
end

x = decodingfun(chrom.gene);
    
fprintf(out,'%e\t%d\t%e\t%e\t',elapsedTime,generation,chrom.f,chrom.phi);
for i = 1 : n_gene
    fprintf(out,'%e\t',x(i));
end
for i = 1 : n_constraint
    fprintf(out,'%e\t',chrom.g(i));
end
fprintf(out,'\n');

fclose(out);
