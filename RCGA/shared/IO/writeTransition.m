function writeTransition(elapsedTime, generation, problem, opts, best)

n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
decodingfun = problem.decodingfun;
out_transition = opts.out_transition;
n_population = opts.n_population;
n_children = opts.n_children;

if  strcmpi('none',out_transition)
    return;
end
    
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

x = decodingfun(best.gene);
neval = n_population + ( generation - 1 ) * n_children;


length_x = length(x);
if length_x ~= n_gene
    error('decodingfun should return x with %d elements but it returned x with %d elements.',n_gene,length_x);
end
    
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
