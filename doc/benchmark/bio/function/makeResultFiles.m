function makeResultFiles(infilename,outfilename,fitnessfun)


load(infilename);

[n_row, n_gene] = size(Results.x);

out = fopen(outfilename,'w');
if out == -1
        warning('cannot open %s!\n',outfilename);
        return;
end

fprintf(out,'Time\tNEval\tGeneration\tf\tphi\t');
for i = 1 : n_gene
    fprintf(out,'x(%d)\t',i);
end
[~, g] = fitnessfun(Results.xbest);
n_constraint = length(g);
for i = 1 : n_constraint
    fprintf(out,'g(%d)\t',i);
end
fprintf(out,'\n');

for i = 1 : n_row
    
    elapsedTime = Results.time(i);
    neval = Results.neval(i);
    generation = nan;
    f = Results.f(i);
    x = Results.x(i,:);
    [~, g] = fitnessfun(x);
    phi = sum( max(0,g) .^2 );
    
    fprintf(out,'%e\t%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,f,phi);
    for j = 1 : n_gene
        fprintf(out,'%e\t',x(j));
    end
    for j = 1 : n_constraint
        fprintf(out,'%e\t',g(j));
    end
    fprintf(out,'\n');
end

fclose(out);
