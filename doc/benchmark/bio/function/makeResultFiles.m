function makeResultFiles(infilename,outfilename)


load(infilename);

[n_row, n_gene] = size(Results.x);

out = fopen(outfilename,'w');
if out == -1
        warning('cannot open %s!\n',outfilename);
        return;
end

fprintf(out,'Time\tNEval\tGeneration\tf\t');
for i = 1 : n_gene
    fprintf(out,'x(%d)\t',i);
end
fprintf(out,'\n');

for i = 1 : n_row
    elapsedTime = Results.time(i);
    neval = Results.neval(i);
    generation = nan;
    f = Results.f(i);
    x = Results.x(i,:);
    fprintf(out,'%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,f);
    for j = 1 : n_gene
        fprintf(out,'%e\t',x(j));
    end
    fprintf(out,'\n');
end

fclose(out);
