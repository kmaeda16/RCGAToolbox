function makeResultFiles_MATLAB_GA(TSVfilename,finalfilename,outfilename)


TSV = load(TSVfilename);
final = table2array(readtable(finalfilename));

[n_row, ~] = size(TSV);
n_gene = length(final) - 5; % Final element is NaN

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
    
    elapsedTime = final(1) * TSV(i,2) / TSV(end,2);
    neval = TSV(i,2);
    generation = TSV(i,1);
    f = TSV(i,3);
    x = nan(1,n_gene);
    if i == n_row
        x = final(5:end-1);
    end
    
    fprintf(out,'%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,f);
    for j = 1 : n_gene
        fprintf(out,'%e\t',x(j));
    end
    fprintf(out,'\n');
    
end

fclose(out);
