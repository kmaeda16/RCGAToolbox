function RCGAwritePopulation(problem, opts, Population)
% RCGAwritePopulation makes a text file with population
% 
% [SYNTAX]
% RCGAwritePopulation(problem, opts, Population)
% 
% [INPUT]
% problem    :  Problem structure
% opts       :  RCGA options. See XXXXXXXXXXX for options.
% Population :  Array of chrom


%% Shortening variable names
decodingfun = problem.decodingfun;
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
n_population = opts.n_population;
out_population = opts.out_population;


%% If out_transition is 'none', nothing done.
if strcmpi('none',out_population)
    return;
end


%% Opening out_transition
out = fopen(out_population,'w');
if out == -1
    warning('cannot open %s!\n',out_transition);
    return;
end

fprintf(out,'f\t');
if n_constraint > 0
    fprintf(out,'phi\t');
end
for j = 1 : n_gene
    fprintf(out,'x(%d)\t',j);
end
for j = 1 : n_constraint
    fprintf(out,'g(%d)\t',j);
end
fprintf(out,'\n');


%% Making outputs
for i = 1 : n_population
    x = decodingfun(Population(i).gene);
    fprintf(out,'%e\t',Population(i).f);
    if n_constraint > 0
        fprintf(out,'%e\t',Population(i).phi);
    end
    for j = 1 : n_gene
        fprintf(out,'%e\t',x(j));
    end
    for j = 1 : n_constraint
        fprintf(out,'%e\t',Population(i).g(j));
    end
    fprintf(out,'\n');
end

fclose(out);
