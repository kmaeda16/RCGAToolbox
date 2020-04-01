function RCGAwriteTransition_PE(elapsedTime, generation, problem, opts, best, modelfun)
% RCGAwriteTransition is executed at the end of RCGA, making two output
% files.
% 
% [SYNTAX]
% RCGAwriteTransition_PE(elapsedTime, generation, problem, opts, chrom)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec)
% generation  :  Generation
% problem     :  Problem structure
% opts        :  RCGA options. See XXXXXXXXXXX for options.
% best        :  Structure of the the best individual
% modelfun    :  Function handle for model (odefun or mex)


%% Shortening variable names
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
decodingfun = problem.decodingfun;
out_transition = opts.out_transition;
n_population = opts.n_population;
n_children = opts.n_children;


%% If out_transition is 'none', nothing done.
if  strcmpi('none',out_transition)
    return;
end


%% Getting parameter names
paramnames = modelfun('parameters');


%% Opening out_transition
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
        fprintf(out,'%s\t',char(paramnames(i)));
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


%% Decoding genes to x
x = decodingfun(best.gene);
length_x = length(x);

if length_x ~= n_gene
    error('decodingfun should return a vector with %d elements but it returned a vector with %d elements.',n_gene,length_x);
end


%% Calculating the number of fitness function evaluations
global RCGA_LOCALNEVAL;
neval = n_population + ( generation - 1 ) * n_children + RCGA_LOCALNEVAL;


%% Making outputs
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
