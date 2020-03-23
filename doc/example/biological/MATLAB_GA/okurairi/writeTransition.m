function writeTransition(elapsedTime, generation, problem, opts, x, neval)
% writeTransition is executed at the end of RCGA, making two output files.
% 
% [SYNTAX]
% RCGAwriteTransition(elapsedTime, generation, problem, opts, chrom)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec)
% generation  :  Generation
% problem     :  Problem structure
% opts        :  RCGA options. See XXXXXXXXXXX for options.
% best        :  Structure of the the best individual


%% Shortening variable names
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
decodingfun = problem.decodingfun;
out_transition = opts.out_transition;
n_population = opts.n_population;
n_children = opts.n_children;
fitnessfun = problem.fitnessfun;


%% If out_transition is 'none', nothing done.
if  strcmpi('none',out_transition)
    return;
end


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


%% Decoding genes to x
% x = decodingfun(best.gene);
length_x = length(x);

if length_x ~= n_gene
    error('decodingfun should return a vector with %d elements but it returned a vector with %d elements.',n_gene,length_x);
end


%% Calculating the number of fitness function evaluations
% global RCGA_LOCALNEVAL;
% neval = n_population + ( generation - 1 ) * n_children + RCGA_LOCALNEVAL;


%% Getting f, g, and phi
if n_constraint == 0
    if nargout(fitnessfun) > 1
        warning('n_constraint was set to %d, but %s returns g.',n_constraint,func2str(fitnessfun));
    end
    f = feval(fitnessfun,x);
    g = 0;
    phi = 0;
else
    if nargout(fitnessfun) == 1
        error('n_constraint was set to %d, but %s does not returns g.',n_constraint,func2str(fitnessfun));
    end
    [f, g] = feval(fitnessfun,x);
    phi = sum( max(0,g) .^ 2 );
    length_g = length(g);
    if n_constraint ~= length_g
        error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
    end
end

%% Making outputs
fprintf(out,'%e\t%e\t%e\t%e\t',elapsedTime,neval,generation,f);
if n_constraint > 0
    fprintf(out,'%e\t',phi);
end
for i = 1 : n_gene
    fprintf(out,'%e\t',x(i));
end
for i = 1 : n_constraint
    fprintf(out,'%e\t',g(i));
end
fprintf(out,'\n');

fclose(out);
