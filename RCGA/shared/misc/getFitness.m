function [f, g, phi] = getFitness(problem,chrom)

n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
fitnessfun = problem.fitnessfun;
decodingfun = problem.decodingfun;
x = feval(decodingfun,chrom.gene);

length_x = length(x);
if length_x ~= n_gene
    error('decodingfun should return x with %d elements but it returned x with %d elements.',n_gene,length_x);
end

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
