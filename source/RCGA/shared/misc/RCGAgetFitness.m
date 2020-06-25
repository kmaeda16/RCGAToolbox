function [f, g, phi] = RCGAgetFitness(problem,chrom)
% RCGAgetFitness returns f, g, and phi for an individual chrom.
% 
% [SYNTAX]
% [f, g, phi] = RCGAgetFitness(problem,chrom)
% 
% [INPUT]
% problem :  Problem structure
% chrom   :  Individual
% 
% [OUTPUT]
% f       :  Value of fitness functions
% g       :  Values of constraint functions
% phi     :  Value of penalty function


%% Shortening variable names
n_gene = problem.n_gene;
n_constraint = problem.n_constraint;
fitnessfun = problem.fitnessfun;
decodingfun = problem.decodingfun;


%% Decoding genes to x
x = feval(decodingfun,chrom.gene);
length_x = length(x);

if length_x ~= n_gene
    error('decodingfun should return a vector with %d elements but it returned a vector with %d elements.',n_gene,length_x);
end


%% Getting f, g, and phi
if n_constraint == 0
    if nargout(fitnessfun) > 1
        warning('n_constraint was set to %d, but %s returns g.',n_constraint,func2str(fitnessfun));
    end
    f = feval(fitnessfun,x);
    if ~isreal(f)
        warning('f is a complex value. The imaginary part was discarded.');
        f = real(f);
    end
    g = 0;
    phi = 0;
else
    if nargout(fitnessfun) == 1
        error('n_constraint was set to %d, but %s does not returns g.',n_constraint,func2str(fitnessfun));
    end
    [f, g] = feval(fitnessfun,x);
    if ~isreal(f) || ~isreal(g)
        warning('f and/or g are complex values. The imaginary parts were discarded.');
        f = real(f);
        g = real(g);
    end
    phi = sum( max(0,g) .^ 2 );
    length_g = length(g);
    if n_constraint ~= length_g
        error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
    end
end
