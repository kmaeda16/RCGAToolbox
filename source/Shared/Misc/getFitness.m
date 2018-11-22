function [f, g, phi] = getFitness(Param,chrom)

n_constraint = Param.n_constraint;
decodingfun = Param.decodingfun;
fitnessfun = Param.fitnessfun;

x = decodingfun(chrom.gene);

if Param.n_constraint <= 0
    f = fitnessfun(x);
    g = 0;
    phi = 0;
    if nargout(fitnessfun) >= 2
        error('%s seems to have f and g as outputs, but g is not used because n_constraint was set to %d.',func2str(fitnessfun),n_constraint);
    end
else
    [f, g] = fitnessfun(x);
    phi = sum( max(0,g) .^ 2 );
    length_g = length(g);
    if n_constraint ~= length_g
        error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
    end
end
