function [f, g, phi] = getFitness(Param,chrom)

n_constraint = Param.n_constraint;
fitnessfun = Param.fitnessfun;
ub = Param.ub;
lb = Param.lb;

x = chrom.gene .* ( ub - lb ) + lb;

n_output = nargout(fitnessfun);
switch n_output
    case 1
        f = feval(fitnessfun,x);
        g = 0;
        phi = 0;
        if n_constraint ~= 0
            error('n_constraint was set to %d, but %s does not return g.',n_constraint,func2str(fitnessfun));
        end
    case 2
        [f, g] = feval(fitnessfun,x);
        phi = sum( max(0,g) .^ 2 );
        length_g = length(g);
        if n_constraint ~= length_g
            error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
        end
    otherwise
        error('fitnessfun should have one or two output arguments.');
end


