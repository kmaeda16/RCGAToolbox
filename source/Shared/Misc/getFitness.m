function [f, g, phi] = getFitness(Param,chrom)

n_gene = Param.n_gene;
n_constraint = Param.n_constraint;
fitnessfun = Param.fitnessfun;
decodingfun = Param.decodingfun;
x = feval(decodingfun,chrom.gene);

length_x = length(x);
if length_x ~= n_gene
    error('decodingfun should return x with %d elements but it returned x with %d elements.',n_gene,length_x);
end

% n_output = nargout(fitnessfun);
% switch n_output
%     case 1
%         f = feval(fitnessfun,x);
%         g = 0;
%         phi = 0;
%         if n_constraint ~= 0
%             error('n_constraint was set to %d, but %s does not return g.',n_constraint,func2str(fitnessfun));
%         end
%     case 2
%         [f, g] = feval(fitnessfun,x);
%         phi = sum( max(0,g) .^ 2 );
%         length_g = length(g);
%         if n_constraint ~= length_g
%             error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
%         end
%     otherwise
%         error('fitnessfun should have one or two output arguments.');
% end


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
