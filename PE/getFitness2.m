function [f, g, phi] = getFitness2(Param,chrom)

Param.n_constraint = 0;
% Param.decodingfun = @(x) x;
Param.fitnessfun = @mySSR;
% Param.model = SBmodel('../app/SBMLexampleLevel2.xml');
% Param.mst = SBmeasurement('../app/MeasurementExample.xls');
% chrom.gene = ( 0.1 : 0.1 : 0.9 );
n_gene = length(chrom.gene);


n_constraint = Param.n_constraint;
% decodingfun = Param.decodingfun;
fitnessfun = Param.fitnessfun;
model = Param.model;
mst = Param.mst;
ub = Param.ub;
lb = Param.lb;
mex_name = Param.mex_name;

% x = decodingfun(chrom.gene);
x = chrom.gene .* ( ub - lb ) + lb;
temp = struct(model);
for i = 1 : n_gene
    temp.parameters(i).value = x(i);
end
model = SBmodel(temp);

if Param.n_constraint <= 0
%     f = fitnessfun(model,mst);
    f = fitnessfun(model,mex_name,mst);
    g = 0;
    phi = 0;
    if nargout(fitnessfun) >= 2
        error('%s seems to have f and g as outputs, but g is not used because n_constraint was set to %d.',func2str(fitnessfun),n_constraint);
    end
else
%     [f, g] =fitnessfun(model,mst);
    [f, g] =fitnessfun(model,mex_name,mst);
    phi = sum( max(0,g) .^ 2 );
    length_g = length(g);
    if n_constraint ~= length_g
        error('n_constraint was set to %d, but the length of g returned from %s is %d.',n_constraint,func2str(fitnessfun),length_g);
    end
end
