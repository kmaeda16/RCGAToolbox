function problem = checkProblem(problem)

C = {
    'n_gene',...         %  1
    'n_constraint',...   %  2
    'fitnessfun',...     %  3
    'decodingfun',...    %  4
    };

tf = isfield(problem,C);

if ~tf(1) % n_gene
    error('n_gene needs to be set!');
end
if ~tf(2) % n_constraint
    warning('n_constraint not provided. Default value used (i.e. n_constraint = 0).');
    problem.n_constraint = 0;
end
if ~tf(3) % fitnessfun
    error('fitnessfun needs to be set!');
end
if ~tf(4) % decodingfun
    warning('decodingfun needs to be set!');
end
