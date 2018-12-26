function Param = checkInputs(Param,RCGA)

C = {
    'n_gene',...         %  1
    'n_constraint',...   %  2
    'fitnessfun',...     %  3
    'decodingfun',...    %  4
    'n_population',...   %  5
    'n_children',...     %  6
    'n_parent',...       %  7
    't_rextar',...       %  8
    'selection_type',... %  9
    'Pf',...             % 10
    'n_generation',...   % 11
    't_limit',...        % 12
    'vtr',...            % 13
    'output_intvl',...   % 14
    'out_transition',... % 15
    'out_solution',...   % 16
    'out_population',... % 17
    'interimreportfun',... % 18
    'finalreportfun',... % 19
    'opts',...           % 20
    };

tf = isfield(Param,C);

if strcmp(RCGA,'UNDXMGG')
    tf(7:9) = 1;
end

if ~tf(1) % n_gene
    error('n_gene needs to be set!');
end
if ~tf(2) % n_constraint
    warning('n_constraint not provided. Default value used (i.e. n_constraint = 0).');
    Param.n_constraint = 0;
end
if ~tf(3) % fitnessfun
    error('fitnessfun needs to be set!');
end
if ~tf(4) % decodingfun
    warning('decodingfun is not provided. Default decodingfun used (i.e. x = gene .* ( ub - lb ) + lb with lb = - 1e+5 and ub = - 1e+5).');
    Param.decodingfun = @defaultdecodingfun;
end
if ~tf(5) % n_population
    warning('n_population not provided. Default value used (i.e. n_population = 20 * n_gene).');
    Param.n_population = 20 * Param.n_gene;
end
if ~tf(6) % n_children
    warning('n_children not provided. Default value used (i.e. n_children = 5 * n_gene).');
    Param.n_children = 5 * Param.n_gene;
end
if ~tf(7) % n_parent
    warning('n_parent not provided. Default value used (i.e. n_parent = n_gene + 1).');
    Param.n_parent = Param.n_gene + 1;
end
if ~tf(8) % t_rextar
    warning('t_rextar not provided. Default value used (i.e. t_rextar = 6.0).');
    Param.t_rextar = 6.0;
end
if ~tf(9) % selection_type
    warning('selection_type not provided. Default value used (i.e. selection_type = 0).');
    Param.selection_type = 0;
end
if ~tf(10) % Pf
    if Param.n_constraint > 0
        warning('Pf not provided. Default value used (i.e. Pf = 0.45).');
        Param.Pf = 0.45;
    else
        warning('Pf not provided. Default value used (i.e. Pf = 0).');
        Param.Pf = 0;
    end
end
if ~tf(11) % n_generation
    warning('n_generation not provided. Default value used (i.e. n_generation = 1000).');
    Param.n_generation = 1000;
end
if ~tf(12) % t_limit
    warning('t_limit not provided. Default value used (i.e. t_limit = 1000).');
    Param.t_limit = 1000;
end
if ~tf(13) % vtr
    warning('vtr not provided. Default value used (i.e. vtr = 0).');
    Param.vtr = 0;
end
if ~tf(14) % output_intvl
    warning('output_intvl not provided. Default value used (i.e. output_intvl = 1).');
    Param.output_intvl = 1;
end
if ~tf(15) % out_transition
    warning('out_transition not provided. Default value used (i.e. out_transition = None).');
    Param.out_transition = 'None';
end
if ~tf(16) % out_solution
    warning('out_solution not provided. Default value used (i.e. out_solution = None).');
    Param.out_solution = 'None';
end
if ~tf(17) % out_population
    warning('out_population not provided. Default value used (i.e. out_population = None).');
    Param.out_population = 'None';
end
if ~tf(18) % interimreportfun
    %     warning('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.interimreportfun = @(elapsedTime,i,Param,best) x;
    Param.interimreportfun = @defaultinterimreportfun;
end
if ~tf(19) % finalreportfun
    %     warning('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.finalreportfun = @(elapsedTime,i,Param,best) x;
    Param.finalreportfun = @defaultfinalreportfun;
end
if ~tf(20) % opts
    %     warning('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.finalreportfun = @(elapsedTime,i,Param,best) x;
    Param.opts = struct;
end
