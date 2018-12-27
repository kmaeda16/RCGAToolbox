function opts = checkOpts(opts,RCGA_Type)

if isempty(strfind(RCGA_Type,'REXstarJGG')) && isempty(strfind(RCGA_Type,'UNDXMGG'))
    error('RCGA_Type must include REXstarJGG or UNDXMGG');
end

C = {
    'n_population',...   %  1
    'n_children',...     %  2
    'n_parent',...       %  3
    't_rextar',...       %  4
    'selection_type',... %  5
    'Pf',...             %  6
    'n_generation',...   %  7
    't_limit',...        %  8
    'vtr',...            %  9
    'output_intvl',...   % 10
    'out_transition',... % 11
    'out_solution',...   % 12
    'out_population',... % 13
    'interimreportfun',... % 14
    'finalreportfun',... % 15
    };

tf = isfield(opts,C);

if ~isempty(strfind(RCGA_Type,'UNDXMGG'))
    tf(3:5) = 1;
end

if ~tf(1) % n_population
    warning('n_population not provided. Default value used (i.e. n_population = 20 * n_gene).');
    opts.n_population = 20 * opts.n_gene;
end
if ~tf(2) % n_children
    warning('n_children not provided. Default value used (i.e. n_children = 5 * n_gene).');
    opts.n_children = 5 * opts.n_gene;
end
if ~tf(3) % n_parent
    warning('n_parent not provided. Default value used (i.e. n_parent = n_gene + 1).');
    opts.n_parent = opts.n_gene + 1;
end
if ~tf(4) % t_rextar
    warning('t_rextar not provided. Default value used (i.e. t_rextar = 6.0).');
    opts.t_rextar = 6.0;
end
if ~tf(5) % selection_type
    warning('selection_type not provided. Default value used (i.e. selection_type = 0).');
    opts.selection_type = 0;
end
if ~tf(6) % Pf
    if opts.n_constraint > 0
        warning('Pf not provided. Default value used (i.e. Pf = 0.45).');
        opts.Pf = 0.45;
    else
        warning('Pf not provided. Default value used (i.e. Pf = 0).');
        opts.Pf = 0;
    end
end
if ~tf(7) % n_generation
    warning('n_generation not provided. Default value used (i.e. n_generation = 1000).');
    opts.n_generation = 1000;
end
if ~tf(8) % t_limit
    warning('t_limit not provided. Default value used (i.e. t_limit = 1000).');
    opts.t_limit = 1000;
end
if ~tf(9) % vtr
    warning('vtr not provided. Default value used (i.e. vtr = 0).');
    opts.vtr = 0;
end
if ~tf(10) % output_intvl
    warning('output_intvl not provided. Default value used (i.e. output_intvl = 1).');
    opts.output_intvl = 1;
end
if ~tf(11) % out_transition
    warning('out_transition not provided. Default value used (i.e. out_transition = None).');
    opts.out_transition = 'None';
end
if ~tf(12) % out_solution
    warning('out_solution not provided. Default value used (i.e. out_solution = None).');
    opts.out_solution = 'None';
end
if ~tf(13) % out_population
    warning('out_population not provided. Default value used (i.e. out_population = None).');
    opts.out_population = 'None';
end
if ~tf(14) % interimreportfun
    %     warning('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.interimreportfun = @(elapsedTime,i,Param,best) x;
    opts.interimreportfun = @defaultinterimreportfun;
end
if ~tf(15) % finalreportfun
    %     warning('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.finalreportfun = @(elapsedTime,i,Param,best) x;
    opts.finalreportfun = @defaultfinalreportfun;
end

