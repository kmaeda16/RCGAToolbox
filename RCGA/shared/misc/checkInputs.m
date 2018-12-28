function [problem, opts] = checkInputs(problem,opts,RCGA_Type)

if isempty(strfind(RCGA_Type,'REXstarJGG')) && isempty(strfind(RCGA_Type,'UNDXMGG'))
    error('RCGA_Type must include REXstarJGG or UNDXMGG');
end

%%

C1 = {
    'n_gene',...         %  1
    'n_constraint',...   %  2
    'fitnessfun',...     %  3
    'decodingfun',...    %  4
    };

tf = isfield(problem,C1);

if ~tf(1) % n_gene
    error('n_gene needs to be set!');
end
if ~tf(2) % n_constraint
    disp('n_constraint not provided. Default value used (i.e. n_constraint = 0).');
    problem.n_constraint = 0;
end
if ~tf(3) % fitnessfun
    error('fitnessfun needs to be set!');
end
if ~tf(4) % decodingfun
    error('decodingfun needs to be set!');
end

%%

C2 = {
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
    'out_best',...       % 12
    'out_population',... % 13
    'out_report',...     % 14
    'interimreportfun',... % 15
    'finalreportfun',... % 16
    'par', ...           % 17
    };

tf = isfield(opts,C2);

if ~isempty(strfind(RCGA_Type,'UNDXMGG'))
    tf(3:5) = 1;
end

if ~tf(1) % n_population
%     disp('n_population not provided. Default value used (i.e. n_population = 20 * n_gene).');
    opts.n_population = 20 * problem.n_gene;
end
if ~tf(2) % n_children
%     disp('n_children not provided. Default value used (i.e. n_children = 5 * n_gene).');
    opts.n_children = 5 * problem.n_gene;
end
if ~tf(3) % n_parent
%     disp('n_parent not provided. Default value used (i.e. n_parent = n_gene + 1).');
    opts.n_parent = problem.n_gene + 1;
end
if ~tf(4) % t_rextar
%     disp('t_rextar not provided. Default value used (i.e. t_rextar = 6.0).');
    opts.t_rextar = 6.0;
end
if ~tf(5) % selection_type
%     disp('selection_type not provided. Default value used (i.e. selection_type = 0).');
    opts.selection_type = 0;
end
if ~tf(6) % Pf
    if problem.n_constraint > 0
%         disp('Pf not provided. Default value used (i.e. Pf = 0.45).');
        opts.Pf = 0.45;
    else
%         disp('Pf not provided. Default value used (i.e. Pf = 0).');
        opts.Pf = 0;
    end
end
if ~tf(7) % n_generation
%     disp('n_generation not provided. Default value used (i.e. n_generation = 1000).');
    opts.n_generation = 1000;
end
if ~tf(8) % t_limit
%     disp('t_limit not provided. Default value used (i.e. t_limit = 1000).');
    opts.t_limit = 1000;
end
if ~tf(9) % vtr
%     disp('vtr not provided. Default value used (i.e. vtr = 0).');
    opts.vtr = 0;
end
if ~tf(10) % output_intvl
%     disp('output_intvl not provided. Default value used (i.e. output_intvl = 1).');
    opts.output_intvl = 1;
end
if ~tf(11) % out_transition
%     disp('out_transition not provided. Default value used (i.e. out_transition = None).');
    opts.out_transition = 'None';
end
if ~tf(12) % out_best
%     disp('out_best not provided. Default value used (i.e. out_best = None).');
    opts.out_best = 'None';
end
if ~tf(13) % out_population
%     disp('out_population not provided. Default value used (i.e. out_population = None).');
    opts.out_population = 'None';
end
if ~tf(14) % out_report
%     disp('out_report not provided. Default value used (i.e. out_report = None).');
    opts.out_report = 'None';
end
if ~tf(15) % interimreportfun
    %     disp('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.interimreportfun = @(elapsedTime,i,Param,best) x;
    opts.interimreportfun = @defaultinterimreportfun;
end
if ~tf(16) % finalreportfun
    %     disp('out_population not provided. Default value used (i.e. out_population = None).');
%     Param.finalreportfun = @(elapsedTime,i,Param,best) x;
    opts.finalreportfun = @defaultfinalreportfun;
end
if ~tf(17) % par
%     disp('par not provided. Default value used (i.e. par = 0).');
    opts.par = 0;
end
