function [problem, opts] = RCGAcheckInputs(problem,opts,RCGA_Type)
% RCGAcheckInputs checks problem and opts and sets default values
% 
% [SYNTAX]
% [problem, opts] = RCGAcheckInputs(problem,opts,RCGA_Type)
% 
% [INPUT]
% problem   :  Problem structure.
% opts      :  RCGA options. See XXXXXXXXXXX for options.
% RCGA_Type :  
% 
% [OUTPUT]
% problem   :  Problem structure with default values.
% opts      :  RCGA options with default values.


%% Checking RCGA_Type
if ~strcmp(RCGA_Type,'RCGA_REXstarJGG') && ~strcmp(RCGA_Type,'RCGA_UNDXMGG')
    error('RCGA_Type must be ''RCGA_REXstarJGG'' or ''RCGA_UNDXMGG''.');
end


%% Checking problem
C1 = {
    'fitnessfun',...     % 1
    'decodingfun',...    % 2
    'n_gene',...         % 3
    'n_constraint',...   % 4
    };

tf = isfield(problem,C1);

if ~tf(1) % fitnessfun
    error('fitnessfun needs to be set!');
end
if ~tf(2) % decodingfun
    error('decodingfun needs to be set!');
end
if ~tf(3) % n_gene
    error('n_gene needs to be set!');
end
if ~tf(4) % n_constraint
    disp('n_constraint not provided. Default value used (i.e. n_constraint = 0).');
    problem.n_constraint = 0;
end



%% Checking opts and seting default values
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
    'n_localoptimind',... % 18
    'localoptimopts',... % 19
    };

tf = isfield(opts,C2);

if strcmpi(RCGA_Type,'UNDXMGG')
    tf(3:5) = 1;
end

if ~tf(1) % n_population
    opts.n_population = 20 * problem.n_gene;
end
if ~tf(2) % n_children
    opts.n_children = 5 * problem.n_gene;
end
if ~tf(3) % n_parent
    opts.n_parent = problem.n_gene + 1;
end
if ~tf(4) % t_rextar
    opts.t_rextar = 6.0;
end
if ~tf(5) % selection_type
    opts.selection_type = 0;
end
if ~tf(6) % Pf
    if problem.n_constraint > 0
        opts.Pf = 0.45;
    else
        opts.Pf = 0;
    end
end
if ~tf(7) % n_generation
    opts.n_generation = 1000;
end
if ~tf(8) % t_limit
    opts.t_limit = 1000;
end
if ~tf(9) % vtr
    opts.vtr = -inf;
end
if ~tf(10) % output_intvl
    opts.output_intvl = 1;
end
if ~tf(11) % out_transition
    opts.out_transition = 'None';
end
if ~tf(12) % out_best
    opts.out_best = 'None';
end
if ~tf(13) % out_population
    opts.out_population = 'None';
end
if ~tf(14) % out_report
    opts.out_report = 'None';
end
if ~tf(15) % interimreportfun
    opts.interimreportfun = @RCGAdefaultinterimreportfun;
end
if ~tf(16) % finalreportfun
    opts.finalreportfun = @RCGAdefaultfinalreportfun;
end
if ~tf(17) % par
    opts.par = 0;
end
if ~tf(18) % n_localoptimind
    opts.n_localoptimind = 0;
end
if ~tf(19) % localoptimopts
    opts.localoptimopts = optimoptions('fmincon',...
        'ConstraintTolerance',0,...
        'MaxFunctionEvaluations',opts.n_children,...
        'Display','off');
end

