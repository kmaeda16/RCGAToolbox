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
    't_rexstar',...       %  4
    'selection_type',... %  5
    'Pf',...             %  6
    'local',...          %  7
    'localopts',...      %  8
    'n_generation',...   %  9
    'maxtime',...        % 10
    'maxeval',...        % 11
    'vtr',...            % 12
    'n_par', ...         % 13
    'output_intvl',...   % 14
    'out_transition',... % 15
    'out_best',...       % 16
    'out_population',... % 17
    'out_report',...     % 18
    'interimreportfun',... % 19
    'finalreportfun',... % 20
};

tf = isfield(opts,C2);

if strcmpi(RCGA_Type,'UNDXMGG')
    tf(3:5) = 1;
end

if ~tf(1) % n_population
    opts.n_population = 20 * problem.n_gene; % Recommended by Kobayashi 2009
%     opts.n_population = 5 * problem.n_gene;
end
if ~tf(2) % n_children
    opts.n_children = 5 * problem.n_gene; % Recommended by Kobayashi 2009
%     opts.n_children = ceil(0.5 * opts.n_population);
%     opts.n_children = opts.n_population;
end
if ~tf(3) % n_parent
    opts.n_parent = problem.n_gene + 1; % Recommended by Kobayashi 2009
%     opts.n_parent = ceil(0.5 * opts.n_population);
end
if ~tf(4) % t_rexstar
    opts.t_rexstar = 6.0;
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
if ~tf(7) % local
    opts.local = 0;
end
if ~tf(8) % localopts
    opts.localopts = optimoptions('fmincon',...
        'ConstraintTolerance',0,...
        'MaxFunctionEvaluations',opts.n_children,...
        'Display','off');
end
if ~tf(9) % n_generation
    opts.n_generation = 1000;
end
if ~tf(10) % maxtime
    opts.maxtime = 10 * 60; % 10 min
end
if ~tf(11) % maxeval
    opts.maxeval = inf;
end
if ~tf(12) % vtr
    opts.vtr = -inf;
end
if ~tf(13) % n_par
    opts.n_par = 1;
end
if ~tf(14) % output_intvl
    opts.output_intvl = 1;
end
if ~tf(15) % out_transition
    opts.out_transition = 'None';
end
if ~tf(16) % out_best
    opts.out_best = 'None';
end
if ~tf(17) % out_population
    opts.out_population = 'None';
end
if ~tf(18) % out_report
    opts.out_report = 'None';
end
if ~tf(19) % interimreportfun
    opts.interimreportfun = @RCGAdefaultinterimreportfun;
end
if ~tf(20) % finalreportfun
    opts.finalreportfun = @RCGAdefaultfinalreportfun;
end


