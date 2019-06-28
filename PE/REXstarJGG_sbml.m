function [ Results, optimizedmodel ] = REXstarJGG_sbml(model,decodingfun,mst,varargin)
% REXstarJGG_sbml(model,decodingfun,mst);
% REXstarJGG_sbml(model,decodingfun,mst,opts); 4
% REXstarJGG_sbml(model,decodingfun,mst,simopts,opts); 5 
% REXstarJGG_sbml(model,decodingfun,mst,fast_flag,simopts,opts); 6
% REXstarJGG_sbml(model,decodingfun,mst,n_gene,n_constraint,fitnessfun,fast_flag,simopts,opts);

switch nargin
    case 3
        fast_flag = 1;
        fitnessfun = @SSR_sbml;
        n_constraint = 0;
        simopts = struct;
        opts = struct;
    case 4
        fast_flag = 1;
        fitnessfun = @SSR_sbml;
        n_constraint = 0;
        simopts = struct;
        opts = varargin{1};
    case 5
        fast_flag = 1;
        fitnessfun = @SSR_sbml;
        n_constraint = 0;
        simopts = varargin{1};
        opts = varargin{2};
    case 6
        fast_flag = varargin{1};
        fitnessfun = @SSR_sbml;
        n_constraint = 0;
        simopts = varargin{2};
        opts = varargin{3};
    case 9
        n_gene = varargin{1};
        n_constraint = varargin{2};
        fitnessfun = varargin{3};
        fast_flag = varargin{4};
        simopts = varargin{5};
        opts = varargin{6};
    otherwise
        error('Incorrect number of input arguments');
end

% If model is the name of a SBML file
if ~isIQMmodel(model)
    if exist(model,'file') == 2
        fprintf('Reading %s ...',model);
        model = IQMmodel(model);
        fprintf(' Finished.\n');
    else
        error('%s not found.',model);
    end
end
% fast_flag = 0;

if ~exist('n_gene','var') || isempty(n_gene)
    st_model = IQMstruct(model);
    n_gene = length(st_model.parameters);
end
if isempty(n_constraint)
    n_constraint = 0;
end
if isempty(fitnessfun)
    fitnessfun = @SSR_sbml;
end
if isempty(fast_flag)
    fast_flag = 1;
end
if isempty(simopts)
    simopts = struct;
end
if isempty(opts)
    opts = struct;
end

if ~isIQMmeasurement(mst)
    fprintf('Reading %s ...',mst);
    mst = IQMmeasurement(mst);
    fprintf(' Finished.\n');
end

% If model is an SBmodel object and fast_flg is one
if fast_flag == 1
    st_model = IQMstruct(model);
    mex_name = strcat(st_model.name,'_mex');
    clear(mex_name);
    fprintf('Making %s ...\n',mex_name);
    IQMmakeMEXmodel(model,mex_name);
%     fprintf(' Finished.\n');
    Simulation = @Simulation_mex;
    
else
    st_model = IQMstruct(model);
    mex_name = strcat(st_model.name,'_odefun');
    fprintf('Making %s ...\n',mex_name);
    IQMcreateODEfile(model,mex_name);
    mex_name = str2func(mex_name);
    Simulation = @Simulation_odexx;
end

if isempty(which('IQMmodel'))
    error('IQMTools is not properly installed.');
end



problem.n_gene = n_gene;
problem.n_constraint = n_constraint;
% problem.fitnessfun = @(x) fitnessfun(x,mex_name,mst,simopts);
problem.decodingfun = decodingfun;
problem.fitnessfun = @(x) fitnessfun(Simulation,x,mex_name,mst,simopts);
if ~isfield(opts,'interimreportfun')
    opts.interimreportfun = @interimreportfun_sbml;
end
interimreportfun = opts.interimreportfun;
opts.interimreportfun = @(elapsedTime,generation,problem,opts,Population,best) ...
    interimreportfun(...
    Simulation, elapsedTime,generation,problem,opts,Population,best,...
    mex_name,mst,simopts,fast_flag);

Results = REXstarJGG(problem,opts);

if isIQMmodel(model)
    st_model = struct(model);
    for i = 1 : length(st_model.parameters)
        st_model.parameters(i).value = Results.Best.x(i);
    end
    optimizedmodel = IQMmodel(st_model);
else
    optimizedmodel = [];
end
