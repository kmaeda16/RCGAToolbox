function Results = REXstarJGG_PE(model,decodingfun,mst,varargin)
% function [ Results, optimizedmodel ] = REXstarJGG_PE(model,decodingfun,mst,varargin)
% REXstarJGG_PE(model,decodingfun,mst);
% REXstarJGG_PE(model,decodingfun,mst,opts); 4
% REXstarJGG_PE(model,decodingfun,mst,simopts,opts); 5 
% REXstarJGG_PE(model,decodingfun,mst,fast_flag,simopts,opts); 6
% REXstarJGG_PE(model,decodingfun,mst,n_constraint,fitnessfun,fast_flag,simopts,opts); 8

%% Handling input arguments
switch nargin
    case 3
        n_constraint = 0;
        fitnessfun = @SSR_sbml;
        fast_flag = 0;
        simopts = struct;
        opts = struct;
    case 4
        n_constraint = 0;
        fitnessfun = @SSR_sbml;
        fast_flag = 0;
        simopts = struct;
        opts = varargin{1};
    case 5
        n_constraint = 0;
        fitnessfun = @SSR_sbml;
        fast_flag = 0;
        simopts = varargin{1};
        opts = varargin{2};
    case 6
        n_constraint = 0;
        fitnessfun = @SSR_sbml;
        fast_flag = varargin{1};
        simopts = varargin{2};
        opts = varargin{3};
    case 8
        n_constraint = varargin{1};
        fitnessfun = varargin{2};
        fast_flag = varargin{3};
        simopts = varargin{4};
        opts = varargin{5};
    otherwise
        error('Incorrect number of input arguments');
end

if isempty(n_constraint)
    n_constraint = 0;
end
if isempty(fitnessfun)
    fitnessfun = @SSR_sbml;
end
if isempty(fast_flag)
    fast_flag = 0;
end
if isempty(simopts)
    simopts = struct;
end
if isempty(opts)
    opts = struct;
end


%% If model is an IQMmodel, it is converted into odefun.m or odefun.c.
if isIQMmodel(model)
    % st_model = IQMstruct(model);
    st_model = struct(model);
    switch fast_flag
        case {0, 1}
            odefun_name = strcat(st_model.name,'_odefun');
            fprintf('Making %s ...\n',odefun_name);
            IQMcreateODEfile(model,odefun_name);
            model = str2func(odefun_name);
        case 2
            mex_name = strcat(st_model.name,'_mex');
            clear(mex_name);
            fprintf('Making %s ...\n',mex_name);
            IQMmakeMEXmodel(model,mex_name,1);
            mexcompileIQM(mex_name);
            model = mex_name;
        otherwise
            error('Unexpected fast_flag!');
    end
end

%% if model is a string, it is converted into a function handle
if ischar(model)
    
    [~, filename, ext ] = fileparts(model);
    
    if strcmp('.m',ext) || strcmp(strcat('.',mexext),ext)
        model = str2func(model);
    elseif strcmp('.c',ext)
        mexcompileIQM(filename);
        model = str2func(filename);
    elseif strcmp('.sbml',ext) || strcmp('.xml',ext)
        sbm = IQMmodel(model);
        st_model = struct(sbm);
        switch fast_flag
            case {0, 1}
                sbm = IQMmodel(model);
                odefun_name = strcat(st_model.name,'_odefun');
                fprintf('Making %s ...\n',odefun_name);
                IQMcreateODEfile(sbm,odefun_name);
                model = str2func(odefun_name);
            case 2
                mex_name = strcat(st_model.name,'_mex');
                clear(mex_name);
                fprintf('Making %s ...\n',mex_name);
                IQMmakeMEXmodel(model,mex_name,1);
                mexcompileIQM(mex_name);
                model = mex_name;
            otherwise
                error('Unexpected fast_flag!');
        end 
    elseif isempty(ext)
        if exist(strcat(model,'.m'),'file') || exist(strcat(model,'.',mexext),'file')
            model = str2func(model);
        elseif exist(strcat(model,'.c'),'file')
            mexcompileIQM(model);
            model = str2func(model);
        elseif exist(strcat(model,'.sbml'),'file')
            sbm = IQMmodel(strcat(model,'.sbml'));
            st_model = struct(sbm);
            switch fast_flag
                case {0, 1}
                    odefun_name = strcat(st_model.name,'_odefun');
                    fprintf('Making %s ...\n',odefun_name);
                    IQMcreateODEfile(model,odefun_name);
                    model = str2func(odefun_name);
                case 2
                    mex_name = strcat(st_model.name,'_mex');
                    clear(mex_name);
                    fprintf('Making %s ...\n',mex_name);
                    IQMmakeMEXmodel(model,mex_name,1);
                    mexcompileIQM(mex_name);
                    model = mex_name;
                otherwise
                    error('Unexpected fast_flag!');
            end
        elseif exist(strcat(model,'.xml'),'file')
            sbm = IQMmodel(strcat(model,'.sbml'));
            st_model = struct(sbm);
            switch fast_flag
                case {0, 1}
                    odefun_name = strcat(st_model.name,'_odefun');
                    fprintf('Making %s ...\n',odefun_name);
                    IQMcreateODEfile(model,odefun_name);
                    model = str2func(odefun_name);
                case 2
                    mex_name = strcat(st_model.name,'_mex');
                    clear(mex_name);
                    fprintf('Making %s ...\n',mex_name);
                    IQMmakeMEXmodel(model,mex_name,1);
                    mexcompileIQM(mex_name);
                    model = mex_name;
                otherwise
                    error('Unexpected fast_flag!');
            end
        else
            error('File %s does not exist!',model);
        end
    end

end

%% Error handling

filename = func2str(model);
file_type = exist(filename,'file');
if  file_type == 0
    error('%s does not exist!',model);
elseif file_type < 2 || 3 < file_type
    error('%s is neither an m file nor a MEX file!',model);
end

try
    output = feval(model,'parametervalues');
    n_param = length(output);
catch ME
    warning(ME.identifier);
    error('%s should return default parameter values when ''parametervalues'' is given as an argument!');
end

%% Set

% odefun.m
switch fast_flag
    case 0
        Simulation = @Simulation_odexx;
    case 1
        Simulation = @Simulation_stb;
    case 2
        Simulation = @Simulation_mex;
    otherwise
        error('Unexpected fast_flag!');
end

filename = func2str(model);
exist_flag = exist(filename,'file');

if exist_flag == 2 && ~( fast_flag == 0 || fast_flag == 1 )
    warning('fast_flg set to 0!');
    fast_flag = 0;
    Simulation = @Simulation_odexx;
end

if exist_flag == 3 && ~( fast_flag == 2 )
    warning('fast_flg set to 2!');
    fast_flag = 2;
    Simulation = @Simulation_mex;
end

%% Reading measurement

if ~isIQMmeasurement(mst)
    fprintf('Reading %s ...',mst);
    mst = IQMmeasurement(mst);
    fprintf(' Finished.\n');
end

%% Prepering inputs for REXstarJGG

problem.n_gene = n_param;
problem.n_constraint = n_constraint;
% problem.fitnessfun = @(x) fitnessfun(x,mex_name,mst,simopts);
problem.decodingfun = decodingfun;
problem.fitnessfun = @(x) fitnessfun(Simulation,x,model,mst,simopts);
if ~isfield(opts,'interimreportfun')
    opts.interimreportfun = @interimreportfun_sbml;
end
interimreportfun = opts.interimreportfun;
opts.interimreportfun = @(elapsedTime,generation,problem,opts,Population,best) ...
    interimreportfun(...
    Simulation, elapsedTime,generation,problem,opts,Population,best,...
    model,mst,simopts,fast_flag);

%% Run parameter estimation

Results = REXstarJGG(problem,opts);

%% If model was a IQMmodel

% if isIQMmodel(model)
%     st_model = struct(model);
%     for i = 1 : length(st_model.parameters)
%         st_model.parameters(i).value = Results.Best.x(i);
%     end
%     optimizedmodel = IQMmodel(st_model);
% else
%     optimizedmodel = [];
% end
