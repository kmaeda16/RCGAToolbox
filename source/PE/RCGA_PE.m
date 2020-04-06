function Results = RCGA_PE(model,decodingfun,mst,n_constraint,fitnessfun,fast_flag,simopts,opts,GenerationAlternation)
% RCGA_PE estimates parameters included in model by fitting it to
% experimental data mst.
% 
% [SYNTAX]
% Results = RCGA_PE(model,decodingfun,mst,n_constraint, ...
%                         fitnessfun,fast_flag,simopts,opts)
% 
% [INPUT]
% model        :  IQMmodel or file name (*.sbml, *.xml, *.m, *.mex, or *.c)
% decodingfun  :  Function handle for decoding function
% mst          :  Experimental data (IQMmeasurement) or filename
% n_constraint :  Number of constraints.
% fitnessfun   :  Function handle for fitness function
% fast_flag    :  Solver flag
%                 * fast_flag = 0 : ODEXX by MATLAB
%                 * fast_flag = 1 : CVODE by SundialsTB
%                 * fast_flag = 2 : CVODE by IQM Tools
% simopts      :  Structure with integrator options. Fields depend on
%                 Simulation_*. See 'help Simulation_*'.
% opts         :  Structure with RCGA options
% 
% [OUTPUT]
% Results      :  Objective function value (scaler)


%% If model is an IQMmodel, it is converted into odefun.m or odefun.c.
%  After this section, model will be a string '*_odefun' or '*_mex'.
if isIQMmodel(model)
    st_model = struct(model);
    switch fast_flag
        case {0, 1}
            odefun_name = strcat(st_model.name,'_odefun');
            odefun_name = regexprep(odefun_name,'\W','');
            fprintf('Making %s.m ...',odefun_name);
            IQMcreateODEfile(model,odefun_name);
            fprintf(' Finished.\n');
            model = odefun_name;
        case 2
            mex_name = strcat(st_model.name,'_mex');
            mex_name = regexprep(mex_name,'\W','');
            clear(mex_name);
            fprintf('Making %s.%s ...',mex_name,mexext);
            IQMmakeMEXmodel(model,mex_name,1);
            mexcompileIQM(mex_name);
            fprintf(' Finished.\n');
            model = mex_name;
        otherwise
            error('Unexpected fast_flag!');
    end
end


%% if model is a string, it is converted into a function handle
if ischar(model)
    
    [ ~, filename, ext ] = fileparts(model);
    
    if strcmp('.m',ext) || strcmp(strcat('.',mexext),ext)
        model = str2func(filename);
    elseif strcmp('.c',ext)
        mexcompileIQM(filename);
        model = str2func(filename);
    elseif strcmp('.sbml',ext) || strcmp('.xml',ext)
        sbm = IQMmodel(model);
        switch fast_flag
            case {0, 1}
                odefun_name = strcat(filename,'_odefun');
                fprintf('Making %s.m ...',odefun_name);
                IQMcreateODEfile(sbm,odefun_name);
                fprintf(' Finished.\n');
                model = str2func(odefun_name);
            case 2
                mex_name = strcat(filename,'_mex');
                clear(mex_name);
                fprintf('Making %s.%s ...',mex_name,mexext);
                IQMmakeMEXmodel(sbm,mex_name,1);
                mexcompileIQM(mex_name);
                fprintf(' Finished.\n');
                model = str2func(mex_name);
            otherwise
                error('Unexpected fast_flag!');
        end
    elseif isempty(ext)
        if exist(strcat(model,'.m'),'file') || exist(strcat(model,'.',mexext),'file')
            model = str2func(filename);
        elseif exist(strcat(model,'.c'),'file')
            mexcompileIQM(filename);
            model = str2func(filename);
        elseif exist(strcat(model,'.sbml'),'file')
            sbm = IQMmodel(strcat(model,'.sbml'));
            switch fast_flag
                case {0, 1}
                    odefun_name = strcat(filename,'_odefun');
                    fprintf('Making %s.m ...',odefun_name);
                    IQMcreateODEfile(sbm,odefun_name);
                    fprintf(' Finished.\n');
                    model = str2func(odefun_name);
                case 2
                    mex_name = strcat(filename,'_mex');
                    clear(mex_name);
                    fprintf('Making %s.%s ...',mex_name,mexext);
                    IQMmakeMEXmodel(sbm,mex_name,1);
                    mexcompileIQM(mex_name);
                    fprintf(' Finished.\n');
                    model = str2func(mex_name);
                otherwise
                    error('Unexpected fast_flag!');
            end
        elseif exist(strcat(model,'.xml'),'file')
            sbm = IQMmodel(strcat(model,'.xml'));
            switch fast_flag
                case {0, 1}
                    odefun_name = strcat(filename,'_odefun');
                    fprintf('Making %s.m ...',odefun_name);
                    IQMcreateODEfile(sbm,odefun_name);
                    fprintf(' Finished.\n');
                    model = str2func(odefun_name);
                case 2
                    mex_name = strcat(filename,'_mex');
                    clear(mex_name);
                    fprintf('Making %s.%s ...',mex_name,mexext);
                    IQMmakeMEXmodel(sbm,mex_name,1);
                    mexcompileIQM(mex_name);
                    fprintf(' Finished.\n');
                    model = str2func(mex_name);
                otherwise
                    error('Unexpected fast_flag!');
            end
        else
            error('File %s does not exist!',model);
        end
    end

end


%% Checking whether model exists as a file
filename = func2str(model);
file_type = exist(filename,'file');
if  file_type == 0
    error('%s does not exist!',filename);
elseif file_type < 2 || 3 < file_type
    error('%s is neither an m file nor a MEX file!',filename);
end


%% Checking whether model returns parameter names
try
    output = feval(model,'parameters');
catch ME
    warning(ME.message);
    error('%s should return parameter names when ''parameters'' is given as an argument!',func2str(model));
end


%% Checking whether model returns the number of parameters
try
    output = feval(model,'parametervalues');
catch ME
    warning(ME.message);
    error('%s should return default parameter values when ''parametervalues'' is given as an argument!',func2str(model));
end


%% Number of parameters
n_param = length(output);


%% Setting Simulation function
switch fast_flag
    case 0
        Simulation = @RCGAsimulateODEXX;
    case 1
        Simulation = @RCGAsimulateSTB;
    case 2
        Simulation = @RCGAsimulateMEX;
    otherwise
        error('Unexpected fast_flag!');
end

filename = func2str(model);
exist_flag = exist(filename,'file');

if exist_flag == 2 && ~( fast_flag == 0 || fast_flag == 1 )
    warning('fast_flg set to 0');
    fast_flag = 0;
    Simulation = @RCGAsimulateODEXX;
end

if exist_flag == 3 && ~( fast_flag == 2 )
    warning('fast_flg set to 2');
    fast_flag = 2;
    Simulation = @RCGAsimulateMEX;
end


%% Reading measurement
if ~isIQMmeasurement(mst)
    if ischar(mst)
        fprintf('Reading %s ...',mst);
        mst = IQMmeasurement(mst);
        fprintf(' Finished.\n');
    else
        error('mst must be an IQMmeasurement or a file name!');
    end
end

if 1 < length(mst)
    warning('Provided experimental data include multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
if iscell(mst)
    mst = mst{1};
end


%% Prepering inputs for REXstarJGG
problem.n_gene = n_param;
problem.n_constraint = n_constraint;
problem.decodingfun = decodingfun;
problem.fitnessfun = @(x) fitnessfun(x,Simulation,model,mst,simopts);

if ~isfield(opts,'interimreportfun')
    opts.interimreportfun = @RCGAinterimreportfun_PE;
end
interimreportfun = opts.interimreportfun;
opts.interimreportfun = @(elapsedTime,generation,problem,opts,Population,best) ...
    interimreportfun(...
    elapsedTime,generation,problem,opts,Population,best,...
    Simulation,model,mst,simopts);

%% Decomment here if you want to use parameter names for out_transition, out_best, and out_population
% if ~isfield(opts,'finalreportfun')
%     opts.finalreportfun = @RCGAfinalreportfun_PE;
% end
% finalreportfun = opts.finalreportfun;
% opts.finalreportfun = @(elapsedTime,generation,problem,opts,Population,best) ...
%     finalreportfun(...
%     elapsedTime,generation,problem,opts,Population,best,...
%     Simulation,model,mst,simopts);


%% Run parameter estimation
Results = GenerationAlternation(problem,opts);


%% Print best
paramnames = model('parameters');

fprintf('\n--- Best parameter set (f = %e) ---\n',Results.Best.f);
for i = 1 : n_param
    fprintf('%s = %e\n',char(paramnames(i)),Results.Best.x(i));
end

