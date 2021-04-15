function S = RCGAsensitivity(model, mst, param, fast_flag, simopts, fitnessfun, delta, norm_flag)
% RCGAsensitivity calculates parameter sensitivity for the objective 
% function (f) provided in fitnessfun.
% 
% [SYNTAX]
% S = RCGAsensitivity(model, mst, param)
% S = RCGAsensitivity(model, mst, param, fast_flag)
% S = RCGAsensitivity(model, mst, param, fast_flag, simopts)
% S = RCGAsensitivity(model, mst, param, fast_flag, simopts, fitnessfun)
% S = RCGAsensitivity(model, mst, param, fast_flag, simopts, fitnessfun, delta)
% S = RCGAsensitivity(model, mst, param, fast_flag, simopts, fitnessfun, delta, norm_flag)
% 
% [INPUT]
% model      :  An IQMmodel object, the name of SBML file, the function 
%               handle for an ODE function (IQM Tools format), the function 
%               handle for a MEXed model, or the C source code.
% mst        :  Experimental data (the name of an IQM measurement file or 
%               an IQMmeasurement object).
% param      :  Parameter value vector. Sensitivity analysis is performed
%               based on these parameter values.
% fast_flag  :  ODE solver flag:
%               - fast_flag = 0: ODEXX by MATLAB built-ins.
%               - fast_flag = 1: CVODE by SundialsTB.
%               - fast_flag = 2: CVODE by IQM Tools.
% simopts    :  Solver option structure. The fields depend on fast_flag. 
%               For fast_flag = 0, 1, and 2, see 'help RCGAsimulateODE', 
%               'help RCGAsimulateSTB', 'help RCGAsimulateMEX', 
%               respectively.
% fitnessfun :  Function handle for a fitness function. The default is 
%               RCGAssr function:
%                f = RCGAssr(param, Simulation, model, mst, simopts)
%               RCGAssr can be used as a template for a user-defined 
%               fitness function.
% delta      :  Magnitude of perturbation used for sensitivity calculation.
%               Default: delta = 0.01.
% norm_flag  :  Normalization flag for sensitivity calculation. 
%               Default: norm_flag = 3.
%               - norm_flag = 0: No normalization for both f and p.
%                  Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std )
%               - norm_flag = 1: Normalization only for p.
%                  Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std ) * param_std
%               - norm_flag = 2: Normalization only for f.
%                  Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std ) / f_std
%               - norm_flag = 3: Normalization for both f and p.
%                  Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std ) / f_std * param_std
% 
% [OUTPUT]
% S          :  Sensitivity table (n_param x 1). The ith element indicates 
%               the objective function (f) sensitivity to changes in the 
%               ith parameter.


%% Handling inputs
if ~exist('model','var') || isempty(model)
    error('model is not provided.');
end
if ~exist('mst','var') || isempty(mst)
    error('mst is not provided.');
end
if ~exist('param','var') || isempty(param)
    error('param is not provided.');
end
if ~exist('fast_flag','var') || isempty(fast_flag)
    fast_flag = 0;
end
if ~exist('simopts','var') || isempty(simopts)
    simopts = [];
end
if ~exist('fitnessfun','var') || isempty(fitnessfun)
    fitnessfun = @RCGAssr;
end
if ~exist('delta','var') || isempty(delta)
    delta = 0.01;
end
if ~exist('norm_flag','var') || isempty(norm_flag)
    norm_flag = 3;
end


%% Checking if IQM Tools are available.
if exist('isIQMmodel','file') == 0 || ...
        exist('IQMcreateODEfile','file') == 0 || ...
        exist('IQMmakeMEXmodel','file') == 0 || ...
        exist('mexcompileIQM','file') == 0 || ...
        exist('IQMmodel','file') == 0 || ...
        exist('isIQMmeasurement','file') == 0 || ...
        exist('IQMmeasurement','file') == 0
    warning('IQM Tools are not properly installed. Run the script RCGAToolbox/install/RCGAToolbox_Diagnosis for diagnosis.');
end


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
        Simulation = @RCGAsimulateODE;
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
    Simulation = @RCGAsimulateODE;
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


%% Calculating parameter sensitivities
if length(param) ~= n_param
    error('The provided param has %d parameters, but the provided model has %d parameters.',length(param),n_param);
end

param_std = param;

f_std = fitnessfun(param_std, Simulation, model, mst, simopts);

Svec = zeros(n_param,1);

for i = 1 : n_param
    
    param_ptb = param_std;
    param_ptb(i) = param_std(i) * ( 1 + delta );
    f_ptb = fitnessfun(param_ptb, Simulation, model, mst, simopts);
    
    switch norm_flag
        case 0
            Svec(i) = ( f_ptb - f_std ) / ( param_ptb(i) - param_std(i) ); % No normalization for both f and p
        case 1
            Svec(i) = ( f_ptb - f_std ) / ( param_ptb(i) - param_std(i) ) * param_std(i); % Normalization only for p
        case 2
            Svec(i) = ( f_ptb - f_std ) / ( param_ptb(i) - param_std(i) ) / f_std; % Normalization only for f
        case 3
            Svec(i) = ( f_ptb - f_std ) / ( param_ptb(i) - param_std(i) ) * param_std(i) / f_std ; % Normalization for both f and p
        otherwise
            error('Unexpected norm_flag!');
    end
    
end


%% Making the output table S
paramnames = model('parameters');
S = table(Svec,'VariableNames',{'Sensitivity_f'},'RowNames',paramnames);
