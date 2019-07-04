function [ T, Y ] = RCGAsimulate(model, tspan, y0, param, fast_flag, options)
% RCGAsimulate simulates model.
% 
% [SYNTAX]
% [ T, Y ] = RCGAsimulate(odefun, tspan, y0)
% [ T, Y ] = RCGAsimulate(odefun, tspan, y0, param)
% [ T, Y ] = RCGAsimulate(odefun, tspan, y0, [], options)
% [ T, Y ] = RCGAsimulate(odefun, tspan, y0, param, options)
% 
% [INPUT]
% odefun :  IQMmodel, File names of SBML, or ODE function, or MEX files.
%           File function of ODE function or MEX function.
% tspan  :  [t0, tf] or [t0, t1, ..., tf]. (default: [0 10])
% y0     :  Initial value vector. (default: Values stored in the model)
% param  :  Parameter value vector. (default: Values stored in the model)
% options:  Structure with integrator options. Fields depend on
%           Simulation_*. See 'help Simulation_*'.
% 
% [OUTPUT]
% T      :  Column vector of timepoints
% Y      :  Variable matrix. Each column corresponds to each variable. 
%           Each row of Y corresponds to each row of T. 


%% Handling inputs
if ~exist('tspan','var')
    tspan = [];
end
if ~exist('y0','var')
    y0 = [];
end
if ~exist('param','var')
    param = [];
end
if ~exist('fast_flag','var')
    fast_flag = 0;
end
if ~exist('options','var')
    options = [];
end


%% If model is an IQMmodel, it is converted into odefun.m or odefun.c.
%  After this section, model will be a string '*_odefun' or '*_mex'.
if isIQMmodel(model)
    st_model = struct(model);
    switch fast_flag
        case {0, 1}
            odefun_name = strcat(st_model.name,'_odefun');
            fprintf('Making %s.m ...',odefun_name);
            IQMcreateODEfile(model,odefun_name);
            fprintf(' Finished.\n');
            model = odefun_name;
        case 2
            mex_name = strcat(st_model.name,'_mex');
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
        st_model = struct(sbm);
        switch fast_flag
            case {0, 1}
                odefun_name = strcat(st_model.name,'_odefun');
                fprintf('Making %s.m ...',odefun_name);
                IQMcreateODEfile(sbm,odefun_name);
                fprintf(' Finished.\n');
                model = str2func(odefun_name);
            case 2
                mex_name = strcat(st_model.name,'_mex');
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
            st_model = struct(sbm);
            switch fast_flag
                case {0, 1}
                    odefun_name = strcat(st_model.name,'_odefun');
                    fprintf('Making %s.m ...',odefun_name);
                    IQMcreateODEfile(sbm,odefun_name);
                    fprintf(' Finished.\n');
                    model = str2func(odefun_name);
                case 2
                    mex_name = strcat(st_model.name,'_mex');
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
            st_model = struct(sbm);
            switch fast_flag
                case {0, 1}
                    odefun_name = strcat(st_model.name,'_odefun');
                    fprintf('Making %s.m ...',odefun_name);
                    IQMcreateODEfile(sbm,odefun_name);
                    fprintf(' Finished.\n');
                    model = str2func(odefun_name);
                case 2
                    mex_name = strcat(st_model.name,'_mex');
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


%% Setting Simulation function
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


%%
[ T, Y ] = feval(Simulation, model, tspan, y0, param, options);
