function [ T, Y ] = Simulation_stb(odefun, tspan, y0, param, options)
% Simulation_stb simulates an ODE model using CVODE from SundialsTB.
% 
% [SYNTAX]
% [ T, Y ] = Simulation_odexx(odefun, param, tspan, options)
% 
% [INPUT]
% odefun :  ODEFUN file.
% tspan  :  [t0, tf] or [t0, t1, ..., tf].
% y0     :  Initial value vector.
% param  :  Parameter value vector.
% options:  Structure with integrator options.
%           * options.LMM: Linear Multistep Method (default: 'BDF')
%           * options.NonlinearSolver: Type of nonlinear solver used
%             (default: 'Newton')
%           * options.AbsTol:  Absolute tolerance (default: 1e-6)
%           * options.RelTol:  Relative tolerance (default: 1e-4)
%           * options.MinStep: Minimum stepsize (default: 0)
%           * options.MaxStep: Maximum stepsize (default: inf)
%           * options.MaxNumSteps: Maximum number of steps (default: 500)
%           For other fields, see SundialsTB documentation.
% 
% [OUTPUT]
% T      :  Column vector of timepoints
% Y      :  Variable matrix. Each column corresponds to each variable. 
%           Each row of Y corresponds to each row of T. 


%% Checking file errors 
odefun_name = func2str(odefun);
flg = exist(odefun_name,'file');

if flg == 0
    error('File "%s" does NOT exist!',odefun_name);
elseif flg ~= 2
    error('%s is NOT an m file',odefun_name);
end


%% Setting the initial condition and parameters
if isempty(tspan)
    tspan = [0 10];
end
if isempty(y0)
    y0 = feval(odefun);
end
if isempty(param)
    param = feval(odefun,'parametervalues');
end

odefun_temp = @(t,y) wrapper_odefun(odefun,t,y,param);


%% Setting solver and initializing CVode
if isempty(fieldnames(options))
    CVodeInit(odefun_temp, 'BDF', 'Newton', tspan(1), y0);
else
    if isfield(options,'LMM')
        LMM = options.LMM;
    else
        LMM = 'BDF';
    end
    if isfield(options,'NonlinearSolver')
        NonlinearSolver = options.NonlinearSolver;
    else
        NonlinearSolver = 'Newton';
    end
    options = CVodeSetOptions(options,'LMM',LMM,'NonlinearSolver',NonlinearSolver);
    CVodeInit(odefun_temp, LMM, NonlinearSolver, tspan(1), y0, options);
end


%% Solving ODEs
try
    [ ~, T, Y ] = CVode(tspan(2:end),'Normal');
    T = [tspan(1) T]';
    Y = [y0 Y]';
catch
    warning('Error in CVode.');
    T = NaN;
    Y = NaN(1,length(y0));
end


%% Deinitializing CVode
CVodeFree;
