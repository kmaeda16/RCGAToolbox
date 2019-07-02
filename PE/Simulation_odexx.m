function [ T, Y ] = Simulation_odexx(odefun, tspan, y0, param, options)
% Simulation_odexx simulates an ODE model using MATLAB odexx solver.
% 
% [SYNTAX]
% [ T, Y ] = Simulation_odexx(odefun, tspan, y0)
% [ T, Y ] = Simulation_odexx(odefun, tspan, y0, param)
% [ T, Y ] = Simulation_odexx(odefun, tspan, y0, [], options)
% [ T, Y ] = Simulation_odexx(odefun, tspan, y0, param, options)
% 
% [INPUT]
% odefun :  ODEFUN file.
% tspan  :  [t0, tf] or [t0, t1, ..., tf].
% y0     :  Initial value vector.
% param  :  Parameter value vector. If param is not given, default
%           parameter values provided in odefun will be used.
% options:  Structure with integrator options.
%           * options.Method: Solver (default: 'ode15s')
%           * options.AbsTol: Absolute tolerance (default: 1e-6)
%           * options.RelTol: Relative tolerance (default: 1e-3)
%           * options.MaxStep: Maximum step-size of integrator (default: 0.1*abs(t0-tf))
%           For other fields, see MATLAB documentation for oddeset.
% 
% [OUTPUT]
% T      :  Column vector of timepoints
% Y      :  Variable matrix. Each column corresponds to each variable. 
%           Each row of Y corresponds to each row of T. 


%% Handling inputs
if nargin == 3
    param = [];
    options = [];
end
if nargin == 4
    options = [];
end


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
if length(param) ~= length(feval(odefun,'parametervalues'))
    error('%s has %d parameters, but %d parameters were provided to %s!',...
        func2str(odefun),length(feval(odefun,'parametervalues')),...
        length(param),mfilename);
end

odefun_temp = @(t,y) odefun(t,y,param);


%% Setting solver
if isfield(options,'Method')
    method = options.Method;
elseif isfield(options,'method')
    method = options.method;
else
    method = 'ode15s';
end


%% Solving ODEs
try
    [T, Y] = feval(method,odefun_temp,tspan,y0,options);
catch ME
    warning('%s',ME.message);
    T = NaN;
    Y = NaN(1,length(y0));
end
    