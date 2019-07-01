function [ T, Y ] = Simulation_odexx(odefun, tspan, y0, param, options)
% Simulation_odexx simulates an ODE model using MATLAB odexx solver.
% 
% [SYNTAX]
% [ T, Y ] = Simulation_odexx(odefun, param, tspan, options)
% 
% [INPUT]
% odefun :  ODEFUN file.
% tspan  :  [t0, tf] or [t0, t1, ..., tf].
% y0     :  Initial value vector.
% param  :  Parameter value vector.
% options:  Options for odexx solvers. Function odeset can be used for
%           setting options. The field 'method' specifies which solver you
%           use. The default is 'ode23s'. For other fields, see MATLAB
%           documentation for oddeset.
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

odefun_temp = @(t,y) odefun(t,y,param);


%% Setting solver
if isfield(options,'method')
    method = options.method;
else
    method = 'ode23s';
end


%% Solving ODEs
try
    [T, Y] = feval(method,odefun_temp,tspan,y0,options);
catch
    warning('Error in %s.',method);
    T = NaN;
    Y = NaN(1,length(y0));
end
    