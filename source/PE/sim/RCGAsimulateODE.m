function [ T, Y ] = RCGAsimulateODE(odefun, tspan, y0, param, options)
% RCGAsimulateODE simulates an odefun (IQM Tools format) using MATLAB
% odexx solver. odefun (IQM Tools format) can be created by
% RCGAcreateODEfile or IQMcreateODEfile.
% 
% [SYNTAX]
% [ T, Y ] = RCGAsimulateODE(odefun)
% [ T, Y ] = RCGAsimulateODE(odefun, tspan)
% [ T, Y ] = RCGAsimulateODE(odefun, tspan, y0)
% [ T, Y ] = RCGAsimulateODE(odefun, tspan, y0, param)
% [ T, Y ] = RCGAsimulateODE(odefun, tspan, y0, [], options)
% [ T, Y ] = RCGAsimulateODE(odefun, tspan, y0, param, options)
% 
% [INPUT]
% odefun  :  Function handle for an ODE function (IQM Tools format).
% tspan   :  [t0, tf] or [t0, t1, ..., tf] (default: [0 10]).
% y0      :  Initial value vector (default: Values stored in odefun).
% param   :  Parameter value vector (default: Values stored in odefun).
% options :  Solver option structure:
%            - options.Method: Solver (default: 'ode23s').
%            - options.AbsTol: Absolute tolerance (default: 1e-6).
%            - options.RelTol: Relative tolerance (default: 1e-3).
%            - options.MaxStep: Maximum step-size of integrator (default:
%               0.1 * abs(t0 - tf)).
%            For other fields, see MATLAB documentation for odeset.
% 
% [OUTPUT]
% T       :  Column vector of timepoints.
% Y       :  Variable matrix. Each column corresponds to each variable. 
%            Each row corresponds to each timepoint.


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
if ~exist('options','var')
    options = [];
end


%% Checking file errors 
odefun_name = func2str(odefun);
flg = exist(odefun_name,'file');

if flg == 0
    error('File "%s" does NOT exist!',odefun_name);
elseif flg ~= 2
    error('%s is NOT an m file.',odefun_name);
end


%% Setting the initial condition and parameters
if isempty(tspan)
    tspan = 0:0.1:10;
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
    method = 'ode23s';
end


%% Solving ODEs
try
    [T, Y] = feval(method,odefun_temp,tspan,y0,options);
catch ME
    warning('%s',ME.message);
    T = NaN;
    Y = NaN(1,length(y0));
end
    