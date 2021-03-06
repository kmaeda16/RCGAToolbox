function [ T, Y ] = RCGAsimulateMEX(mexfun, tspan, y0, param, options)
% RCGAsimulateMEX simulates a mexfun (a MEXed ODE model) using using CVODE
% from IQM Tools. mexfun can be created by RCGAmakeMEXmodel,
% IQMmakeMEXmodel, or mexcompileIQM.
% 
% [SYNTAX]
% [ T, Y ] = RCGAsimulateMEX(mexfun)
% [ T, Y ] = RCGAsimulateMEX(mexfun, tspan)
% [ T, Y ] = RCGAsimulateMEX(mexfun, tspan, y0)
% [ T, Y ] = RCGAsimulateMEX(mexfun, tspan, y0, param)
% [ T, Y ] = RCGAsimulateMEX(mexfun, tspan, y0, [], options)
% [ T, Y ] = RCGAsimulateMEX(mexfun, tspan, y0, param, options)
% 
% [INPUT]
% mexfun    :  Function handle for a MEXed ODE model.
% tspan     :  [t0, tf] or [t0, t1, ..., tf] (default: [0 10]).
% y0        :  Initial value vector (default: Values stored in mexfun).
% param     :  Parameter value vector (default: Values stored in mexfun).
% options   :  Solver option structure:
%              - options.abstol: Absolute tolerance (default: 1e-6).
%              - options.reltol: Relative tolerance (default: 1e-6).
%              - options.minstep: Minimum step-size of integrator (default: 
%                 0).
%              - options.maxstep: Maximum step-size of integrator (default:
%                 inf).
%              - options.maxnumsteps: Maximum number of steps to be taken 
%                 by the solver in its attempt to reach the next output 
%                 time (default: 100000).
%              For other fields, see IQM Tools documentation.
% 
% [OUTPUT]
% T         :  Column vector of timepoints.
% Y         :  Variable matrix. Each column corresponds to each variable. 
%              Each row corresponds to each timepoint.


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
mexfun_name = func2str(mexfun);
flg = exist(mexfun_name,'file');

if flg == 0
    error('File "%s" does NOT exist!',mexfun_name);
elseif flg ~= 3
    error('%s is NOT a MEX file.',mexfun_name);
end


%% Setting the initial condition and parameters
if isempty(tspan)
    tspan = 0:0.1:10;
end
if isempty(y0)
    y0 = feval(mexfun);
end
if isempty(param)
    param = feval(mexfun,'parametervalues');
end
if length(param) ~= length(feval(mexfun,'parametervalues'))
    error('%s has %d parameters, but %d parameters were provided to %s!',...
        func2str(mexfun),length(feval(mexfun,'parametervalues')),...
        length(param),mfilename);
end

% tspan must be a column vector to get T as a column vector
[n_row, n_col] = size(tspan);
if n_col > n_row
    tspan = tspan';
end


%% Running MEX
try
    output = feval(mexfun,tspan,y0,param,options);
    T = output.time;
    Y = output.statevalues;
catch ME
    warning('%s',ME.message);
    T = NaN;
    Y = NaN(1,length(y0));
end
