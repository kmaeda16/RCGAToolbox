function [ T, Y ] = RCGAsimulateSTB(odefun, tspan, y0, param, options)
% RCGAsimulateSTB simulates an odefun (IQM Tools format) using CVODE from
% SundialsTB. odefun (IQM Tools format) can be created by RCGAcreateODEfile
% or IQMcreateODEfile.
% 
% [SYNTAX]
% [ T, Y ] = RCGAsimulateSTB(odefun)
% [ T, Y ] = RCGAsimulateSTB(odefun, tspan)
% [ T, Y ] = RCGAsimulateSTB(odefun, tspan, y0)
% [ T, Y ] = RCGAsimulateSTB(odefun, tspan, y0, param)
% [ T, Y ] = RCGAsimulateSTB(odefun, tspan, y0, [], options)
% [ T, Y ] = RCGAsimulateSTB(odefun, tspan, y0, param, options)
% 
% [INPUT]
% odefun  :  Function handle for an ODE function (IQM Tools format).
% tspan   :  [t0, tf] or [t0, t1, ..., tf] (default: [0 10]).
% y0      :  Initial value vector (default: Values stored in odefun).
% param   :  Parameter value vector (default: Values stored in odefun).
% options :  Solver option structure:
%            - options.LMM: Linear Multistep Method (default: 'BDF').
%            - options.NonlinearSolver: Type of nonlinear solver used
%               (default: 'Newton').
%            - options.AbsTol:  Absolute tolerance (default: 1e-6).
%            - options.RelTol:  Relative tolerance (default: 1e-4).
%            - options.MinStep: Minimum stepsize (default: 0).
%            - options.MaxStep: Maximum stepsize (default: inf).
%            - options.MaxNumSteps: Maximum number of steps (default: 500).
%            For other fields, see SundialsTB documentation.
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


%% Solving ODEs
options.Method = @odestb;
[ T, Y ] = RCGAsimulateODE(odefun, tspan, y0, param, options);
