function [ T, Y ] = odestb(odefun, tspan, y0, options)
% odestb is the wrapper function that enables to use CVode provided by
% SundialsTB in the same way as MATLAB ODE solvers.
% 
% [SYNTAX]
% [ T, Y ] = odestb(odefun, tspan, y0)
% [ T, Y ] = odestb(odefun, tspan, y0, options)
% 
% [INPUT]
% odefun  :  Function handle for odefun (an ODE function).
% tspan   :  [t0, tf] or [t0, t1, ..., tf].
% y0      :  Initial value vector.
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


%% Checking if SundialsTB is available.
if exist('CVodeInit','file') == 0 || ...
        exist('CVodeSetOptions','file') == 0 || ...
        exist('CVodeInit','file') == 0 || ...
        exist('CVode','file') == 0 || ...
        exist('CVodeFree','file') == 0
    warning('SundialsTB is not properly installed. Run the script RCGAToolbox/install/RCGAToolbox_Diagnosis for diagnosis.');
end


%% Handling inputs
if ~exist('tspan','var')
    error('tspan should be provided!');
end
if ~exist('y0','var')
    error('y0 should be provided!');
end
if ~exist('options','var')
    options = [];
end


%% Setting the initial condition
% y0 must be a column vector for CVodeInit
[n_row, n_col] = size(y0);
if n_col > n_row
    y0 = y0';
end

if length(tspan) <= 1
    error('tspan must have at least two elements: start time and end time.');
elseif length(tspan) == 2
    tspan = linspace(tspan(1),tspan(2),100);
end

odefun_temp = @(t,y) wrapper_odefun(odefun,t,y);


%% Setting solver and initializing CVode
if isempty(options) || isempty(fieldnames(options))
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
    T = nan;
    Y = nan;
end


%% Deinitializing CVode
CVodeFree;
