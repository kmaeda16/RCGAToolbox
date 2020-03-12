function [ T, Y ] = Simulation_stb(odefun, tspan, y0, param, options)
% Simulation_stb simulates an ODE model using CVODE from SundialsTB.
% 
% [SYNTAX]
% [ T, Y ] = Simulation_stb(odefun, tspan, y0)
% [ T, Y ] = Simulation_stb(odefun, tspan, y0, param)
% [ T, Y ] = Simulation_stb(odefun, tspan, y0, [], options)
% [ T, Y ] = Simulation_stb(odefun, tspan, y0, param, options)
% 
% [INPUT]
% odefun :  ODEFUN file.
% tspan  :  [t0, tf] or [t0, t1, ..., tf].
% y0     :  Initial value vector.
% param  :  Parameter value vector. If param is not given, default
%           parameter values provided in odefun will be used.
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

options.Method = @odestb;
[ T, Y ] = Simulation_odexx(odefun, tspan, y0, param, options);
