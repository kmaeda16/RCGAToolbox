function [ dydt, flag ] = wrapper_odefun(odefun,t,y,param)
% wrapper_odefun is a wrapper for odefun.
% 
% [SYNTAX]
% [ dxdt, flag ] = wrapper_odefun(odefun,t,y)
% [ dxdt, flag ] = wrapper_odefun(odefun,t,y,param)
% 
% [INPUT]
% odefun :  IQM Tools ODEFUN file.
% t      :  Time.
% y      :  State variable vector.
% param  :  Parameter value vector.
% 
% [OUTPUT]
% dydt      :  Time derivative of y.
% flag      :  Error flag.
%              * flag = 0: Successful
%              * flag < 0: Unrecoverable failure
%              * flag > 0: Recoverable error

dydt = feval(odefun,t,y,param);
flag = 0;
