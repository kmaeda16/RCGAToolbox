function [ dydt, flag ] = wrapper_odefun(odefun,t,y)
% wrapper_odefun is a wrapper for odefun.
% 
% [SYNTAX]
% [ dxdt, flag ] = wrapper_odefun(odefun,t,y)
% 
% [INPUT]
% odefun :  ODEFUN file.
% t      :  Time.
% y      :  State variable vector.
% 
% [OUTPUT]
% dydt   :  Time derivative of y.
% flag   :  Error flag.
%           * flag = 0: Successful
%           * flag < 0: Unrecoverable failure
%           * flag > 0: Recoverable error


dydt = feval(odefun,t,y);
flag = 0;
