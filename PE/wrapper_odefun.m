function [ dxdt, flag, new_data ] = wrapper_odefun(t,x,p,odefun)

dxdt = odefun(t,x,p);
flag = 0;
new_data = [];
