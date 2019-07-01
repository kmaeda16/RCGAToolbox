function [ dxdt, flag, new_data ] = wrapper_odefun(odefun,t,y,param)
% [ dxdt, flag, new_data ] = wrapper_odefun(odefun,t,y)
% [ dxdt, flag, new_data ] = wrapper_odefun(odefun,t,y,param)

% switch nargin
%     case 3
%         dxdt = feval(odefun,t,y);
%     case 4
%         dxdt = feval(odefun,t,y,param);
%     otherwise
%         error('Unexpected number of inputs for wrapper_odefun!');
% end
dxdt = feval(odefun,t,y,param);
flag = 0;
new_data = [];
