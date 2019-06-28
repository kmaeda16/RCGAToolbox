function [ T, X ] = Simulation_odexx(x, odefun, tspan, opts)

odefun_name = func2str(odefun);

flg = exist(odefun_name,'file');

if flg == 0
    error('File "%s" does NOT exist!',odefun_name);
elseif flg ~= 2
    error('%s is NOT an m file',odefun_name);
end

t0 = 0;
x0 = feval(odefun);
param = x';

if ~( t0 <= tspan(1) )
    error('t0 <= mst.time(1) must be satisfied!');
end
if t0 < tspan(1)
    tspan = [ t0;  tspan ];
end


odefun_temp = @(t,y) odefun(t,y,param);
if isfield(opts,'method')
    method = opts.method;
else
   % method = 'ode15s';
    method = 'ode23s';
end

try
    [T, X] = feval(method,odefun_temp,tspan,x0,opts);
catch
    warning('Error in %s.',method);
    T = NaN;
    X = NaN(1,length(x0));
end
    