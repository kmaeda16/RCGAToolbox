function [ T, X ] = Simulation_mex(x, mexfun, tspan, opts)

mexfun_name = func2str(mexfun);

flg = exist(mexfun_name,'file');

if flg == 0
    error('File "%s" does NOT exist!',mexfun_name);
elseif flg ~= 3
    error('%s is NOT a MEX file',mexfun_name);
end

t0 = 0;
x0 = feval(mexfun);
param = x';

if ~( t0 <= tspan(1) )
    error('t0 <= mst.time(1) must be satisfied!');
end
if t0 < tspan(1)
    tspan = [ t0;  tspan ];
end


try
    output = feval(mexfun,tspan,[],param,opts);
    T = output.time;
    X = output.statevalues;
catch
    warning('Error in excuting mexed %s',mexfun_name);
    T = NaN;
    X = NaN(1,length(x0));
end
