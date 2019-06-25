function  [ T, X ] = Simulation_sbml(x, mex_name, tspan, opts)

t0 = 0;
x0 = feval(mex_name);
p = x;

if t0 < tspan(1)
    tspan = [ t0  tspan ];
elseif t0 == tspan(1)
else
    error('t0 <= mst.time(1) must be satisfied!');
end

if exist(mex_name,'file') == 3

    try
        output = feval(mex_name,tspan,[],p',opts);
        T = output.time;
        X = output.statevalues;
    catch ME
        warning(ME.identifier);
        warning('Error in IQMPsimulate');
        T = NaN;
        X = NaN(1,length(x0));
    end
    
elseif exist(mex_name,'file') == 2
    
    odefun = @(t,y) odefun_temp(t,y,p);
    if isfield(opts,'method')
        method = opts.method;
    else
%         method = 'ode15s';
        method = 'ode23s';
    end
    
    try
        [T, X] = feval(method,odefun,tspan,x0,opts);
    catch
        warning('Error in %s.',method);
        T = NaN;
        X = NaN(1,length(x0));
    end
    
else
    
    error('No MEX file or odefun file provided!');
    
end
