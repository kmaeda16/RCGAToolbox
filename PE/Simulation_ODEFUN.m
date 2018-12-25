function  [ T, X ] = Simulation_ODEFUN(x, odefun, IC, tspan, fast_flag, opts)

[ t0, x0 ] = IC(x);
p = x;

if t0 < tspan(1)
    tspan = [ t0  tspan ];
elseif t0 == tspan(1)
else
    error('t0 <= mst.time(1) must be satisfied!');
end

if fast_flag == 0
    
    odefun_temp = @(t,x) odefun(t,x,p);
    [ T, X ] = ode15s(odefun_temp,tspan,x0,opts);
    
else
    
    odefun_temp = @(t,x) wrapper_odefun(t,x,p,odefun);
    if isempty(fieldnames(opts))
        CVodeInit(odefun_temp, 'BDF', 'Newton', tspan(1), x0);
    else
        CVodeInit(odefun_temp, 'BDF', 'Newton', tspan(1), x0, opts);
    end
    [ status, T, X ] = CVode(tspan(2:end),'Normal');
    CVodeFree;
    T = [t0 T]';
    X = [x0 X]';
    
end
