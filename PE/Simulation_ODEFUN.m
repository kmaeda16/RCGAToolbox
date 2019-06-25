function  [ T, X ] = Simulation_odefun(x, odefun, icfun, tspan, fast_flag, opts)

[ t0, x0 ] = icfun(x);
p = x;

if t0 < tspan(1)
    tspan = [ t0  tspan ];
elseif t0 == tspan(1)
else
    error('t0 <= mst.time(1) must be satisfied!');
end

if fast_flag == 0
    
    odefun_temp = @(t,x) odefun(t,x,p);
    if isfield(opts,'method')
        method = opts.method;
    else
%         method = 'ode15s';
        method = 'ode23s';
    end
    
    try
        [ T, X ] = feval(method,odefun_temp,tspan,x0,opts);
    catch
        warning('Error in %s.',method);
        T = NaN;
        X = NaN(1,length(x0));
    end
    
else
    
    odefun_temp = @(t,x) wrapper_odefun(t,x,p,odefun);
    if isempty(fieldnames(opts))
        CVodeInit(odefun_temp, 'BDF', 'Newton', tspan(1), x0);
    else
        if isfield(opts,'LMM')
            LMM = opts.LMM;
        else
            LMM = 'BDF';
        end
        if isfield(opts,'NonlinearSolver')
            NonlinearSolver = opts.NonlinearSolver;
        else
            NonlinearSolver = 'Newton';
        end
        opts = CVodeSetOptions(opts,'LMM',LMM,'NonlinearSolver',NonlinearSolver);
        CVodeInit(odefun_temp, LMM, NonlinearSolver, tspan(1), x0, opts);
    end
    
    try
        [ ~, T, X ] = CVode(tspan(2:end),'Normal');
        T = [t0 T]';
        X = [x0 X]';
    catch
        warning('Error in CVode.');
        T = NaN;
        X = NaN(1,length(x0));
    end
    CVodeFree;
    
end
