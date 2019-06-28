function [ T, X ] = Simulation_stb(x, odefun, tspan, opts)

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


odefun_temp = @(t,x) wrapper_odefun(t,x,param,odefun);
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
    