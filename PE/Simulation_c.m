function  [ T, X ] = Simulation_c(x, tspan, mex_name, n_var, opts)

if exist(mex_name,'file') == 3

    try
        output = feval(mex_name,tspan,[],x',opts); % Alternative way of simulation
        T = output.time;
        X = output.statevalues;
    catch ME
        warning(ME.identifier);
        warning('Error in excuting mexed %s',mex_name);
        T = NaN;
        X = NaN(1,n_var);
    end
    
else
    error('No MEX file provided!');
    
end
