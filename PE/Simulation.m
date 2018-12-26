function  [ T, X ] = Simulation(x, model, tspan, mex_name, opts)


if exist(mex_name,'file') == 3
    
    st_model = struct(model);
    [ ~, n_param ] = size(st_model.parameters);
    for i = 1 : n_param
        param_name{i} = st_model.parameters(i).name;
    end
    
    try
        output = feval(mex_name,tspan,[],x',opts);
        T = output.time;
        X = output.statevalues;
    catch
        warning('Error in MEXed ODEs.');
        T = NaN;
        X = NaN(1,length(x0));
    end
    
elseif isIQMmodel(model)
    
    st_model = struct(model);
    [ ~, n_param ] = size(st_model.parameters);
    for i = 1 : n_param
        st_model.parameters(i).value = x(i);
    end
    new_model = IQMmodel(st_model);
    if isfield(opts,'method')
        method = opts.method;
    else
        method = 'ode23s';
    end
    
    try
        output = IQMsimulate(new_model,method,tspan,[],opts);
        T = output.time;
        X = output.statevalues;
    catch
        warning('Error in IQMsimulate.');
        T = NaN;
        X = NaN(1,length(x0));
    end
    
else
    
    error('No MEX file or IQMmodel provided!');
    
end
