function  [ T, X ] = Simulation_sbml(x, model, tspan, mex_name, opts)

if exist(mex_name,'file') == 3

    try
        st_model = struct(model);
        [ ~, n_param ] = size(st_model.parameters);
        for i = 1 : n_param
            param_name{i} = st_model.parameters(i).name;
        end
        output = IQMPsimulate(mex_name,tspan,[],param_name,x,opts);
%         output = feval(mex_name,tspan,[],x',opts); % Alternative way of simulation
        T = output.time;
        X = output.statevalues;
    catch ME
        warning(ME.identifier);
        warning('Error in IQMPsimulate');
        T = NaN;
        X = NaN(1,length(st_model.states));
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
        X = NaN(1,length(st_model.states));
    end
    
else
    
    error('No MEX file or IQMmodel provided!');
    
end
