function  [ T, X ] = Simulation(x, model, tspan, mex_name, opts)


if exist(mex_name,'file') == 3
    
    st_model = struct(model);
    [ ~, n_param ] = size(st_model.parameters);
    for i = 1 : n_param
        param_name{i} = st_model.parameters(i).name;
    end
    
    try
%         output = feval(mex_name,tspan,[],x');
        output = SBPDsimulate(mex_name,tspan,[],param_name,x',opts);
    catch
        T = NaN;
        X = NaN;
        return;
    end
    
elseif isSBmodel(model)
    
    st_model = struct(model);
    [ ~, n_param ] = size(st_model.parameters);
    for i = 1 : n_param
        st_model.parameters(i).value = x(i);
    end
    new_model = SBmodel(st_model);
    output = SBsimulate(new_model,'ode15s',tspan,[],opts);
    
else
    
    error('No MEX file or SBmodel provided!');
    
end

T = output.time;
X = output.statevalues;
