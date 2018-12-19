function f = mySSR(x, model, mst, mex_name)

if length(mst) > 2
    warning('Measurement has multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
mst = struct(mst{1});

if exist(mex_name,'file') == 3
    
    try
        output = feval(mex_name,mst.time,[],x');
    catch
        f = 1e+10;
        return;
    end
    
elseif isSBmodel(model)
    
    st_model = struct(model);
    [ ~, n_param ] = size(st_model.parameters);
    for i = 1: n_param
        st_model.parameters(i).value = x(i);
    end
    new_model = SBmodel(st_model);
    output = SBsimulate(new_model,mst.time);
    
else
    
    error('No MEX file or SBmodel provided!');
    
end


if max(max(isnan(output.statevalues)))
    f = 1e+10;
    return;
end

x_sim = output.statevalues;
x_exp = zeros(size(x_sim));
for i = 1 : length(mst.data)
    x_exp(:,i) = mst.data(i).values;
end


[n_row, n_col] = size(x_exp);
f = 0;
for i = 1 : n_row
    for j = 1 : n_col
        if ~isnan(x_exp(i,j))
%             f = f + abs( x_sim(i,j) - x_exp(i,j) );
%             f = f + abs( ( x_sim(i,j) - x_exp(i,j) ) / x_exp(i,j) );
            f = f + ( x_sim(i,j) - x_exp(i,j) ) ^ 2;
%             f = f + ( ( x_sim(i,j) - x_exp(i,j) ) / x_exp(i,j) ) ^ 2;
        end
    end
end
