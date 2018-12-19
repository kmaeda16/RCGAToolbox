function interimreportfun_PE(x,f,model,mst,mex_name)

if length(mst) > 2
    warning('Measurement has multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
mst = struct(mst{1});
time = linspace(mst.time(1),mst.time(end),100);

st_model = struct(model);

if exist(mex_name,'file') == 3
    
    try
        output = feval(mex_name,time,[],x');
    catch
%         f = 1e+10;
        return;
    end
    
elseif isSBmodel(model)
    
    
    [ ~, n_param ] = size(st_model.parameters);
    for i = 1: n_param
        st_model.parameters(i).value = x(i);
    end
    new_model = SBmodel(st_model);
    output = SBsimulate(new_model,time);
    
else
    
    error('No MEX file or SBmodel provided!');
    
end


x_sim = output.statevalues;
% x_exp = zeros(mst.data);
for i = 1 : length(mst.data)
    x_exp(:,i) = mst.data(i).values;
end

[~, n_state] = size(output.statevalues);
if isSBmodel(model)
    for i = 1 : n_state;
        statename{i} = st_model.states(i).name;
    end
else
    for i = 1 : n_state;
        statename{i} = sprintf('state variable %d',i);
    end
end

plot(mst.time,x_exp,'o');
legend(statename);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(output.time,x_sim);
hold off;

xlim([time(1) time(end)]);
ylim([0 max(max(x_exp))]);

drawnow;

