function interimreportfun_PE(Simulation,elapsedTime,generation,problem,opts,Population,best,mex_name,mst,simopts,fast_flag)

n_point = 100;
printTransition(elapsedTime,generation,problem,best);
writeTransition(elapsedTime,generation,problem,opts,best);

x = problem.decodingfun(best.gene);

if length(mst) > 2
    warning('Measurement has multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
mst = struct(mst{1});

%% Checking time errors
t0 = 0;
if t0 < mst.time(1)
    tspan = linspace(t0 ,mst.time(end),n_point)';
elseif t0 == mst.time(1)
    tspan = linspace(mst.time(1),mst.time(end),n_point)';
else
    error('Time of the first experimental datapoint should be AFTER or EQUAL TO time 0');
end

%%
y0 = mex_name();
[ T, X ] = feval(Simulation, mex_name, tspan, y0, x', opts);

for i = 1 : length(mst.data)
    x_exp(:,i) = mst.data(i).values;
end

[~, n_state] = size(X);
for i = 1 : n_state;
    statename{i} = sprintf('state variable %d',i);
end

t_sim = T;
x_sim = X;
t_exp = mst.time;
plotter(generation,t_sim,x_sim,t_exp,x_exp,'Time','AU',statename);
