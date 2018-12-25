function interimreportfun_PE_ODEFUN(elapsedTime,generation,Param,Population,best,odefun,IC,mst,fast_flag)

printTransition(elapsedTime,generation,best);
writeTransition(elapsedTime,generation,Param,best);

x = Param.decodingfun(best.gene);

if length(mst) > 2
    warning('Measurement has multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
mst = struct(mst{1});

t0 = IC(x);
if t0 < mst.time(1)
    tspan = linspace(t0 ,mst.time(end),50);
elseif t0 == mst.time(1)
    tspan = linspace(mst.time(1),mst.time(end),50);
else
    error('t0 <= mst.time(1) must be satisfied!');
end

[ T, X ] = Simulation_ODEFUN(x, odefun, IC, tspan, fast_flag);

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
