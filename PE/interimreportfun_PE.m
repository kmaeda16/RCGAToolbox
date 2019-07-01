function interimreportfun_PE(Simulation,elapsedTime,generation,problem,opts,Population,best,modelfun,mst,simopts,fast_flag)


%% Making outputs
printTransition(elapsedTime,generation,problem,best);
writeTransition(elapsedTime,generation,problem,opts,best);


%% Decoding gene to x
x = problem.decodingfun(best.gene);


%% Setting experimental data
if length(mst) > 2
    warning('Measurement has multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
mst = struct(mst{1});


%% Checking time errors
n_point = 100; % You can change this line

t0 = 0;
if t0 < mst.time(1)
    tspan = linspace(t0 ,mst.time(end),n_point)';
elseif t0 == mst.time(1)
    tspan = linspace(mst.time(1),mst.time(end),n_point)';
else
    error('Time of the first experimental datapoint should be AFTER or EQUAL TO time 0');
end


%% Running simulation
y0 = modelfun();
[ T, X ] = feval(Simulation, modelfun, tspan, y0, x', opts);


%% Making plot
for i = 1 : length(mst.data)
    x_exp(:,i) = mst.data(i).values;
end

statename = modelfun('states');

t_sim = T;
x_sim = X;
t_exp = mst.time;
plotter(generation,t_sim,x_sim,t_exp,x_exp,'Time','AU',statename);
