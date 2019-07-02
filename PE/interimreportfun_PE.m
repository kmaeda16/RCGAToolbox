function interimreportfun_PE(elapsedTime,generation,problem,opts,Population,best,Simulation,modelfun,mst,simopts)
% interimreportfun_PE shows fitting.
% 
% [SYNTAX]
% interimreportfun_PE(elapsedTime,generation,problem,opts,Population, ...
%                     best,Simulation,modelfun,mst,simopts)
% 
% [INPUT]
% elapsedTime:  Elaplsed time (sec)
% generation :  Generation
% problem    :  Problem structure
% opts       :  RCGA options. See XXXXXXXXXXX for options.
% Population :  Array of individuals in the current population
% best       :  Structure of the best individua
% Simulation :  Function handle for Simulation_*
% modelfun   :  Function handle for model (odefun or mex)
% mst        :  Experimental data (IQMmeasurement)
% simopts    :  Structure with integrator options. Fields depend on
%               Simulation_*. See 'help Simulation_'.


%% Making outputs
printTransition(elapsedTime,generation,problem,best);
writeTransition(elapsedTime,generation,problem,opts,best);


%% Decoding gene to x
param = problem.decodingfun(best.gene);


%% Converting IQMmeasuremen into structure
mst = struct(mst);


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
[ T, Y ] = feval(Simulation, modelfun, tspan, y0, param, simopts);


%% Making plot
T_exp = mst.time;
n_col = length(mst.data);
n_row = length(mst.data(1).values);
Y_exp = zeros(n_row,n_col);
for i = 1 : n_col
    Y_exp(:,i) = mst.data(i).values;
end
T_sim = T;
Y_sim = Y;
statename = modelfun('states');
plotter(T_sim,Y_sim,T_exp,Y_exp,'Time','AU',statename);
title(sprintf('Generation = %d',generation));
drawnow;
