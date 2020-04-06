clearvars;

addpath(genpath('../../../../RCGA'));
addpath('../../../../PE');
addpath('../Common');


IQMmodel('BIOMD0000000295_url')
model = IQMmodel('BIOMD0000000295_url');
IQMcreateODEfile(model,'model_odefun');
tic;
[ T, Y ] = Simulation_odexx(@model_odefun,0:1000);
% [ T, Y ] = Simulation_stb(@model_odefun,0:1000);
toc
plot(T,Y)