clearvars;

addpath(genpath('../../../../RCGA'));
addpath('../../../../PE');
addpath('../Common');

simopts.abstol = 1e-6;
simopts.reltol = 1e-6;

% IQMmodel('BIOMD0000000295_url')
% model = IQMmodel('BIOMD0000000295_url');
% IQMcreateODEfile(model,'model_odefun');
tic;
% [ T, Y ] = RCGAsimulate('BIOMD0000000295_url',0:1000);
% [ T, Y ] = RCGAsimulate('BIOMD0000000295_url',0:1000,[],[],2);
% [ T, Y ] = RCGAsimulate('Akman2008_Circadian_Clock_Model1_odefun.m',0:1000,[],[],0,simopts);
% [ T, Y ] = RCGAsimulate('Akman2008_Circadian_Clock_Model1_odefun.m',0:1000,[],[],1,simopts);
% [ T, Y ] = RCGAsimulate('Akman2008_Circadian_Clock_Model1_mex.c',0:1000,[],[],1,simopts);
% [ T, Y ] = RCGAsimulate('Akman2008_Circadian_Clock_Model1_mex.mexw64',0:1000,[],[],1,simopts);
% [ T, Y ] = RCGAsimulate('BIOMD0000000217_url');
% [ T, Y ] = RCGAsimulate('BIOMD0000000217_url',[],[],[],1);
% [ T, Y ] = RCGAsimulate('Maeda2019_AmmoniumTransportAssimilation_odefun.m',0:100,[],[],0);
% [ T, Y ] = RCGAsimulate('Maeda2019_AmmoniumTransportAssimilation_odefun.m',0:100,[],[],1);
% [ T, Y ] = RCGAsimulate('Maeda2019_AmmoniumTransportAssimilation_mex.mexw64',0:100,[],[],2);

[ T, Y ] = RCGAsimulate('EcoliCentralCarbonMetabolism_odefun',0:100,[],[],1);

toc
plot(T,Y)