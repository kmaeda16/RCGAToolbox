clearvars;

N_TRIAL = 11;

addpath(genpath('../../../../../RCGA'));
addpath('../../../../../PE');
addpath('../../Common');

simopts.abstol = 1e-6;
simopts.reltol = 1e-6;

SBMLfilename = 'BIOMD0000000016_url.xml';


%%
display('Making SBML file ...');
SBMLfilename = replaceWords(SBMLfilename);
fprintf('%s.xml created.\n',SBMLfilename);

%%
display('Making M file ...');
odefilename = RCGAmakeODEmodel(SBMLfilename);
fprintf('%s.m created.\n',odefilename);

%% Making SBMLfile
display('Making MEX file ...');
mexfilename = RCGAmakeMEXmodel(SBMLfilename);
fprintf('%s.c created.\n',mexfilename);
fprintf('%s.h created.\n',mexfilename);
fprintf('%s.%s created.\n',mexfilename,mexext);

%%
simtime = zeros(N_TRIAL,3);


%%
for i = 1 : N_TRIAL
    tic;
    [ T, Y ] = RCGAsimulate(odefilename,0:100,[],[],0,simopts);
    simtime(i,1) = toc;
    tic;
    [ T, Y ] = RCGAsimulate(odefilename,0:100,[],[],1,simopts);
    simtime(i,2) = toc;
    tic;
    [ T, Y ] = RCGAsimulate(mexfilename,0:100,[],[],2,simopts);
    simtime(i,3) = toc;
end

simtime(1,:) = [];
avg = mean(simtime);
stdev = std(simtime);

fprintf('                 \tMean\tSD\n');
fprintf('ODE15s           \t%.2e\t%.2e\n',avg(1),stdev(1));
fprintf('CVODE(SundialsTB)\t%.2e\t%.2e\n',avg(2),stdev(2));
fprintf('CVODE(IQMTools)  \t%.2e\t%.2e\n',avg(3),stdev(3));

