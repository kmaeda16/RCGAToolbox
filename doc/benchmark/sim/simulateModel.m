function simulateModel(SBMLfilename)
% simulateModel repeats simulations with different solvers and print
% computational time.
% 
% [SYNTAX]
% simulateModel(SBMLfilename)
% 
% [INPUT]
% SBMLfilename : Name of SBML file.

N_TRIAL = 100;

simopts = struct;
AbsTol = 1e-6;
RelTol = 1e-6;


%% Print current directory
pwd


%% Remove reserved words
display('Removing reserved words in the SBML file ...');
SBMLfilename = RCGAreplaceWords(SBMLfilename);
fprintf('%s created.\n',SBMLfilename);


%% Make an ODE file (IQM Tools format)
display('Making an ODE file (IQM Tools format) ...');
odefilename = RCGAcreateODEfile(SBMLfilename);
fprintf('%s.m created.\n',odefilename);


%% Make a MEX model
display('Making a MEX model ...');
mexfilename = RCGAmakeMEXmodel(SBMLfilename);
fprintf('%s.c created.\n',mexfilename);
fprintf('%s.h created.\n',mexfilename);
fprintf('%s.%s created.\n',mexfilename,mexext);


%% Simulation
simtime = zeros(N_TRIAL,3);

for i = 1 : N_TRIAL
    
    % ODE15s(MATLABbuilt-in)
    simopts.AbsTol = AbsTol;
    simopts.RelTol = RelTol;
    simopts.Method = 'ode15s';
    simopts.BDF = 'on';
    tic;
    [ T1, Y1 ] = RCGAsimulate(odefilename,0:100,[],[],0,simopts);
    simtime(i,1) = toc;
    
    % CVODE(SundialsTB)
    simopts.AbsTol = AbsTol;
    simopts.RelTol = RelTol;
    tic;
    [ T2, Y2 ] = RCGAsimulate(odefilename,0:100,[],[],1,simopts);
    simtime(i,2) = toc;
    
    % CVODE(IQMTools)
    simopts.abstol = AbsTol;
    simopts.reltol = RelTol;
    tic;
    [ T3, Y3 ] = RCGAsimulate(mexfilename,0:100,[],[],2,simopts);
    simtime(i,3) = toc;
end

%% Check if models were simulated correctly
fprintf('Diff1: %e\n',max(max(abs((Y2-Y1)./Y1))));
fprintf('Diff2: %e\n',max(max(abs((Y3-Y1)./Y1))));


%% Print results
avg = mean(simtime);
stdev = std(simtime);

fprintf('                      \tMean\tSD\n');
fprintf('ODE15s(MATLABbuilt-in)\t%.2e\t%.2e\n',avg(1),stdev(1));
fprintf('CVODE(SundialsTB)     \t%.2e\t%.2e\n',avg(2),stdev(2));
fprintf('CVODE(IQMTools)       \t%.2e\t%.2e\n',avg(3),stdev(3));
