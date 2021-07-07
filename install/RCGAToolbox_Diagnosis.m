function [] = RCGAToolbox_Diagnosis
% Diagnosis script.
% This script must be executed under the directory RCGAToolbox/install/.

diary('RCGAToolbox_Diagnosis_Log.txt');


Flag = zeros(1,6);


fprintf('################################################################\n');
fprintf('#                  RCGAToolbox Diagnosis Tool                  #\n');
fprintf('################################################################\n');
fprintf('\n');


%% Checking for System
fprintf('####### Checking for System... #######\n');

archstr = computer('arch');
fprintf('System Architecture: %s\n',archstr);

str =  version;
fprintf('MATLAB Version: %s\n',str);

fprintf('\n\n');


%% Checking for core functions
fprintf('####### Checking for core functions... #######\n');

flg1 = exist('RCGA_Main','file');
if flg1 > 0
    fprintf('RCGA_Main found.\n');
else
    warning('RCGA_Main NOT found.');
end

flg2 = 1;
try
    %========= TEST =========%
    problem.n_gene = 1;
    problem.n_constraint = 0;
    problem.fitnessfun = @(x) x^2; % Fitness Function
    problem.decodingfun = @(gene) 2*gene-1; % Decoding Function
    opts.n_population = 10; % Population Size
    opts.n_children = 10; % Number of Children per Generation
    opts.output_intvl = 500; % Output Interval Generation
    Results = RCGA_UNDXMGG(problem,opts); % UNDX/MGG
    Results = RCGA_REXstarJGG(problem,opts); % REXstar/JGG
    %========================%
catch ME
    warning(ME.message);
    flg2 = 0;
end

if flg2 > 0
    fprintf('Passed Test.\n');
    fprintf('RESULT: Core functions are available.\n');
else
    warning('Failed Test.');
    warning('RCGAToolbox is NOT installed.');
    warning('RESULT: Core functions are NOT available.');
    Flag(1) = 1;
end

fprintf('\n\n');


%% Checking for SBML-related functions
fprintf('####### Checking for SBML-related functions... #######\n');

flg1 = exist('IQMmodel','file');
if flg1 > 0
    fprintf('IQMmodel found.\n');
else
    warning('IQMmodel NOT found.');
end

flg2 = 1;
try
    %========= TEST =========%
    model = RCGAreadConciseODEfile('Model_Example_conciseOdefun.m');
    IQMexportSBML(model,'Model_Example_conciseOdefun_SBML.xml');
    model = IQMmodel('Model_Example_conciseOdefun_SBML.xml');
    %========================%
catch ME
    warning(ME.message);
    flg2 = 0;
end

if flg2 > 0
    fprintf('Passed Test.\n');
    fprintf('RESULT: SBML-related functions are available.\n');
else
    warning('Failed Test.');
    warning('IQM Tools and/or libSBML are NOT installed.');
    warning('RESULT: SBML-related functions are NOT available.');
    Flag(2) = 1;
end

fprintf('\n\n');


%% Checking for fast ODE solver option (fast_flag = 2)
fprintf('####### Checking for fast ODE solver option (fast_flag = 2)... #######\n');
if Flag(2) == 0
    
    flg1 = exist('IQMPsimulate','file');
    if flg1 > 0
        fprintf('IQMPsimulate found.\n');
    else
        warning('IQMPsimulate NOT found.');
    end
    
    flg2 = 1;
    try
        %========= TEST =========%
        model = IQMmodel('Model_Example_SBML.xml');
        simdata = IQMPsimulate(model);
        %========================%
    catch ME
        warning(ME.message);
        flg2 = 0;
    end
    
    if flg2 > 0
        fprintf('Passed Test.\n');
        fprintf('RESULT: Fast ODE solver option (fast_flag = 2) is available.\n');
    else
        warning('Failed Test.');
        warning('IQM Tools is NOT installed.');
        warning('RESULT: Fast ODE solver option (fast_flag = 2) is NOT available.');
        Flag(3) = 1;
    end
    
else
    warning('Test for the fast ODE solver option (fast_flag = 2) could not be run because SBML-related functions are not available.');
    Flag(3) = 1;
end

fprintf('\n\n');


%% Checking for fast ODE solver option (fast_flag = 1)
fprintf('####### Checking for fast ODE solver option (fast_flag = 1)... #######\n');

flg1 = exist('CVode','file');
if flg1 > 0
    fprintf('CVode found.\n');
else
    warning('CVode NOT found.');
end

flg2 = 1;
try
    %========= TEST =========%
    [ T, Y ] = odestb(@(t,y) -y, [0 10], 1);
    %========================%
catch ME
    warning(ME.message);
    flg2 = 0;
end

if flg2 > 0
    fprintf('Passed Test.\n');
    fprintf('RESULT: Fast ODE solver option (fast_flag = 1) is available.\n');
else
    warning('Failed Test.');
    warning('SundialsTB is NOT installed.');
    warning('RESULT: Fast ODE solver option (fast_flag = 1) is NOT available.');
    Flag(4) = 1;
end

fprintf('\n\n');


%% Checking for parallel computation option
fprintf('####### Checking for parallel computation option... #######\n');

flg1 = exist('parfor','file');
if flg1 > 0
    fprintf('parfor found.\n');
else
    warning('parfor NOT found.');
end

flg2 = 1;
try
    %========= TEST =========%
    p = gcp('nocreate');
    if isempty(p)
        parflg = 0;
        p = parpool(2);
    else
        parflg = 1;
    end
    parfor i = 1 : 5
        fprintf('%d\n',i);
    end
    if parflg == 0
        delete(p);
    end
    %========================%
catch ME
    warning(ME.message);
    flg2 = 0;
end

if flg2 > 0
    fprintf('Passed Test.\n');
    fprintf('RESULT: Perallel computation option is available.\n');
else
    warning('Failed Test.');
    warning('Parallel Computation Toolbox is NOT installed.');
    warning('RESULT: Perallel computation option is NOT available.');
    Flag(5) = 1;
end

fprintf('\n\n');


%% Checking for local solver option
fprintf('####### Checking for local solver option... #######\n');

flg1 = exist('fmincon','file');
if flg1 > 0
    fprintf('fmincon found.\n');
else
    warning('fmincon NOT found.');
end

flg2 = 1;
try
    %========= TEST =========%
    fmincon(@(x) x^2, 1, 1, 1);
    %========================%
catch ME
    warning(ME.message);
    flg2 = 0;
end

if flg2 > 0
    fprintf('Passed Test.\n');
    fprintf('RESULT: Local solver option is available.\n');
else
    warning('Failed Test.');
    warning('Optimization Toolbox is NOT installed.');
    warning('RESULT: Local solver option is NOT available.');
    Flag(6) = 1;
end

fprintf('\n\n');


%% Summary

fprintf('############ RCGAToolbox Diagnosis Summary ############\n');

fprintf('Core functions                         : ');
if Flag(1) == 0
    fprintf('Available\n');
else
    fprintf('Not Available\n');
end

fprintf('SBML-related functions                 : ');
if Flag(2) == 0
    fprintf('Available\n');
else
    fprintf('Not Available\n');
end

fprintf('Fast ODE solver option (fast_flag = 2) : ');

if Flag(3) == 0
    fprintf('Available\n');
else
    fprintf('Not Available\n');
end

fprintf('Fast ODE solver option (fast_flag = 1) : ');
if Flag(4) == 0
    fprintf('Available\n');
else
    fprintf('Not Available\n');
end

fprintf('Parallel computation option            : ');
if Flag(5) == 0
    fprintf('Available\n');
else
    fprintf('Not Available\n');
end

fprintf('Local solver option                    : ');
if Flag(6) == 0
    fprintf('Available\n');
else
    fprintf('Not Available\n');
end


diary off;
