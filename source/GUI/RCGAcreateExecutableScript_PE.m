function RCGAcreateExecutableScript_PE(app)
% RCGAcreateExecutableScript_PE creates an executable script for parameter
% estimation using real-coded genetic algorithms (RCGAs). This function
% must be called by the GUI RCGA_MissionControl_PE.

filename = app.ExecutableScriptFileName.Value;

%% Opening Output File
out = fopen(filename,'w');
if out == -1
    warning('cannot open %s!\n',filename);
    return;
end

%% Print Information
fprintf(out,'%% This script was created by RCGAToolbox Mission Control PE\n');
fprintf(out,'%% Generated: %s\n',date);

%% Clear variables
fprintf(out,'clearvars;\n');
fprintf(out,'\n');

%% Problem Settings
fprintf(out,'%% ========= Problem Settings ========= %%\n');

[modelpath, modelfile, modelext] = fileparts(app.Model.Value);
fprintf(out,'modelpath = ''%s''; %% Path to Model File\n',modelpath);
fprintf(out,'addpath(modelpath);\n');
addpath(modelpath)
if exist(modelfile)
    fprintf(out,'modelfile = @%s; %% Model File\n',modelfile);
else
    modelfile = strcat(modelpath,'/',modelfile,modelext);
    fprintf(out,'modelfile = ''%s''; %% Model File\n',modelfile);
end


[decodingpath, decodingfile, ~] = fileparts(app.DecodingFunction.Value);
fprintf(out,'decodingpath = ''%s''; %% Path to Decoding Function File\n',decodingpath);
fprintf(out,'addpath(decodingpath);\n');
fprintf(out,'decodingfun = @%s; %% Decoding Function\n',decodingfile);

measurement = app.Measurement.Value;
fprintf(out,'measurement = ''%s''; %% Measurement File\n',measurement);

fprintf(out,'\n');

%% Option Settings
fprintf(out,'%% ========= Option Settings ========== %%\n');

fprintf(out,'opts.n_population = %d; %% Population Size\n',app.PopulationSize.Value);
fprintf(out,'opts.n_children = %d; %% # Children per Generation\n',app.ChildrenSize.Value);
if strcmp(app.AlgorithmSwitch.Value,'REXstar/JGG')
    fprintf(out,'%% opts.n_parent = n_gene + 1; %% # Parents Used for REXstar\n');
    fprintf(out,'opts.t_rexstar = 6.0; %% Step-size Parameter for REXstar\n');
    fprintf(out,'opts.selection_type = 0; %% Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)\n');
end
fprintf(out,'opts.maxgen = %d; %% Max # Generations\n',app.MaxGenerations.Value);
fprintf(out,'opts.maxtime = 60 * %e; %% Max Time (sec)\n',app.MaxTime.Value);
fprintf(out,'opts.maxeval = inf; %% Max # fitnessfun Evaluations\n');
fprintf(out,'opts.vtr = %e; %% Value To Be Reached\n',app.ValueToBeReached.Value);
fprintf(out,'opts.output_intvl = %d; %% Output Interval Generation\n',app.OutputIntervalGeneration.Value);
fprintf(out,'opts.out_transition = ''%s''; %% Transition File Name\n',app.TransitionFileName.Value);
fprintf(out,'opts.out_best = ''%s''; %% Best Individual File Name\n',app.BestIndividualFileName.Value);
fprintf(out,'opts.out_population = ''%s''; %% Final Population File Name\n',app.FinalPopulationFileName.Value);
fprintf(out,'opts.out_report = ''%s''; %% Report File Name\n',app.ReportFileName.Value);
fprintf(out,'opts.n_par = %d; %% # Workers for Parallel Computation\n',app.N_Workers.Value);
fprintf(out,'fast_flag = %d; %% fast_flag\n',app.FastFlag.Value);
if strcmp(app.LocalSwitch.Value,'Off')
    fprintf(out,'opts.local = 0; %% Local Optimizer\n');
else
    fprintf(out,'opts.local = 1; %% Local Optimizer\n');
end
fprintf(out,'\n');

%% Setting Random Seed
fprintf(out,'%% ======= Setting Random Seed ======== %%\n');
fprintf(out,'rng(%d); %% Random Seed\n',app.Seed.Value);
fprintf(out,'\n');

%% Executing RCGA
fprintf(out,'%% ========== Executing RCGA ========== %%\n');
fprintf(out,'clear RCGAssr;\n');
if strcmp(app.AlgorithmSwitch.Value,'UNDX/MGG')
    fprintf(out,'Results = RCGA_UNDXMGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts);\n');
else
    fprintf(out,'Results = RCGA_REXstarJGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts);\n');
end
fprintf(out,'\n');

%% Removing Path
fprintf(out,'%% ========== Removing Path =========== %%\n');
fprintf(out,'rmpath(modelpath);\n');
fprintf(out,'if ~strcmp(modelpath,decodingpath)\n');
fprintf(out,'    rmpath(decodingpath);\n');
fprintf(out,'end\n');

%% Print Path and File name
fprintf('Executable script "%s" was created in "%s"\n',filename,pwd);

%% Closing Output File
fclose(out);
