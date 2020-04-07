function CreateExecutableScript(app)

filename = app.ExecutableScriptFileName.Value;

%% Opening Output File
out = fopen(filename,'w');
if out == -1
    warning('cannot open %s!\n',filename);
    return;
end

%% Print Information
fprintf(out,'%% This script was created by RCGAToolbox Mission Control\n');
fprintf(out,'%% %s\n\n',date);

%% Problem Settings
fprintf(out,'%% ========= Problem Settings ========= %%\n');

fprintf(out,'problem.n_gene = %d; %% # Variables\n',app.N_Variables.Value);
fprintf(out,'problem.n_constraint = %d; %% # Constraints\n',app.N_Constraints.Value);

[fitnesspath, fitnessfile, ~] = fileparts(app.FitnessFunction.Value);
fprintf(out,'fitnesspath = ''%s''; %% Path to Fitness Function File\n',fitnesspath);
fprintf(out,'addpath(fitnesspath);\n');
fprintf(out,'problem.fitnessfun = @%s; %% Fitness Function\n',fitnessfile);

[decodingpath, decodingfile, ~] = fileparts(app.DecodingFunction.Value);
fprintf(out,'decodingpath = ''%s''; %% Path to Decoding Function File\n',decodingpath);
fprintf(out,'addpath(decodingpath);\n');
fprintf(out,'problem.decodingfun = @%s; %% Decoding Function\n',decodingfile);

fprintf(out,'\n');

%% Option Settings
fprintf(out,'%% ========= Option Settings ========== %%\n');

fprintf(out,'opts.n_population = %d; %% Population Size\n',app.PopulationSize.Value);
fprintf(out,'opts.n_children = %d; %% # Children per Generation\n',app.ChildrenSize.Value);
if strcmp(app.AlgorithmSwitch.Value,'REXstar/JGG')
    fprintf(out,'opts.n_parent = problem.n_gene + 1; %% # Parents used for REXstar\n');
    fprintf(out,'opts.t_rexstar = 6.0; %% Step-size parameter for REXstar\n');
    fprintf(out,'opts.selection_type = 0; %% Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)\n');
end
fprintf(out,'opts.Pf = %e; %% Pf\n',app.Pf.Value);
fprintf(out,'opts.n_generation = %d; %% Max # Generations\n',app.Max_N_Generation.Value);
fprintf(out,'opts.maxtime = 60 * %e; %% Max Time (sec)\n',app.MaxTime.Value);
fprintf(out,'opts.maxeval = inf; %% Max # fitnessfun Evaluations\n');
fprintf(out,'opts.vtr = %e; %% Value To Be Reached\n',app.ValueToBeReached.Value);
fprintf(out,'opts.output_intvl = %d; %% Output Interval Generation\n',app.OutputIntervalGeneration.Value);
fprintf(out,'opts.out_transition = ''%s''; %% Transition File Name\n',app.TransitionFileName.Value);
fprintf(out,'opts.out_best = ''%s''; %% Best Individual File Name\n',app.BestIndividualFileName.Value);
fprintf(out,'opts.out_population = ''%s''; %% Final Population File Name\n',app.FinalPopulationFileName.Value);
fprintf(out,'opts.out_report = ''%s''; %% Report File Name\n',app.ReportFileName.Value);
fprintf(out,'opts.n_par = %d; %% # Workders\n',app.N_Workers.Value);
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
if strcmp(app.AlgorithmSwitch.Value,'UNDX/MGG')
    fprintf(out,'RCGA_UNDXMGG(problem,opts);\n');
else
    fprintf(out,'RCGA_REXstarJGG(problem,opts);\n');
end
fprintf(out,'\n');

%% Removing Path
fprintf(out,'%% ========== Removing Path =========== %%\n');
fprintf(out,'rmpath(fitnesspath);\n');
fprintf(out,'if ~strcmp(fitnesspath,decodingpath)\n');
fprintf(out,'    rmpath(decodingpath);\n');
fprintf(out,'end\n');

%% Print Path and File name
fprintf('Executable script "%s" was created in "%s"\n',filename,pwd);

%% Closing Output File
fclose(out);
