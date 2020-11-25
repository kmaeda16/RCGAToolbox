
n_repeat = 5;


%% Start Calculation
for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
        dirname = sprintf('result_%s_%d',name,i);
        mkdir(dirname);
        cd(dirname);
        addpath('..');
%         Biological_eSS(name,i); % For normal calculation
        batch(@Biological_eSS,0,{name,i}); % For batch calculation
        rmpath('..');
        cd('..');
    end
end


%% Make Result files
% Make the directory Results and execute below in it.
PROBLEM = {'hiv','threestep'};
n_repeat = 5;

for i = 1 : length(PROBLEM)
    name = char(PROBLEM(i));
    for j = 1 : n_repeat
        infilename = sprintf('../result_%s_%d/ess_report.mat',name,j);
        outfilename = sprintf('eSS_%s_transition_%d.dat',name,j);
        makeResultFiles(infilename,outfilename);
    end
end
