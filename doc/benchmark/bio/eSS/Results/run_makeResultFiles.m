
PROBLEM = {'hiv','threestep'};
n_repeat = 5;

for i = 1 : length(PROBLEM)
    problem = char(PROBLEM(i));
    for j = 1 : n_repeat
        infilename = sprintf('../result_%s_%d/ess_report.mat',problem,j);
        outfilename = sprintf('eSS_%s_transition_%d.dat',problem,j);
        makeResultFiles(infilename,outfilename);
    end
end
