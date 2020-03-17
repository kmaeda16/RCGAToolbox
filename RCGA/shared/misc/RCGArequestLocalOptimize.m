function Group = RCGArequestLocalOptimize(problem,opts,Group)


%% Shortening variable names
par = opts.par;
n_group = length(Group);


%% Calculating f, g, and phi
if 0 < par
    parfor i = 1 : n_group
        Group(i) = RCGAlocalOptimize(problem, opts, Group(i));
    end
else
    for i = 1 : n_group
        Group(i) = RCGAlocalOptimize(problem, opts, Group(i));
    end
end
