function Group = requestFitnessCalc(problem,opts,Group)

n_constraint = problem.n_constraint;
par = opts.par;
n_group = length(Group);

f_temp = zeros(n_group,1);
g_temp = zeros(n_group,max(1,n_constraint));
phi_temp = zeros(n_group,1);

if 0 < par
    parfor i = 1 : n_group
        [f_temp(i), g_temp(i,:), phi_temp(i)] = getFitness(problem,Group(i));
    end
else
    for i = 1 : n_group
        [f_temp(i), g_temp(i,:), phi_temp(i)] = getFitness(problem,Group(i));
    end
end

for i = 1 : n_group
    Group(i).f = f_temp(i);
    Group(i).g = g_temp(i,:);
    Group(i).phi = phi_temp(i);
end
