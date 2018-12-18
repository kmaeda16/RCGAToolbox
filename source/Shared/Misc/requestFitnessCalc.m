function Group = requestFitnessCalc(Param,Group)

par = Param.par;
n_constraint = Param.n_constraint;
n_group = length(Group);

f_temp = zeros(n_group,1);
g_temp = zeros(n_group,max(1,n_constraint));
phi_temp = zeros(n_group,1);

if 0 < par
    parfor i = 1 : n_group
        [f_temp(i), g_temp(i,:), phi_temp(i)] = getFitness(Param,Group(i));
%         [f_temp(i), g_temp(i,:), phi_temp(i)] = getFitness2(Param,Group(i));
    end
else
    for i = 1 : n_group
        [f_temp(i), g_temp(i,:), phi_temp(i)] = getFitness(Param,Group(i));
%         [f_temp(i), g_temp(i,:), phi_temp(i)] = getFitness2(Param,Group(i));
    end
end

for i = 1 : n_group
    Group(i).f = f_temp(i);
    Group(i).g = g_temp(i,:);
    Group(i).phi = phi_temp(i);
end
