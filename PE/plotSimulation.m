function plotSimulation(Param, best)

n_gene = Param.n_gene;
model = Param.model;
mst = struct(Param.mst{1});
ub = Param.ub;
lb = Param.lb;
x = best.gene .* ( ub - lb ) + lb;
temp = struct(model);
for i = 1 : n_gene
    temp.parameters(i).value = x(i);
end
model = SBmodel(temp);

output = SBsimulate(model,mst.time);

x_sim = output.statevalues;

x_exp = zeros(size(x_sim));
for i = 1 : length(mst.data)
    x_exp(:,i) = mst.data(i).values;
end

st_model = SBstruct(model);
for i = 1 : length(st_model.states)
    statename{i} = st_model.states(i).name;
end

plot(mst.time,x_exp,'o');
legend(statename);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(output.time,x_sim);
hold off;
drawnow;
