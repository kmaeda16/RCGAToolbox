function Population = getInitPopulation(Param)

n_gene = Param.n_gene;
n_population = Param.n_population;
n_constraint = Param.n_constraint;
Pf = Param.Pf;

Population(1,1:n_population) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',0,'phi',0);

for i = 1:n_population
   Population(i).gene = rand(1,n_gene);
end

Population = requestFitnessCalc(Param,Population);

Population = SRsort(Population,Pf);
