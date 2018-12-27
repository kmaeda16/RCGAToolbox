function Population = getInitPopulation(problem,opts)

n_gene = problem.n_gene;
n_population = opts.n_population;
n_constraint = problem.n_constraint;
Pf = opts.Pf;

Population(1,1:n_population) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',0,'phi',0);

for i = 1 : n_population
    for j = 1 : n_gene
        Population(i).gene(j) = rand();
    end
end

Population = requestFitnessCalc(problem,opts,Population);

Population = SRsort(Population,Pf);
