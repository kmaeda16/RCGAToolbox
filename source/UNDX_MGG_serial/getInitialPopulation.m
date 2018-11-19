function population = getInitialPopulation(n_population, n_gene, SearchRegion)

for i = 1:n_population
   population(i).gene = rand(1,n_gene);
   x = decodeGene2Variable(population(i),SearchRegion);
   population(i).fitness = getFitness(x);
end
population = mysort(population);
