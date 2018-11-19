function [x fitness] = myGA(...
    n_generation, n_population, n_children, n_gene, allowable_error, SearchRegion, FileNameTransition)

population = getInitialPopulation(n_population,n_gene,SearchRegion);
fprintf('%d: Best Fitness = %e\n',1,population(1).fitness);
writeFitnessTransition(1,population(1),SearchRegion,FileNameTransition);
if population(1).fitness <= allowable_error
    x = decodeGene2Variable(population(1),SearchRegion);
    fitness = population(1).fitness;
    return;
end

for i = 2:n_generation
    population = MGGvariant(population,n_children,SearchRegion);
    fprintf('%d: Best Fitness = %e\n',i,population(1).fitness);
    writeFitnessTransition(i,population(1),SearchRegion,FileNameTransition);
    if population(1).fitness <= allowable_error
        x = decodeGene2Variable(population(1),SearchRegion);
        fitness = population(1).fitness;
        return;
    end
end

x = decodeGene2Variable(population(1),SearchRegion);
fitness = population(1).fitness;
