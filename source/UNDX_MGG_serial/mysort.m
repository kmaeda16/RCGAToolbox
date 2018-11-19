function population = mysort(population)

n_population = length(population);
n_gene = length(population(1).gene);

temp = zeros(n_population,n_gene+1);
for i = 1:n_population
    temp(i,:) = [population(i).gene population(i).fitness];
end
temp = sortrows(temp,n_gene+1);
for i = 1:n_population
    for j = 1:n_gene
        population(i).gene(j) = temp(i,j);
    end
    population(i).fitness = temp(i,n_gene+1);
end
