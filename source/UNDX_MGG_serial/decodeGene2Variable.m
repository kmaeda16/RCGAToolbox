function x = decodeGene2Variable(individual, SearchRegion)

x = individual.gene.*(SearchRegion(2,:) - SearchRegion(1,:)) + SearchRegion(1,:);
