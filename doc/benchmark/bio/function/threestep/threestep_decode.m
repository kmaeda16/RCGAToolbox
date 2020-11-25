function x = threestep_decode(gene)

n_gene = length(gene);

lb(1:n_gene) = -12;
ub(1:n_gene) =   6;
lb([3,5,9,11,15,17]) = -1;
ub([3,5,9,11,15,17]) =  1;

x = gene .* ( ub - lb ) + lb;
