function x = Schwefel_decode(gene)

lb = -512;
ub =  512;

x = gene * ( ub - lb ) + lb;
