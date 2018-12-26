function x = Ackley_decode(gene)

lb = -32.768;
ub =  32.768;

x = gene * ( ub - lb ) + lb;
