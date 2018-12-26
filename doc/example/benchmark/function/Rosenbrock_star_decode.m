function x = Rosenbrock_star_decode(gene)

lb = -2.048;
ub =  2.048;

x = gene * ( ub - lb ) + lb;
