function x = Ellipsoid_decode(gene)

lb = -5.12;
ub =  5.12;

x = gene * ( ub - lb ) + lb;
