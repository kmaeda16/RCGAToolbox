function x = k_tablet_decode(gene)

lb = -5.12;
ub =  5.12;

x = gene * ( ub - lb ) + lb;
