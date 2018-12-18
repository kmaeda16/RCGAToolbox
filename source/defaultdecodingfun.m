function x = defaultdecodingfun(gene)

lb = -1e+3;
ub =  1e+3;

x = gene * ( ub - lb ) + lb;
