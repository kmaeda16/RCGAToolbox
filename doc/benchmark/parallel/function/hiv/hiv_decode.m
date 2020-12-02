function x = hiv_decode(gene)

lb(1 : 5) = 1e-6;  ub( 1: 5) = 1e+6;
lb(6 :10) = 10;    ub( 6:10) = 40;
lb(11:15) = 0.002; ub(11:15) = 0.006;
lb(16:20) = -0.1;  ub(16:20) = 0.1;

x = gene .* ( ub - lb ) + lb;
