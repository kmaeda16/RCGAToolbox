function x = g10_decode(gene)

lb = zeros(1,8);
ub = zeros(1,8);

lb(1) = 100; 
ub(1) = 10000;

lb(2:3) = 1000;
ub(2:3) = 10000;

lb(4:8) = 10;
ub(4:8) = 1000;

x = gene .* ( ub - lb ) + lb;
