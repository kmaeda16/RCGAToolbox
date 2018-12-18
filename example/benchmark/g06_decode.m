function x = g06_decode(gene)

lb = zeros(1,2);
ub = zeros(1,2);

lb(1) = 13; 
ub(1) = 100;

lb(2) = 0;
ub(2) = 100;

x = gene .* ( ub - lb ) + lb;
