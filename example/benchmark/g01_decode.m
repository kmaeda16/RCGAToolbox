function x = g01_decode(gene)

lb = zeros(1,13);
ub = zeros(1,13);

lb(1:9) = 0; 
ub(1:9) = 1;

lb(10:12) = 0;
ub(10:12) = 100;

lb(13) = 0;
ub(13) = 1;

x = gene .* ( ub - lb ) + lb;
