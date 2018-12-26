function x = g04_decode(gene)

lb = zeros(1,5);
ub = zeros(1,5);

lb(1) = 78; 
ub(1) = 102;

lb(2) = 33;
ub(2) = 45;

lb(3:5) = 27;
ub(3:5) = 45;

x = gene .* ( ub - lb ) + lb;
