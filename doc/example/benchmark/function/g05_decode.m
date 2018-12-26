function x = g05_decode(gene)

lb = zeros(1,4);
ub = zeros(1,4);

lb(1:2) = 0; 
ub(1:2) = 1200;

lb(3:4) = -0.55;
ub(3:4) =  0.55;

x = gene .* ( ub - lb ) + lb;
