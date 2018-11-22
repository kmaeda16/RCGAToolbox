function x = g13_decode(gene)

lb = zeros(1,5);
ub = zeros(1,5);

lb(1:2) = -2.3;
ub(1:2) =  2.3;

lb(3:5) = -3.2;
ub(3:5) =  3.2;

x = gene .* ( ub - lb ) + lb;
