function x = g10_decode(gene)

lb = zeros(1,8);
ub = zeros(1,8);

lb(1) = 2; 
ub(1) = 4;

lb(2:3) = 3;
ub(2:3) = 4;

lb(4:8) = 1;
ub(4:8) = 3;

x = 10 .^ ( gene .* ( ub - lb ) + lb );
