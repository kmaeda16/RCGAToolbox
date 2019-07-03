function x = mydecodingfun(gene)

lb = -1;
ub = 1;

x = gene' * ( ub - lb ) + lb;

param = model_odefun('parametervalues');
x = param .* 10 .^ x;
