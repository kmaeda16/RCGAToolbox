function x = BioPreDynBench_decode(gene)

global lb_global ub_global;

lb = lb_global;
ub = ub_global;

x = gene .* ( ub - lb ) + lb;
