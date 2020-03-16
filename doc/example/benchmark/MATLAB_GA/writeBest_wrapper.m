function writeBest_wrapper(elapsedTime, n_gene, n_constraint, fitnessfun, decodingfun, out_best, gene, output)

problem.n_gene = n_gene;
problem.n_constraint = n_constraint;
problem.decodingfun = decodingfun;
opts.out_best = out_best;
opts.n_population = nan;
opts.n_children = nan;

best.gene = gene;
ObjectiveFunction  = @(gene) obj_wrapper(fitnessfun, decodingfun, gene);
best.f = ObjectiveFunction(gene);
if n_constraint > 0
    ConstraintFunction = @(gene) cst_wrapper(fitnessfun, decodingfun, gene);
    best.g = ConstraintFunction(gene);
else
    best.g = 0;
end
best.phi = sum( max(0,best.g) .^2 );

RCGAwriteBest(elapsedTime, output.generations, problem, opts, best)
