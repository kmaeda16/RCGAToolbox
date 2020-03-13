% Script StoreBestAndFinalPopulation
% This script is called by RCGA_Main

Results.Best.time = elapsedTime;
Results.Best.neval = n_population + ( i - 1 ) * n_children;
Results.Best.generation = i;
Results.Best.f = best.f;
Results.Best.x = decodingfun(best.gene);
if n_constraint > 0
    Results.Best.phi = best.phi;
    Results.Best.g = best.g;
end
for j = 1 : n_population
    Results.FinalPopulation.f(j,1) = Population(j).f;
    Results.FinalPopulation.x(j,:) = decodingfun(Population(j).gene);
    if n_constraint > 0
        Results.FinalPopulation.phi(j,1) = Population(j).phi;
        Results.FinalPopulation.g(j,:) = Population(j).g;
    end
end
if best.phi == 0 && best.f <= vtr
    Results.end_crit = 0;
    fprintf('Value to reach achieved.\n');
end
if elapsedTime >= t_limit
    Results.end_crit = 1;
    fprintf('Maximum allowed CPU time achieved.\n');
end
if i >= n_generation
    Results.end_crit = 2;
    fprintf('Maximal number of generations achieved.\n');
end
if ~strcmpi('none',out_report)
    save(out_report,'Results');
end
