function [ best, Population ] = RCGA_Main(problem, opts, GenerationAlternation)

tic;

interimreportfun = opts.interimreportfun;
finalreportfun = opts.finalreportfun;
n_generation = opts.n_generation;
n_constraint = problem.n_constraint;
output_intvl = opts.output_intvl;
vtr = opts.vtr;
t_limit = opts.t_limit;

flg_printed = 0;


i = 1;
Population = getInitPopulation(problem,opts);
index = findBest(Population);
best = Population(index);

if 0 < output_intvl
    elapsedTime = toc;
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
    flg_printed = 1;
end

elapsedTime = toc;
if ( best.phi == 0 && best.f <= vtr) || elapsedTime >= t_limit
    if 0 < output_intvl && flg_printed == 0
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
    end
    finalreportfun(elapsedTime,problem,opts,Population,best);
    return;
end

while i < n_generation
    i = i + 1;
    flg_printed = 0;
    Population = GenerationAlternation(problem,opts,Population);
    index = findBest(Population);
    if Population(index).phi < best.phi || ( Population(index).phi == best.phi && Population(index).f < best.f )
        best = Population(index);
    end
    if 0 < output_intvl && mod(i,output_intvl) == 0
        elapsedTime = toc;
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
        flg_printed = 1;
    end
    if (best.phi == 0 && best.f <= vtr) || toc >= t_limit
        break;
    end
end

if 0 < output_intvl && flg_printed == 0
    elapsedTime = toc;
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
end

finalreportfun(elapsedTime,problem,opts,Population,best);
