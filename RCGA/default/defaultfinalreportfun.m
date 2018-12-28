function defaultfinalreportfun(elapsedTime,generation,problem,opts,Population,best)

writeFinalPopulation(problem,opts,Population);
writeBest(elapsedTime,generation,problem,opts,best);
