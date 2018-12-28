Results.Transition.time = elapsedTime;
Results.Transition.neval = n_population;
Results.Transition.generation = i;
Results.Transition.f = best.f;
Results.Transition.x = decodingfun(best.gene);
if n_constraint > 0
    Results.Transition.phi = best.phi;
    Results.Transition.g = best.g;
end
if ~strcmpi('none',out_report)
    save(out_report,'Results');
end