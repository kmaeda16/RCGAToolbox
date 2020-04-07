% Script StoreTransition
% This script is called by RCGA_Main

Results.Transition.time(end+1,1) = elapsedTime;
Results.Transition.neval(end+1,1) = neval;
Results.Transition.generation(end+1,1) = i;
Results.Transition.f(end+1,1) = best.f;
Results.Transition.x(end+1,:) = decodingfun(best.gene);
if n_constraint > 0
    Results.Transition.phi(end+1,1) = best.phi;
    Results.Transition.g(end+1,:) = best.g;
end
if ~isempty(out_report) && ~strcmpi('none',out_report)
    save(out_report,'Results');
end
