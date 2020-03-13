function RCGAprintTransition(elapsedTime, generation, problem, best)
% RCGAprintTransition prints elapsed time, generation, fitness, and phi of the
% best individual.
% 
% [SYNTAX]
% RCGAprintTransition(elapsedTime, generation, problem, best)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec)
% generation  :  Generation
% problem     :  Problem structure
% best        :  Structure of the the best individual


fprintf('Elapsed Time = %e, Generation = %d, f = %e',elapsedTime,generation,best.f);

if problem.n_constraint > 0
    fprintf(', phi = %e\n',best.phi);
else
    fprintf('\n');
end
