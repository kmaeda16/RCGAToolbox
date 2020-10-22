function RCGAprintTransition(elapsedTime, generation, problem, best)
% RCGAprintTransition prints elapsed time, generation, fitness, and phi of the
% best individual.
% 
% [SYNTAX]
% RCGAprintTransition(elapsedTime, generation, problem, best)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec).
% generation  :  Generation.
% problem     :  Problem structure:
%                - problem.n_gene: Number of decision variables.
%                - problem.n_constraint: Number of constraint functions. 
%                   For unconstained problems, this must be zero.
%                - problem.fitnessfun: Function handle for a fitness 
%                   function.
%                - problem.decodingfun: Function handle for a decoding 
%                   function.
% best        :  Structure of the the best individual.
%                - best.gene: Decision variable vector.
%                - best.g: Constraint function value vector.
%                - best.f: Fitness function value.
%                - best.phi: Penalty function value.


fprintf('Elapsed Time = %e, Generation = %d, f = %e',elapsedTime,generation,best.f);

if problem.n_constraint > 0
    fprintf(', phi = %e\n',best.phi);
else
    fprintf('\n');
end
