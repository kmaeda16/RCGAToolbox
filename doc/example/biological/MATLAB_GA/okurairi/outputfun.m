function [state,options,optchanged] = outputfun(options,state,flag,problem,opts)
optchanged = false;

elapsedTime = toc;
generation = state.Generation;
x = state.Population(1,:);
neval = state.FunEval;

writeTransition(elapsedTime, generation, problem, opts, x, neval);
