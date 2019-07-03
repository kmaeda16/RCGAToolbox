function defaultinterimreportfun(elapsedTime,generation,problem,opts,Population,best)
% defaultinterimreportfun print the fitness etc at every opts.output_intvl.
% 
% [SYNTAX]
% defaultinterimreportfun(elapsedTime,generation,problem,opts,Population,best)
% 
% [INPUT]
% elapsedTime:  Elaplsed time (sec)
% generation :  Generation
% problem    :  Problem structure
% opts       :  RCGA options. See XXXXXXXXXXX for options.
% Population :  Array of individuals in the current population
% best       :  Structure of the best individual


printTransition(elapsedTime,generation,problem,best);
writeTransition(elapsedTime,generation,problem,opts,best);
