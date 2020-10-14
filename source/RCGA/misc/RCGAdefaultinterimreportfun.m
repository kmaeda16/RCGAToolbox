function RCGAdefaultinterimreportfun(elapsedTime,generation,problem,opts,Population,best)
% RCGAdefaultinterimreportfun print the fitness etc at every opts.output_intvl.
% 
% [SYNTAX]
% RCGAdefaultinterimreportfun(elapsedTime,generation,problem,opts,Population,best)
% 
% [INPUT]
% elapsedTime :  Elaplsed time (sec)
% generation  :  Generation
% problem     :  Problem structure
% opts        :  RCGA options. See XXXXXXXXXXX for options.
% Population  :  Array of individuals in the current population
% best        :  Structure of the best individual


RCGAprintTransition(elapsedTime,generation,problem,best);
RCGAwriteTransition(elapsedTime,generation,problem,opts,best);
