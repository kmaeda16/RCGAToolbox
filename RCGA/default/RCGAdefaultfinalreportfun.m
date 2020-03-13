function RCGAdefaultfinalreportfun(elapsedTime,generation,problem,opts,Population,best)
% RCGAdefaultfinalreportfun is executed at the end of RCGA, making two output
% files.
% 
% [SYNTAX]
% RCGAdefaultfinalreportfun(elapsedTime,generation,problem,opts,Population,best)
% 
% [INPUT]
% elapsedTime:  Elaplsed time (sec)
% generation :  Generation
% problem    :  Problem structure
% opts       :  RCGA options. See XXXXXXXXXXX for options.
% Population :  Array of individuals in the current population
% best       :  Structure of the best individual


RCGAwritePopulation(problem,opts,Population);
RCGAwriteBest(elapsedTime,generation,problem,opts,best);
