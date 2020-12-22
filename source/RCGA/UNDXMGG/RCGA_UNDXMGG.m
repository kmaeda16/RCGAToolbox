function Results = RCGA_UNDXMGG(problem, varargin)
% RCGA_UNDXMGG is the main function of UNDX/MGG.
% 
% [SYNTAX]
% Results = RCGA_UNDXMGG(problem)
% Results = RCGA_UNDXMGG(problem, opts)
% 
% [INPUT]
% problem :  Problem structure:
%            - problem.n_gene: Number of decision variables.
%            - problem.n_constraint: Number of constraint functions. For 
%               unconstained problems, this must be zero.
%            - problem.fitnessfun: Function handle for a fitness function.
%            - problem.decodingfun: Function handle for a decoding 
%               function.
% opts    :  Option structure:
%            - opts.n_population: Population size.
%            - opts.n_children: Number of children.
%            - opts.n_parent: Number of parents.
%            - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%            - opts.selection_type: Selection type for REXstar/JGG 
%               (0 or 1).
%            - opts.Pf: Probability that only the objective function f is 
%               used in comparisons of individuals in the stochastic 
%               ranking.
%            - opts.local: Local optimizer (0 or 1). If it is 1, the local 
%               optimizer is used.
%            - opts.localopts: Options for the local optimizer.
%            - opts.maxgen: Maximum number of generations.
%            - opts.maxtime: Maximum time (sec).
%            - opts.maxeval: Maximum number of fitnessfun evaluations.
%            - opts.vtr: Value to be reached.
%            - opts.n_par: Number of workers in parallel computation.
%            - opts.output_intvl: Interval generation for updating the 
%               transition file and the report file.
%            - opts.out_transition: Name of an output file called the 
%               transition file.
%            - opts.out_best: Name of an output file called the best 
%               individual file.
%            - opts.out_population: Name of an output file called the 
%               final population file.
%            - opts.out_report: Name of an output file called the report 
%               file.
%            - opts.interimreportfun: Function handle for the interim 
%               report function.
%            - opts.finalreportfun: Function handle for the final report 
%               function.
% 
% [OUTPUT]
% Results :  Results structure:
%            - Results.Transition: Information on the fitness transition.
%            - Results.Best: Information on the best individual.
%            - Results.FinalPopulation: Information on the final 
%               population.
%            - Results.end_crit: Exit flag: Success (0), maxgen reached 
%               (1), maxtime reached (2), maxeval reached (3).


%% Handling inputs
switch nargin
    case 1
        opts = struct;
    case 2
        opts = varargin{1};
    otherwise
        error('Incorrect number of input arguments');
end


%% Checking inputs
[problem, opts] = RCGAcheckInputs(problem,opts,str2func(mfilename));


%% Printing welcome messages
RCGAprintWelcomeMessage(problem,opts,str2func(mfilename));


%% Executing UNDX/MGG
Results = RCGA_Main(problem,opts,@RCGA_MGG);
