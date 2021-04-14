function Results = RCGA_UNDXMGG_PE(model,decodingfun,mst,varargin)
% RCGA_UNDXMGG_PE estimates parameters included in model by fitting it to 
% experimental data mst, using UNDX/MGG.
% 
% [SYNTAX]
% Results = RCGA_UNDXMGG_PE(model, decodingfun, mst)
% Results = RCGA_UNDXMGG_PE(model, decodingfun, mst, opts)
% Results = RCGA_UNDXMGG_PE(model, decodingfun, mst, simopts, opts)
% Results = RCGA_UNDXMGG_PE(model, decodingfun, mst, fast_flag, simopts,
%                           opts)
% Results = RCGA_UNDXMGG_PE(model, decodingfun, mst, n_constraint, ...
%                           fitnessfun, fast_flag, simopts, opts)
% 
% [INPUT]
% model        :  An IQMmodel object, the name of SBML file, the function 
%                 handle for an ODE function (IQM Tools format), the 
%                 function handle for a MEXed model, or the C source code.
% decodingfun  :  Function handle for a decoding function.
% mst          :  Experimental data (IQMmeasurement) or file name (*.xls).
% n_constraint :  Number of constraint functions.
% fitnessfun   :  Function handle for a fitness function.
% fast_flag    :  ODE solver flag:
%                 - fast_flag = 0: ODEXX by MATLAB built-ins.
%                 - fast_flag = 1: CVODE by SundialsTB.
%                 - fast_flag = 2: CVODE by IQM Tools.
% simopts      :  Solver option structure. The fields depend on fast_flag. 
%                 For fast_flag = 0, 1, and 2, see 
%                 'help RCGAsimulateODE', 'help RCGAsimulateSTB', 
%                 'help RCGAsimulateMEX', respectively.
% opts         :  Option structure:
%                 - opts.n_population: Population size.
%                 - opts.n_children: Number of children.
%                 - opts.n_parent: Number of parents.
%                 - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%                 - opts.selection_type: Selection type for REXstar/JGG 
%                    (0 or 1).
%                 - opts.Pf: Probability that only the objective function f
%                    is used in comparisons of individuals in the 
%                    stochastic ranking.
%                 - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                    local optimizer is used.
%                 - opts.localopts: Options for the local optimizer.
%                 - opts.maxgen: Maximum number of generations.
%                 - opts.maxtime: Maximum time (sec).
%                 - opts.maxeval: Maximum number of fitnessfun evaluations.
%                 - opts.vtr: Value to be reached.
%                 - opts.n_par: Number of workers in parallel computation.
%                 - opts.output_intvl: Interval generation for updating the 
%                    transition file and the report file.
%                 - opts.out_transition: Name of an output file called the 
%                    transition file.
%                 - opts.out_best: Name of an output file called the best 
%                    individual file.
%                 - opts.out_population: Name of an output file called the 
%                    final population file.
%                 - opts.out_report: Name of an output file called the 
%                    report file.
%                 - opts.interimreportfun: Function handle for the interim 
%                    report function.
%                 - opts.finalreportfun: Function handle for the final 
%                    report function.
% 
% [OUTPUT]
% Results      :  Results structure:
%                 - Results.Transition: Information on the fitness 
%                    transition.
%                 - Results.Best: Information on the best individual.
%                 - Results.FinalPopulation: Information on the final 
%                    population.
%                 - Results.end_crit: Exit flag: Success (0), maxgen 
%                    reached (1), maxtime reached (2), maxeval reached (3).


%% Handling input arguments
switch nargin
    case 3
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = 0;
        simopts = struct;
        opts = struct;
    case 4
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = 0;
        simopts = struct;
        opts = varargin{1};
    case 5
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = 0;
        simopts = varargin{1};
        opts = varargin{2};
    case 6
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = varargin{1};
        simopts = varargin{2};
        opts = varargin{3};
    case 8
        n_constraint = varargin{1};
        fitnessfun = varargin{2};
        fast_flag = varargin{3};
        simopts = varargin{4};
        opts = varargin{5};
    otherwise
        error('Incorrect number of input arguments.');
end

if isempty(n_constraint)
    n_constraint = 0;
end
if isempty(fitnessfun)
    fitnessfun = @RCGAssr;
end
if isempty(fast_flag)
    fast_flag = 0;
end
if isempty(simopts)
    simopts = struct;
end
if isempty(opts)
    opts = struct;
end


%% Run parameter estimation
Results = RCGA_PE(model,decodingfun,mst,n_constraint,fitnessfun,fast_flag,simopts,opts,@RCGA_UNDXMGG);
