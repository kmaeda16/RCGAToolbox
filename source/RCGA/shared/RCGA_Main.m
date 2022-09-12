function Results = RCGA_Main(problem, opts, GenerationAlternation)
% RCGA_Main is the main function of real-coded genetic algorithms.
% 
% [SYNTAX]
% Results = RCGA_Main(problem, opts, GenerationAlternation)
% 
% [INPUT]
% problem               :  Problem structure:
%                          - problem.n_gene: Number of decision variables.
%                          - problem.n_constraint: Number of constraint 
%                             functions. For unconstained problems, this 
%                             must be zero.
%                          - problem.fitnessfun: Function handle for 
%                             a fitness function.
%                          - problem.decodingfun: Function handle for 
%                             a decoding function.
% opts                  :  Option structure:
%                          - opts.n_population: Population size.
%                          - opts.n_children: Number of children.
%                          - opts.n_parent: Number of parents.
%                          - opts.t_rexstar: Step-size parameter for
%                             REXstar/JGG.
%                          - opts.selection_type: Selection type for
%                             REXstar/JGG (0 or 1).
%                          - opts.Pf: Probability that only the objective 
%                             function f is used in comparisons of 
%                             individuals in the stochastic ranking.
%                          - opts.local: Local optimizer (0 or 1). If it is
%                             1, the local optimizer is used.
%                          - opts.localopts: Options for the local 
%                             optimizer.
%                          - opts.maxgen: Maximum number of generations.
%                          - opts.maxtime: Maximum time (sec).
%                          - opts.maxeval: Maximum number of fitnessfun 
%                             evaluations.
%                          - opts.maxstallgen: Maximum number of stall 
%                             generations for early stopping.
%                          - opts.vtr: Value to be reached.
%                          - opts.n_par: Number of workers in parallel 
%                             computation.
%                          - opts.initial_population: n x n_gene matrix in
%                             which each row represents an individual. Note
%                             that each gene should be 0 ~ 1. The first
%                             n_population individuals of the designated 
%                             initial population are used, and others are 
%                             ignored. If n < n_population, 
%                             n_population - n individuals are randomly 
%                             generated. 
%                          - opts.output_intvl: Interval generation for 
%                             updating the transition file and the report 
%                             file.
%                          - opts.out_transition: Name of an output file 
%                             called the transition file.
%                          - opts.out_best: Name of an output file called 
%                             the best individual file.
%                          - opts.out_population: Name of an output file 
%                             called the final population file.
%                          - opts.out_report: Name of an output file called
%                             the report file.
%                          - opts.interimreportfun: Function handle for the
%                             interim report function.
%                          - opts.finalreportfun: Function handle for the 
%                             final report function.
% GenerationAlternation :  Function handle for a generation alternation
%                          function.
% 
% [OUTPUT]
% Results               :  Results structure:
%                          - Results.Transition: Information on the fitness
%                             transition.
%                          - Results.Best: Information on the best
%                             individual.
%                          - Results.FinalPopulation: Information on the 
%                             final population.
%                          - Results.end_crit: Exit flag: Success (0),
%                             maxgen reached (1), maxtime reached (2), 
%                             maxeval reached (3), maxstallgen (4).


%% Getting the time RCGA starts
tic;


%% Shortening variable names
decodingfun = problem.decodingfun;
n_constraint = problem.n_constraint;
interimreportfun = opts.interimreportfun;
finalreportfun = opts.finalreportfun;
maxgen = opts.maxgen;
output_intvl = opts.output_intvl;
n_population = opts.n_population;
out_report = opts.out_report;
n_children = opts.n_children;
vtr = opts.vtr;
maxtime = opts.maxtime;
maxeval = opts.maxeval;
maxstallgen = opts.maxstallgen;
n_par = opts.n_par;
local = opts.local;


global RCGA_LOCALNEVAL;
RCGA_LOCALNEVAL = 0;


%% Making parallel pool
if 1 < n_par
    p = gcp('nocreate');
    if isempty(p)
        parpool(n_par);
    elseif ~isempty(p)
        if p.NumWorkers ~= n_par
            delete(p);
            parpool(n_par);
        end
    end
end


%% Setting flg_printed
flg_printed = 0; % flg_printed == 1 indicates output is made


%% Stall generations
stallgenerations = 0;


%% First generation
i = 1;
Population = RCGAgetInitPopulation(problem,opts);
index = RCGAfindBest(Population);
if local > 0
    [Population(index), localneval] = RCGAlocalOptimize(problem, opts, Population(index));
    RCGA_LOCALNEVAL = RCGA_LOCALNEVAL + localneval;
end
best = Population(index);

elapsedTime = toc;
neval = n_population + ( i - 1 ) * n_children + RCGA_LOCALNEVAL;

if 0 < output_intvl
    elapsedTime = toc;
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptInitTransition; % RCGA/shared/misc/ScriptInitTransition
    ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
    flg_printed = 1;
end

if ( best.phi == 0 && best.f <= vtr) || elapsedTime >= maxtime || neval >= maxeval || stallgenerations >= maxstallgen
    if 0 < output_intvl && flg_printed == 0
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
        ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
    end
    finalreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptStoreBestAndFinalPopulation; % RCGA/shared/misc/ScriptStoreBestAndFinalPopulation
    return;
end


%% Second to final generations
while i < maxgen
    i = i + 1;
    flg_printed = 0;
    Population = GenerationAlternation(problem,opts,Population);
    index = RCGAfindBest(Population);
    if Population(index).phi < best.phi || ( Population(index).phi == best.phi && Population(index).f < best.f )
        stallgenerations = 0;
        if local > 0
            [Population(index), localneval] = RCGAlocalOptimize(problem, opts, Population(index));
            RCGA_LOCALNEVAL = RCGA_LOCALNEVAL + localneval;
        end
        best = Population(index);
    else
        stallgenerations = stallgenerations + 1;
    end
    
    elapsedTime = toc;
    neval = n_population + ( i - 1 ) * n_children + RCGA_LOCALNEVAL;

    if 0 < output_intvl && mod(i,output_intvl) == 0
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
        ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
        flg_printed = 1;
    end

    if (best.phi == 0 && best.f <= vtr) || elapsedTime >= maxtime || neval >= maxeval || stallgenerations >= maxstallgen
        break;
    end
end

if 0 < output_intvl && flg_printed == 0
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
end

finalreportfun(elapsedTime,i,problem,opts,Population,best);
ScriptStoreBestAndFinalPopulation; % RCGA/shared/misc/ScriptStoreBestAndFinalPopulation

