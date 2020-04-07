function Results = RCGA_Main(problem, opts, GenerationAlternation)
% RCGA_Main is the main function of real-coded genetic algorithms.
% 
% [SYNTAX]
% Results = RCGA_Main(problem, opts, GenerationAlternation)
% 
% [INPUT]
% problem               :  Problem structure.
% opts                  :  RCGA options. See XXXXXXXXXXX for options.
% GenerationAlternation :  Function handle to REXstarJGG or UNDXMGG
% 
% [OUTPUT]
% Results               :  Structure with results


%% Getting the time RCGA starts
tic;


%% Shortening variable names
decodingfun = problem.decodingfun;
n_constraint = problem.n_constraint;
interimreportfun = opts.interimreportfun;
finalreportfun = opts.finalreportfun;
n_generation = opts.n_generation;
output_intvl = opts.output_intvl;
n_population = opts.n_population;
out_report = opts.out_report;
n_children = opts.n_children;
vtr = opts.vtr;
maxtime = opts.maxtime;
maxeval = opts.maxeval;
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

if ( best.phi == 0 && best.f <= vtr) || elapsedTime >= maxtime || neval >= maxeval
    if 0 < output_intvl && flg_printed == 0
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
        ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
    end
    finalreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptStoreBestAndFinalPopulation; % RCGA/shared/misc/ScriptStoreBestAndFinalPopulation
    return;
end


%% Second to final generations
while i < n_generation
    i = i + 1;
    flg_printed = 0;
    Population = GenerationAlternation(problem,opts,Population);
    index = RCGAfindBest(Population);
    if Population(index).phi < best.phi || ( Population(index).phi == best.phi && Population(index).f < best.f )
        if local > 0
            [Population(index), localneval] = RCGAlocalOptimize(problem, opts, Population(index));
            RCGA_LOCALNEVAL = RCGA_LOCALNEVAL + localneval;
        end
        best = Population(index);
    end
    
    elapsedTime = toc;
    neval = n_population + ( i - 1 ) * n_children + RCGA_LOCALNEVAL;

    if 0 < output_intvl && mod(i,output_intvl) == 0
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
        ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
        flg_printed = 1;
    end

    if (best.phi == 0 && best.f <= vtr) || elapsedTime >= maxtime || neval >= maxeval
        break;
    end
end

if 0 < output_intvl && flg_printed == 0
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
end

finalreportfun(elapsedTime,i,problem,opts,Population,best);
ScriptStoreBestAndFinalPopulation; % RCGA/shared/misc/ScriptStoreBestAndFinalPopulation

