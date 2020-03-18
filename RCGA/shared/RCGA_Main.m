function Results = RCGA_Main(problem, opts, GenerationAlternation)
% checkInputs checks problem and opts and sets default values
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
t_limit = opts.t_limit;
n_localoptimind = opts.n_localoptimind;


%% Setting flg_printed
flg_printed = 0; % flg_printed == 1 indicates output is made


%% First generation
i = 1;
Population = RCGAgetInitPopulation(problem,opts);
% if n_localoptimind > 0
%     Population = RCGArequestLocalOptimize(problem,opts,Population);
% end
index = RCGAfindBest(Population);
if n_localoptimind > 0
    Population(index) = RCGAlocalOptimize(problem, opts, Population(index));
end
best = Population(index);

if 0 < output_intvl
    elapsedTime = toc;
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptInitTransition; % RCGA/shared/misc/ScriptInitTransition
    ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
    flg_printed = 1;
end

elapsedTime = toc;
if ( best.phi == 0 && best.f <= vtr) || elapsedTime >= t_limit
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
%     if n_localoptimind > 0
%         Population = RCGArequestLocalOptimize(problem,opts,Population);
%     end
    index = RCGAfindBest(Population);
    if Population(index).phi < best.phi || ( Population(index).phi == best.phi && Population(index).f < best.f )
        if n_localoptimind > 0
            Population(index) = RCGAlocalOptimize(problem, opts, Population(index));
        end
        best = Population(index);
    end
    if 0 < output_intvl && mod(i,output_intvl) == 0
        elapsedTime = toc;
        interimreportfun(elapsedTime,i,problem,opts,Population,best);
        ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
        flg_printed = 1;
    end
    if (best.phi == 0 && best.f <= vtr) || toc >= t_limit
        break;
    end
end

if 0 < output_intvl && flg_printed == 0
    elapsedTime = toc;
    interimreportfun(elapsedTime,i,problem,opts,Population,best);
    ScriptStoreTransition; % RCGA/shared/misc/ScriptStoreTransition
end

finalreportfun(elapsedTime,i,problem,opts,Population,best);
ScriptStoreBestAndFinalPopulation; % RCGA/shared/misc/ScriptStoreBestAndFinalPopulation

