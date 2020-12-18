function RCGAinterimreportfun_PE(elapsedTime, generation, problem, opts, Population, best, Simulation, model, mst, simopts)
% RCGAinterimreportfun_PE is called every output_intvl generations,
% printing messages on the display, making the transition file, showing a
% plot with simulation and experimental data.
% 
% [SYNTAX]
% RCGAinterimreportfun_PE(elapsedTime, generation, problem, opts, ...
%                         Population, best, Simulation, model, mst, ...
%                         simopts)
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
% opts        :  Option structure:
%                - opts.n_population: Population size.
%                - opts.n_children: Number of children.
%                - opts.n_parent: Number of parents.
%                - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%                - opts.selection_type: Selection type for REXstar/JGG 
%                   (0 or 1).
%                - opts.Pf: Probability that only the objective function f 
%                   is used in comparisons of individuals in the stochastic 
%                   ranking.
%                - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                   local optimizer is used.
%                - opts.localopts: Options for the local optimizer.
%                - opts.n_generation: Number of maximum generations.
%                - opts.maxtime: Maximum time (sec).
%                - opts.maxeval: Maximum number of fitnessfun evaluations.
%                - opts.vtr: Value to be reached.
%                - opts.n_par: Number of workers in parallel computation.
%                - opts.output_intvl: Interval generation for updating the 
%                   transition file and the report file.
%                - opts.out_transition: Name of an output file called the 
%                   transition file.
%                - opts.out_best: Name of an output file called the best 
%                   individual file.
%                - opts.out_population: Name of an output file called the 
%                   final population file.
%                - opts.out_report: Name of an output file called the 
%                   report file.
%                - opts.interimreportfun: Function handle for the interim 
%                   report function.
%                - opts.finalreportfun: Function handle for the final 
%                   report function.
% Population  :  Population (Array of individuals).
% best        :  Structure of the the best individual.
%                - best.gene: Decision variable vector.
%                - best.g: Constraint function value vector.
%                - best.f: Fitness function value.
%                - best.phi: Penalty function value.
% Simulation  :  Function handle for RCGAsimulateODEXX, RCGAsimulateSTB, or
%                RCGAsimulateMEX.
% model       :  Function handle for an ODE function (IQM Tools format) or 
%                a MEXed model.
% mst         :  Experimental data (An IQMmeasurement object).
% simopts     :  Solver option structure. The fields depend on fast_flag. 
%                For fast_flag = 0, 1, and 2, see 'help RCGAsimulateODEXX', 
%                'help RCGAsimulateSTB', 'help RCGAsimulateMEX', 
%                respectively.


%% Making outputs
RCGAprintTransition(elapsedTime,generation,problem,best);
RCGAwriteTransition(elapsedTime,generation,problem,opts,best);


%% Decoding gene to x
param = problem.decodingfun(best.gene);


%% Converting IQMmeasurement into structure
mst = struct(mst);


%% Checking time errors
n_point = 100; % You can change this line

t0 = 0;
if t0 < mst.time(1)
    tspan = linspace(t0 ,mst.time(end),n_point)';
elseif t0 == mst.time(1)
    tspan = linspace(mst.time(1),mst.time(end),n_point)';
else
    error('Time of the first experimental datapoint should be AFTER or EQUAL TO time 0');
end


%% Running simulation
y0 = model();
[ T, Y ] = feval(Simulation, model, tspan, y0, param, simopts);


%% Making plot
T_exp = mst.time;
n_col = length(mst.data);
n_row = length(mst.data(1).values);
Y_exp = zeros(n_row,n_col);
for i = 1 : n_col
    Y_exp(:,i) = mst.data(i).values;
end
T_sim = T;
Y_sim = Y;
statename = model('states');
RCGAplotter(T_sim,Y_sim,T_exp,Y_exp,'Time','Value',statename);
title(sprintf('Generation = %d',generation));
drawnow;
