function Interimreportfun_Example(elapsedTime, generation, problem, opts, Population, best)
% Interimreportfun_Example is called every output_intvl generations,
% printing messages on the display and making the transition file.
% 
% [SYNTAX]
% Interimreportfun_Example(elapsedTime, generation, problem, opts, Population, best)
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


global experimentaldata;


%% Making outputs
RCGAprintTransition(elapsedTime,generation,problem,best);
RCGAwriteTransition(elapsedTime,generation,problem,opts,best);


%% Decoding gene to x
param = problem.decodingfun(best.gene);


%% Checking time errors
n_point = 100; % You can change this line
tspan = linspace(experimentaldata(1,1),experimentaldata(end,1),n_point)';


%% Running simulation
y0 = [0 0];
[ T, Y ] = RCGAsimulateODEXX(@Model_Example_odefun, tspan, y0, param);


%% Making plot
T_exp = experimentaldata(:,1);
Y_exp = experimentaldata(:,2:end);
T_sim = T;
Y_sim = Y;
statename = {'X1','X2'};
RCGAplotter(T_sim,Y_sim,T_exp,Y_exp,'Time','Value',statename);
title(sprintf('Generation = %d',generation));
drawnow;
