function Population = RCGA_MGG(problem, opts, Population)
% RCGA_MGG updates population by using Minimal Generation Gap (MGG).
% 
% [SYNTAX]
% Population = RCGA_MGG(problem, opts, Population)
% 
% [INPUT]
% problem    :  Problem structure:
%               - problem.n_gene: Number of decision variables.
%               - problem.n_constraint: Number of constraint functions. 
%                  For unconstained problems, this must be zero.
%               - problem.fitnessfun: Function handle for a fitness 
%                  function.
%               - problem.decodingfun: Function handle for a decoding 
%                  function.
% opts       :  Option structure:
%               - opts.n_population: Population size.
%               - opts.n_children: Number of children.
%               - opts.n_parent: Number of parents.
%               - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%               - opts.selection_type: Selection type for REXstar/JGG 
%                  (0 or 1).
%               - opts.Pf: Probability that only the objective function f 
%                  is used in comparisons of individuals in the stochastic 
%                  ranking.
%               - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                  local optimizer is used.
%               - opts.localopts: Options for the local optimizer.
%               - opts.n_generation: Number of maximum generations.
%               - opts.maxtime: Maximum time (sec).
%               - opts.maxeval: Maximum number of fitnessfun evaluations.
%               - opts.vtr: Value to be reached.
%               - opts.n_par: Number of workers in parallel computation.
%               - opts.output_intvl: Interval generation for updating the 
%                  transition file and the report file.
%               - opts.out_transition: Name of an output file called the 
%                  transition file.
%               - opts.out_best: Name of an output file called the best 
%                  individual file.
%               - opts.out_population: Name of an output file called the 
%                  final population file.
%               - opts.out_report: Name of an output file called the report
%                  file.
%               - opts.interimreportfun: Function handle for the interim 
%                  report function.
%               - opts.finalreportfun: Function handle for the final report 
%                  function.
% Population :  Population (Array of individuals).
% 
% [OUTPUT]
% Population :  Updated population (Array of individuals).
% 
% 
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p234, 2000 


%% Shortening variable names
n_population = opts.n_population;
n_children = opts.n_children;
n_constraint = problem.n_constraint;
n_gene = problem.n_gene;
Pf = opts.Pf;


%% Default number of trials of generating children
maxitr = 10; % You can change this line


%% Generating new children
ip = randperm(n_population,2);

c(1,1:n_children) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',zeros,'phi',zeros);

for i = 1 : n_children
    for j = 1 : maxitr
        ip(3) = randi(n_population);
        if ip(3) ~= ip(1) && ip(3) ~= ip(2)
            break;
        end
    end
    c(i) = RCGAgetNewChild(Population(ip(1)),Population(ip(2)),Population(ip(3)));
end
c = RCGArequestFitnessCalc(problem,opts,c);


%% Making a family
f = [ c Population(ip(1:2)) ];
f = RCGAsrsort(f,Pf);


%% Making a new population
Population(ip(1:2)) = [ f(1),  f(randi([2 n_children+2])) ];
Population = RCGAsrsort(Population,Pf);
