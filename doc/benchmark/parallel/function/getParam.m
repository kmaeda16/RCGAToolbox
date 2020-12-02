function [problem, opts] = getParam(problem_name,opts)
% Based on problem_name, getParam sets the fields of problem and opts.
% 
% [SYNTAX]
% [problem, opts] = getParam(problem_name,opts)
% 
% [INPUT]
% problem_name :  Name of problem.
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
%                 - opts.n_generation: Number of maximum generations.
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
% problem      :  Problem structure:
%                 - problem.n_gene: Number of decision variables.
%                 - problem.n_constraint: Number of constraint functions. 
%                    For unconstained problems, this must be zero.
%                 - problem.fitnessfun: Function handle for a fitness 
%                    function.
%                 - problem.decodingfun: Function handle for a decoding 
%                    function.
% opts         :  Modified option structure.


switch problem_name
        
    case 'threestep'
        problem.fitnessfun   = @wrapper_threestep_con_mex;
        problem.decodingfun  = @threestep_decode;
        problem.n_gene       = 36;
        problem.n_constraint = 24;
        opts.vtr             = -inf; % 1e-3;
        
    case 'hiv'
        problem.fitnessfun   = @wrapper_hiv_con_mex;
        problem.decodingfun  = @hiv_decode;
        problem.n_gene       = 20;
        problem.n_constraint = 12;
        opts.vtr             = -inf; % 2e-2;
        
    otherwise
        error('Unexpected Problem_Name!');
end

opts.n_population = 350;
opts.n_children = 350;
opts.n_generation = 1000; % 1e+8;
opts.output_intvl = 10;
opts.maxtime = 24 * 60 * 60; % 1 day
% opts.maxtime = 12 * 60 * 60; % 12 hr
% opts.maxtime = 2 * 60 * 60; % 2 hr
% opts.maxtime = 5; % 5 sec
opts.maxeval = 1e+8;
opts.local = 1;
