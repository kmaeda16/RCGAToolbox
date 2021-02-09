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
        
    case 'B2'
        problem.fitnessfun   = @wrapper_b2_obj;
        problem.decodingfun  = @BioPreDynBench_decode;
        problem.n_gene       = 116;
        problem.n_constraint = 0;
        opts.output_intvl    = 10;
        opts.vtr             = -inf; % 2.3390e+2 * 2 = 4.678000e+02 (Villaverde, BMC Syst Biol, 2015, Table 2 Jf)
        load('b2_bounds.mat','lb','ub');
        
    case 'B4'
        problem.fitnessfun   = @wrapper_b4_obj;
        problem.decodingfun  = @BioPreDynBench_decode;
        problem.n_gene       = 117;
        problem.n_constraint = 0;
        opts.output_intvl    = 10;
        opts.vtr             = -inf; % 4.5718e+1 * 2 = 9.143600e+01 (Villaverde, BMC Syst Biol, 2015, Table 2 Jf)
        load('b4_bounds.mat','lb','ub');
        
    case 'B5'
        problem.fitnessfun   = @wrapper_b5_obj;
        problem.decodingfun  = @BioPreDynBench_decode;
        problem.n_gene       = 86;
        problem.n_constraint = 0;
        opts.output_intvl    = 1;
        opts.vtr             = -inf; % 3.0725e+3 * 2 = 6.145000e+03 (Villaverde, BMC Syst Biol, 2015, Table 2 Jf)
        load('b5_bounds.mat','lb','ub');
        
    case 'B6'
        problem.fitnessfun   = @wrapper_b6_obj;
        problem.decodingfun  = @BioPreDynBench_decode;
        problem.n_gene       = 37;
        problem.n_constraint = 0;
        opts.output_intvl    = 10;
        opts.vtr             = -inf; % 1.0833e+5 * 2 = 2.166600e+05 (Villaverde, BMC Syst Biol, 2015, Table 2 Jf)
        load('b6_bounds.mat','lb','ub');
        
    otherwise
        error('Unexpected Problem_Name!');
end

opts.maxgen = 1e+8;
opts.maxtime = 24 * 60 * 60; % 1 day
% opts.maxtime = 5; % 5 sec
opts.maxeval = 1e+8;
opts.local = 1;
opts.Pf = 0;

global lb_global ub_global;
lb_global = lb;
ub_global = ub;
