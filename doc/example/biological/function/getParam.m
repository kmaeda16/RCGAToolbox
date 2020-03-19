function [problem, opts] = getParam(problem_name,opts)
% Based on problem_name, getParam set fields in problem and opts.
% 
% [SYNTAX]
% [problem, opts] = getParam(problem_name,opts)
% 
% [INPUT]
% problem_name :  Name of problem
% opts         :  RCGA options. See XXXXXXXXXXX for options.
% 
% [OUTPUT]
% problem      :  Problem structure
% opts         :  RCGA options. See XXXXXXXXXXX for options.

switch problem_name
        
    case 'threestep'
        problem.fitnessfun   = @wrapper_threestep_con_mex;
        problem.decodingfun  = @threestep_decode;
        problem.n_gene       = 36;
        problem.n_constraint = 24;
        opts.vtr             = 1e-3;
        
    case 'hiv'
        problem.fitnessfun   = @wrapper_hiv_con_mex;
        problem.decodingfun  = @hiv_decode;
        problem.n_gene       = 20;
        problem.n_constraint = 12;
        opts.vtr             = 2e-2;
        
    otherwise
        error('Unexpected Problem_Name!');
end

opts.n_population = 350;
opts.n_children = 350;
opts.n_generation = 1e+8;
opts.output_intvl = 10;
opts.t_limit = 24 * 60 * 60; % 1 day
opts.localoptim = 1;
