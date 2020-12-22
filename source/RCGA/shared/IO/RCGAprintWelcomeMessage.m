function RCGAprintWelcomeMessage(problem, opts, RCGAfun)
% RCGAprintWelcomeMessage prints the information on the problem and the
% options.
% 
% [SYNTAX]
% RCGAprintWelcomeMessage(problem, opts, RCGA_Type)
% 
% [INPUT]
% problem   :  Problem structure:
%              - problem.n_gene: Number of decision variables.
%              - problem.n_constraint: Number of constraint functions. For 
%                 unconstained problems, this must be zero.
%              - problem.fitnessfun: Function handle for a fitness 
%                 function.
%              - problem.decodingfun: Function handle for a decoding 
%                 function.
% opts      :  Option structure:
%              - opts.n_population: Population size.
%              - opts.n_children: Number of children.
%              - opts.n_parent: Number of parents.
%              - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%              - opts.selection_type: Selection type for REXstar/JGG 
%                 (0 or 1).
%              - opts.Pf: Probability that only the objective function f is
%                 used in comparisons of individuals in the stochastic 
%                 ranking.
%              - opts.local: Local optimizer (0 or 1). If it is 1, the 
%                 local optimizer is used.
%              - opts.localopts: Options for the local optimizer.
%              - opts.maxgen: Maximum number of generations.
%              - opts.maxtime: Maximum time (sec).
%              - opts.maxeval: Maximum number of fitnessfun evaluations.
%              - opts.vtr: Value to be reached.
%              - opts.n_par: Number of workers in parallel computation.
%              - opts.output_intvl: Interval generation for updating the 
%                 transition file and the report file.
%              - opts.out_transition: Name of an output file called the 
%                 transition file.
%              - opts.out_best: Name of an output file called the best 
%                 individual file.
%              - opts.out_population: Name of an output file called the 
%                 final population file.
%              - opts.out_report: Name of an output file called the report 
%                 file.
%              - opts.interimreportfun: Function handle for the interim 
%                 report function.
%              - opts.finalreportfun: Function handle for the final report 
%                 function.
% RCGAfun   :  Function handle for RCGA_UNDXMGG, RCGA_REXstarJGG,
%              RCGA_CustomRCGA.


RCGA_Type = func2str(RCGAfun);


%% Checking RCGA_Type
if ~strcmp(RCGA_Type,'RCGA_UNDXMGG') && ~strcmp(RCGA_Type,'RCGA_REXstarJGG') && ~strcmp(RCGA_Type,'RCGA_CustomRCGA')
    error('RCGA_Type must be ''RCGA_UNDXMGG'', ''RCGA_REXstarJGG'', or ''RCGA_CustomRCGA''.');
end


%% Prenting details
fprintf('\n');
fprintf('================================================\n\n');
fprintf('           %s by RCGAToolbox\n\n',RCGA_Type);
fprintf('================================================\n\n');
fprintf('\n');


disp('------------------- Problem --------------------');
fprintf('          n_gene :  %d\n',problem.n_gene);
fprintf('    n_constraint :  %d\n',problem.n_constraint);
fprintf('      fitnessfun :  %s\n',func2str(problem.fitnessfun));
fprintf('     decodingfun :  %s\n',func2str(problem.decodingfun));
disp('------------------- Options --------------------');
fprintf('    n_population :  %d\n',opts.n_population);
fprintf('      n_children :  %d\n',opts.n_children);
if strcmp(RCGA_Type,'RCGA_REXstarJGG')
    fprintf('        n_parent :  %d\n',opts.n_parent);
    fprintf('       t_rexstar :  %g\n',opts.t_rexstar);
    fprintf('  selection_type :  %d\n',opts.selection_type);
end
fprintf('              Pf :  %g\n',opts.Pf);
fprintf('           local :  %g\n',opts.local);
fprintf('          maxgen :  %d\n',opts.maxgen);
fprintf('         maxtime :  %g\n',opts.maxtime);
fprintf('         maxeval :  %g\n',opts.maxeval);
fprintf('             vtr :  %g\n',opts.vtr);
fprintf('           n_par :  %d\n',opts.n_par);
fprintf('    output_intvl :  %g\n',opts.output_intvl);
fprintf('  out_transition :  %s\n',opts.out_transition);
fprintf('        out_best :  %s\n',opts.out_best);
fprintf('  out_population :  %s\n',opts.out_population);
fprintf('      out_report :  %s\n',opts.out_report);
fprintf('interimreportfun :  %s\n',func2str(opts.interimreportfun));
fprintf('  finalreportfun :  %s\n',func2str(opts.finalreportfun));
disp('------------------------------------------------');
fprintf('\n');

fprintf('%s starts.\n',RCGA_Type);
