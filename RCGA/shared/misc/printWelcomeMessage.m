function printWelcomeMessage(problem, opts, RCGA_Type)

if isempty(strfind(RCGA_Type,'REXstarJGG')) && isempty(strfind(RCGA_Type,'UNDXMGG'))
    error('RCGA_Type must include REXstarJGG or UNDXMGG');
end

fprintf('\n');
fprintf('==========================================\n\n');
fprintf('        %s by RCGAToolbox\n\n',RCGA_Type);
fprintf('==========================================\n\n');
fprintf('\n');


disp('---------------- Problem -----------------');
fprintf('          n_gene: %d\n',problem.n_gene);
fprintf('    n_constraint: %d\n',problem.n_constraint);
fprintf('      fitnessfun: %s\n',func2str(problem.fitnessfun));
fprintf('     decodingfun: %s\n',func2str(problem.decodingfun));
disp('--------------- GA Options ---------------');
fprintf('    n_population: %d\n',opts.n_population);
fprintf('      n_children: %d\n',opts.n_children);
if isempty(strfind(RCGA_Type,'UNDXMGG'))
    fprintf('        n_parent: %d\n',opts.n_parent);
    fprintf('        t_rextar: %g\n',opts.t_rextar);
    fprintf('  selection_type: %d\n',opts.selection_type);
end
fprintf('              Pf: %g\n',opts.Pf);
fprintf('    n_generation: %d\n',opts.n_generation);
fprintf('         t_limit: %g\n',opts.t_limit);
fprintf('             vtr: %g\n',opts.vtr);
fprintf('    output_intvl: %g\n',opts.output_intvl);
fprintf('  out_transition: %s\n',opts.out_transition);
fprintf('        out_best: %s\n',opts.out_best);
fprintf('  out_population: %s\n',opts.out_population);
fprintf('      out_report: %s\n',opts.out_report);
fprintf('interimreportfun: %s\n',func2str(opts.interimreportfun));
fprintf('  finalreportfun: %s\n',func2str(opts.finalreportfun));
fprintf('             par: %d\n',opts.par);
disp('------------------------------------------');
fprintf('\n');

fprintf('%s starts.\n',RCGA_Type);