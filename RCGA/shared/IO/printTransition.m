function printTransition(elapsedTime, generation, problem, chrom)

fprintf('Elapsed Time = %e, Generation = %d, f = %e',elapsedTime,generation,chrom.f);

if problem.n_constraint > 0
    fprintf(', phi = %e\n',chrom.phi);
else
    fprintf('\n');
end
