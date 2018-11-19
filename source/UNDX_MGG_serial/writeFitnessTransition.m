function writeFitnessTransition(generation, individual, SearchRegion, FileNameTransition)

fid = fopen(FileNameTransition,'a');

fprintf(fid,'%d\t',generation);
x = decodeGene2Variable(individual,SearchRegion);
fprintf(fid,'%e\t',x);
fprintf(fid,'%e\n',individual.fitness);

fclose(fid);
