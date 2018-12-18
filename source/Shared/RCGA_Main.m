function [ best, Population ] = RCGA_Main(Param, GenerationAlternation)

tic;

n_generation = Param.n_generation;
n_constraint = Param.n_constraint;
output_intvl = Param.output_intvl;
fitnessfun = Param.fitnessfun;
vtr = Param.vtr;
t_limit = Param.t_limit;
out_transition = Param.out_transition;
out_population = Param.out_population;
out_solution = Param.out_solution;

flg_printed = 0;



i = 1;
Population = getInitPopulation(Param);
index = findBest(Population);
best = Population(index);

if 0 < output_intvl
    elapsedTime = toc;
    printTransition(elapsedTime,i,best);
    writeTransition(elapsedTime,i,Param,best);
    flg_printed = 1;
end

elapsedTime = toc;
if ( best.phi == 0 && best.f <= vtr) || elapsedTime >= t_limit
    if 0 < output_intvl && flg_printed == 0
        printTransition(elapsedTime,i,best);
        writeTransition(elapsedTime,i,Param,best);
    end
    writePopulation(Param,Population);
    writeSolution(elapsedTime,Param,best);
    return;
end

while i < n_generation
    i = i + 1;
    flg_printed = 0;
    Population = GenerationAlternation(Param,Population);
    index = findBest(Population);
    if Population(index).phi < best.phi || ( Population(index).phi == best.phi && Population(index).f < best.f )
        best = Population(index);
    end
    if 0 < output_intvl && mod(i,output_intvl) == 0
        elapsedTime = toc;
        printTransition(elapsedTime,i,best);
        writeTransition(elapsedTime,i,Param,best);
        flg_printed = 1;
    end
    if (best.phi == 0 && best.f <= vtr) || toc >= t_limit
        break;
    end
end

if 0 < output_intvl && flg_printed == 0
    elapsedTime = toc;
    printTransition(elapsedTime,i,best);
    writeTransition(elapsedTime,i,Param,best);
end

% for j = 1 : Param.n_population
%     fprintf('%e\t%e\t%e\t%e\t%e\t%e\t%e\n',i,Population(j).gene(1),Population(j).gene(2),Population(j).g(1),Population(j).g(2),Population(j).f,Population(j).phi);
% end

writePopulation(Param,Population);
writeSolution(elapsedTime,Param,best);
