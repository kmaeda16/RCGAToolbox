function doBenchmark(idum)

addpath('benchmark');
addpath('../source');
addpath('../source/UNDXMGG');
addpath('../source/REXstarJGG');
addpath('../source/Shared');
addpath('../source/Shared/IO');
addpath('../source/Shared/Misc');
addpath('../source/Shared/Sort');

BENCHMARK1 = {'Sphere','ScaledSphere','Ellipsoid','Cigar','k_tablet', ...
    'MMbenchmark','Rosenbrock_star','Rosenbrock_chain','Ackley','Bohachevsky', ...
    'Rastrigin','Schaffer','Schwefel'};
BENCHMARK2 = {'g01','g02','g03','g04','g05', ...
    'g06','g07', 'g08','g09','g10', ...
    'g11','g12','g13'};


rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);

Param.Pf = 0;
for Problem_Name = BENCHMARK1
    
    fprintf('%s\n',char(Problem_Name));
    Param = getParam(Param,char(Problem_Name));
    
    Param.n_population = 200;
    Param.n_children = 200;
    Param.n_generation = 1e+5;
    Param.n_parent = Param.n_gene + 1;
    Param.t_rextar = 6.0;
    Param.output_intvl = 1e+8;
    Param.selection_type = 0;
    Param.t_limit = 24 * 60 * 60;
    Param.par = 0;
    Param.out_population = 'None'; % 'Population.dat';
    Param.out_solution = 'None'; % 'Solution.dat';
    
    Param.out_transition = sprintf('UNDXMGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    UNDXMGG(Param);
    
    Param.out_transition = sprintf('REXstarJGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    REXstarJGG(Param);
    
end

for Problem_Name = BENCHMARK2
    
    fprintf('%s\n',char(Problem_Name));
    Param = getParam(Param,char(Problem_Name));
    
    Param.n_population = 200;
    Param.n_children = 200;
    Param.n_generation = 1e+5;
    Param.n_parent = Param.n_gene + 1;
    Param.t_rextar = 6.0;
    Param.output_intvl = 1e+8;
    Param.selection_type = 0;
    Param.t_limit = 24 * 60 * 60;
    Param.par = 0;
    Param.out_population = 'None'; % 'Population.dat';
    Param.out_solution = 'None'; % 'Solution.dat';
    
    Param.Pf = 0;
    Param.out_transition = sprintf('UNDXMGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    UNDXMGG(Param);
    Param.out_transition = sprintf('REXstarJGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    REXstarJGG(Param);
    
    Param.Pf = 0.45;
    Param.out_transition = sprintf('UNDXMGG_%s_SR_transition_%d.dat',char(Problem_Name),idum);
    UNDXMGG(Param);
    Param.out_transition = sprintf('REXstarJGG_%s_SR_transition_%d.dat',char(Problem_Name),idum);
    REXstarJGG(Param);
    
end
