clearvars;
clear randn rand;

addpath('benchmark');
addpath('../app');
addpath('../source');
addpath('../source/UNDXMGG');
addpath('../source/REXstarJGG');
addpath('../source/Shared');
addpath('../source/Shared/IO');
addpath('../source/Shared/Misc');
addpath('../source/Shared/Sort');

Param = struct;
Param = getParam(Param,'Sphere');

Param.par = 0;
Param.output_intvl = 1;
Param.vtr = 1e-6;
Param.t_limit = 1e+4;
Param.t_rextar = 6.0;
Param.selection_type = 0;
Param.out_transition = 'Transition.dat';
Param.out_population = 'Population.dat';
Param.out_solution = 'Solution.dat';

Param.n_generation = 1000;
Param.n_population = 500;
Param.n_children = 200;
Param.output_intvl = 20;
Param.n_constraint = 0;
Param.n_parent = Param.n_gene + 1;
Param.selection_type = 0;
Param.Pf = 0.45;
% Param.Pf = 0;

% Param.n_gene = 2;
% Param.n_generation = 1e+5;
% Param.n_population = 200;
% Param.n_children = 200;
% Param.output_intvl = 100;


% Param = getParam(Param,'g10');
% Param.fitnessfun   = @ConstrainedSphere;
% Param.decodingfun  = @Sphere_decode;
        
% if 0 < poolsize
%     parpool(poolsize);
% end

rng(3);
% [ best, Population ] = UNDXMGG(Param);
[ best, Population ] = REXstarJGG(Param);

%%
% clear randn_test;
% Population(1).gene(1) = 0.1;
% Population(1).gene(2) = 0.2;
% Population(2).gene(1) = 0.2;
% Population(2).gene(2) = 0.3;
% Population(3).gene(1) = 0.3;
% Population(3).gene(2) = 0.4;
% for i = 1 : 100
%     c(i) = getNewChild(Population(1), Population(2), Population(3));
%     fprintf('%e\t%e\n',c(i).gene(1),c(i).gene(2));
% end
% 
% %%
% clear randn_test rand_test;
% for i = 1 : 20
%     Population(i).gene(1) = 0.1;
%     Population(i).gene(2) = 0.2;
%     Population(i).f = rand_test();
%     Population(i).phi = 1 - rand_test();
%     Population(i).g = rand_test();
% end
% Population = SRsort(Population, 0.45);
