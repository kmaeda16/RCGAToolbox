clearvars;

addpath('benchmark');
addpath('../source');
addpath('../source/UNDXMGG');
addpath('../source/REXstarJGG');
addpath('../source/Shared');
addpath('../source/Shared/IO');
addpath('../source/Shared/Misc');
addpath('../source/Shared/Sort');

% Param.fitnessfun = @k_tablet;
% Param.decodingfun = @k_tablet_decode;
% Param.fitnessfun = @Rosenbrock_chain;
% Param.decodingfun = @Rosenbrock_chain_decode;
% Param.fitnessfun = @Rosenbrock_star;
% Param.decodingfun = @Rosenbrock_star_decode;
% Param.fitnessfun = @Schwefel;
% Param.decodingfun = @Schwefel_decode;
Param.fitnessfun = @g01;
Param.decodingfun = @g01_decode;



Param.par = 0;
Param.output_intvl = 1;
Param.vtr = 1e-6;
Param.t_limit = 1e+4;
Param.t_rextar = 6.0;
Param.selection_type = 0;
Param.out_transition = 'Transition.dat';
Param.out_population = 'Population.dat';
Param.out_solution = 'Solution.dat';

Param.n_population = 5;
Param.n_parent = 3;
Param.n_children = 5;
Param.n_generation = 2;
Param.n_gene = 2;
Param.output_intvl = 1;

Param.n_population = 200;
Param.n_children = 200;
Param.n_generation = 1e+5;
Param.n_gene = 13;
Param.n_constraint = 9;
Param.n_parent = Param.n_gene + 1;
Param.output_intvl = 1;
Param.selection_type = 0;
Param.vtr = -15 * ( 1 - 1e-2 );
% Param.Pf = 0.45;
Param.Pf = 0;

% Param = getParam(Param,'Sphere');
Param = getParam(Param,'g05');




% if 0 < poolsize
%     parpool(poolsize);
% end

rng(3);
[ best, Population ] = UNDXMGG(Param);
% [ best, Population ] = REXstarJGG(Param);
