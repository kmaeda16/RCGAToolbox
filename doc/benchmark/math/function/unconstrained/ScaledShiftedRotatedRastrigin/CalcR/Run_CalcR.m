% This script calculates an n x n rotation matrix R.

clearvars;

% n = 2;
n = 10;
% n = 25; % R25.mat was obtained from lsgo_2013_benchmarks_improved/matlab/datafiles/f04.mat

rng(1); % For Reproducibility

% An n x n rotation matrix R has n^2 (unknown) elements. So, you need n^2
% randam vectors to uniquely determine R.
global random_vs;
random_vs = rand(n,n^2);

problem.fitnessfun   = @RotationMatrix;
problem.decodingfun  = @RotationMatrix_decode;
problem.n_gene       = n ^ 2;
problem.n_constraint = 1;
opts.vtr          = 1e-9;
opts.maxgen       = 1e+8;
opts.output_intvl = 10;
opts.maxtime      = 10 * 60; % 10 min
opts.local        = 1;

Results = RCGA_REXstarJGG(problem,opts);

% Check results
R = reshape(Results.Best.x,n,n);
v = rand(n,1);
fprintf('(norm(R*v) - norm(v))^2 = %e\n',(norm(R*v) - norm(v))^2);
fprintf('Det(R)                  = %e\n',det(R));

save(sprintf('R%d.mat',n),'R');
