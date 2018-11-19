clear all;

n = 2;

% Sphere Function (-5.12 < x < 5.12)
SearchRegion(1,:) = -5.12*ones(1,n);
SearchRegion(2,:) = 5.12*ones(1,n);

% % Rosenbrock Function (-2.048 < x < 2.048)
% SearchRegion(1,:) = -2.048*ones(1,n);
% SearchRegion(2,:) = 2.048*ones(1,n);

% % Rastrigin Function (-5.12 < x < 5.12)
% SearchRegion(1,:) = -5.12*ones(1,n);
% SearchRegion(2,:) = 5.12*ones(1,n);

% % Schwefel Function (-512 < x < 512)
% SearchRegion(1,:) = -512*ones(1,n);
% SearchRegion(2,:) = 512*ones(1,n);

[x fitness]= myGA(100,5,5,n,0,SearchRegion,'FitnessTransition.dat');
