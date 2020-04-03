function x = decodingExample(gene)

lb = 0;
ub = 10;
x = gene * ( ub - lb ) + lb;

% 
% lb = -2;
% ub = 1;
% x = 10 .^ ( gene * ( ub - lb ) + lb );

% x(1) = 0;
% x(end) = 1;
% x = [0, 5, 5.5, 4, 0.5, 0.1, 0.1, 0.01, 1];
