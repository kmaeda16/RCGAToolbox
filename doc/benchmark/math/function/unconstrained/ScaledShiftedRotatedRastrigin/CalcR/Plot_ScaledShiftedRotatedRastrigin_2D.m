% This script makes a surface plot to show how 2-dimensional
% ScaledShiftedRotatedRastrigin looks like.

clearvars;

addpath('..');

n = 80; % Number of grids
lb = -5.12;
ub =  5.12;

space = linspace(lb,ub,n);

[X,Y] = meshgrid(space,space);

Z = zeros(n,n);
for i = 1 : n
    for j = 1 : n
        x = [ X(i,j), Y(i,j) ];
        Z(i,j) = ScaledShiftedRotatedRastrigin(x);
    end
end

figure;
surf(X,Y,Z);
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');
xlim([-6, 6]);
ylim([-6, 6]);

rmpath('..');
