function f = ScaledShiftedRotatedRastrigin(x)
% Unconstrained benchmark function ScaledShiftedRotatedRastrigin.
% 
% The optimum is located at x* = (4, 2, ..., 4/i, ..., 4/n) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = ScaledShiftedRotatedRastrigin(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);

% Load the rotation matrix R
% Make sure det(R) = 1 and norm(R*x) = norm(x) for any x.
switch n
    case 2
        load('R2.mat')
    case 10
        load('R10.mat')
    case 25
        load('R25.mat');
    otherwise
        error('Number of Decision Variables Must be 2, 10, or 25.');
end

% For debug. No rotation.
% R = eye(n);

if ~iscolumn(x)
    x = x';
end

M = diag(1:n);
y = M * x - 4;
z = R * y;

f = Rastrigin(z);
