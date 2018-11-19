function fitness = getFitness(x)

n = length(x);

% Sphere Function (-5.12 < x < 5.12)
fitness = sum(x.^2);

% % Rosenbrock Function (-2.048 < x < 2.048)
% fitness = 0;
% for i = 1:n-1
%     fitness = fitness + 100*(x(i+1)-x(i)^2)^2 + (x(i)-1)^2;
% end

% % Rastrigin Function (-5.12 < x < 5.12)
% fitness = 10*n;
% for i = 1:n
%     fitness = fitness + x(i)^2 - 10*cos(2*pi*x(i));
% end

% % Schwefel Function (-512 < x < 512)
% fitness = 418.9828873*n;
% for i = 1:n
%     fitness = fitness + x(i)*sin(abs(x(i))^0.5);
% end
