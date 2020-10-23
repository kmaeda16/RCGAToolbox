% This script demonstrates how to use CVODE by SundialsTB.


clearvars;


% =============== Time ================ %
tspan = 0 : 0.2 : 20;

% ========= Initial Condition ========= %
y0(1) = 0; % X1
y0(2) = 0; % X2

% ============ Simulation ============= %
tic;
% [ T, Y ] = ode45(@Example_Model, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23(@Example_Model, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode113(@Example_Model, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode15s(@Example_Model, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23s(@Example_Model, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23t(@Example_Model, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23tb(@Example_Model, tspan, y0); % MATLAB Built-in
[ T, Y ] = odestb(@ExampleModel_original, tspan, y0); % SundialsTB is required
toc

% =============== Figure ============== %
figure;
plot(T,Y,'-','LineWidth',2);
legend('X_1','X_2','Location','best');
xlabel('Time');
ylabel('Concentration');
