% This script demonstrates how to use CVODE by SundialsTB. SundialsTB is
% required to use the function odestb.


clearvars;


% =============== Time ================ %
tspan = 0 : 0.1 : 10;

% ========= Initial Condition ========= %
y0(1) = 0; % X1
y0(2) = 0; % X2

% ============ Simulation ============= %
tic;
% [ T, Y ] = ode45(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode113(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode15s(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23s(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23t(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
% [ T, Y ] = ode23tb(@Model_Example_conciseOdefun, tspan, y0); % MATLAB Built-in
[ T, Y ] = odestb(@Model_Example_conciseOdefun, tspan, y0); % CVODE provided by SundialsTB
toc

% =============== Figure ============== %
figure;
plot(T,Y,'-','LineWidth',2);
legend('X_1','X_2','Location','best');
xlabel('Time');
ylabel('Concentration');
