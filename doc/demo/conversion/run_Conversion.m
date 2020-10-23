% This script demonstrates how to convert an ODE file into an SBML file.
% For the functions RCGAconvertODE2IQMmodel and IQMexportSBML, IQM Tools is
% required


% Simualtion settings
tspan = 0 : 0.2 : 20;
y0 = zeros(1,2);

% Simulate the ODE model in the original file (ExampleModel_original.m).
[T1, Y1] = ode15s(@Model_Example, tspan, y0);

% Convert the original file into the SBML file.
model = RCGAreadConciseODEfile('Model_Example.m');
IQMexportSBML(model,'Model_Example_converted.xml');

% Simulate the ODE model in the generated SBML file
% (Model_Example_converted.xml). RCGAsimulate generates the intermediate
% file (Model_Example_converted_odefun.m). The intermediate file is
% actually used for the simulation.
[T2, Y2] =  RCGAsimulate('Model_Example_converted.xml', tspan, y0);

% Comfirm that the two simulations gave the same results.
figure;
plot(T1,Y1,'-','LineWidth',2);
hold on;
plot(T2,Y2,'--','LineWidth',2);
legend('X_1 (Model\_Example.m)','X_2 (Model\_Example.m)',...
    'X_1 (Model\_Example\_converted.xml)','X_2 (Model\_Example\_converted.xml)'...
    ,'Location','best');
xlabel('Time');
ylabel('Concentration');
hold off;
