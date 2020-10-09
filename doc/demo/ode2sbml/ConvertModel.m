% This script demonstrates how to convert a ODE file into a SBML file.
% For the functions RCGAconvertODE2IQMmodel and IQMexportSBML, IQM Tools is
% required


% Simualtion settings
tspan = 0 : 0.2 : 20;
y0 = zeros(1,2);

% Simulate the ODE model in the original file (ExampleModel_original.m).
[T1, Y1] = ode15s(@ExampleModel_original, tspan, y0);

% Convert the original file into the SBML file.
model = RCGAconvertODE2IQMmodel('ExampleModel_original.m');
IQMexportSBML(model,'ExampleModel_converted.xml');

% Simulate the ODE model in the generated SBML file
% (ExampleModel_converted.xml). RCGAsimulate generates the intermediate
% file (ExampleModel_converted_odefun.m). The intermediate file is actually
% used for the simulation.
[T2, Y2] =  RCGAsimulate('ExampleModel_converted.xml', tspan, y0);

% Comfirm that the two simulations gave the same results.
figure;
plot(T1,Y1,'-','LineWidth',2);
hold on;
plot(T2,Y2,'--','LineWidth',2);
xlabel('Time','FontSize',11,'FontName','Arial');
ylabel('Conc.','FontSize',11,'FontName','Arial');
set(gca,'FontSize',11,'FontName','Arial');
legend('X_1 (ExampleModel\_original.m)','X_2 (ExampleModel\_original.m)',...
    'X_1 (ExampleModel\_converted.xml)','X_2 (ExampleModel\_converted.xml)'...
    ,'Location','best');
hold off;
