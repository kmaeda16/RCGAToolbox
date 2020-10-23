% This script demonstrates how to convert an ODE file into an SBML file.
% For the functions RCGAconvertODE2IQMmodel and IQMexportSBML, IQM Tools is
% required


% % Simulation settings
% tspan = 0 : 0.2 : 20;
% y0 = zeros(1,2);
% 
% % Simulate the ODE model in the original file (ExampleModel_original.m).
% [T1, Y1] = ode15s(@Model_Example, tspan, y0);
% 
% % Convert the original file into the SBML file.
% model = RCGAreadConciseODEfile('Model_Example.m');
% IQMexportSBML(model,'Model_Example_converted.xml');
% 
% % Simulate the ODE model in the generated SBML file
% % (Model_Example_converted.xml). RCGAsimulate generates the intermediate
% % file (Model_Example_converted_odefun.m). The intermediate file is
% % actually used for the simulation.
% [T2, Y2] =  RCGAsimulate('Model_Example_converted.xml', tspan, y0);
% 
% % Comfirm that the two simulations gave the same results.
% figure;
% plot(T1,Y1,'-','LineWidth',2);
% hold on;
% plot(T2,Y2,'--','LineWidth',2);
% legend('X_1 (Model\_Example.m)','X_2 (Model\_Example.m)',...
%     'X_1 (Model\_Example\_converted.xml)','X_2 (Model\_Example\_converted.xml)'...
%     ,'Location','best');
% xlabel('Time');
% ylabel('Concentration');
% hold off;

figure;
tspan = [0 10];
y0(1) = 0.000000e+00;
y0(2) = 0.000000e+00;
[T, Y] = ode15s(@Model_Example_conciseOdefun,tspan,y0);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% Concise ODE file (RCGAToolbox format) --> IQMmodel object
% Model_Example_conciseOdefun.m --> iqmmodel1
iqmmodel1 = RCGAreadConciseODEfile('Model_Example_conciseOdefun.m');

figure;
[T, Y] = RCGAsimulate(iqmmodel1);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% IQMmodel object --> Concise ODE file (RCGAToolbox format)
% iqmmodel1 --> Model_Example_conciseOdefun_1.m
RCGAcreateConciseODEfile(iqmmodel1,'Model_Example_conciseOdefun_1.m');

figure;
tspan = [0 10];
y0(1) = 0.000000e+00;
y0(2) = 0.000000e+00;
[T, Y] = ode15s(@Model_Example_conciseOdefun_1,tspan,y0);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% IQM model object --> SBML file
% iqmmodel1 -> Model_Example_SBML.xml
IQMexportSBML(iqmmodel1,'Model_Example_SBML.xml');

figure;
[T2, Y2] =  RCGAsimulate('Model_Example_SBML.xml'); % Model_Example_SBML_odefun.m automatically created.
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% SBML file --> Concise ODE file (RCGAToolbox format)
% Model_Example_SBML.xml --> Model_Example_conciseOdefun_2.m
RCGAcreateConciseODEfile('Model_Example_SBML.xml','Model_Example_conciseOdefun_2.m');

figure;
tspan = [0 10];
y0(1) = 0.000000e+00;
y0(2) = 0.000000e+00;
[T, Y] = ode15s(@Model_Example_conciseOdefun_2,tspan,y0);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% SBML file --> IQMmodel object
% Model_Example_SBML.xml --> iqmmodel2
iqmmodel2 = IQMmodel('Model_Example_SBML.xml');

figure;
[T, Y] = RCGAsimulate(iqmmodel2); % Model_Example_SBML_odefun.m automatically created.
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% IQMmodel object --> ODE file (IQM Tools format)
% iqmmodel1 --> Model_Example_odefun_1.m
RCGAcreateODEfile(iqmmodel1,'Model_Example_odefun_1.m');
% IQMcreateODEfile(iqmmodel1,'Model_Example_odefun_1');

figure;
[T, Y] = RCGAsimulate(@Model_Example_odefun_1);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% SBML file --> ODE file (IQM Tools format)
% Model_Example_SBML.xml --> Model_Example_odefun_2.m
RCGAcreateODEfile('Model_Example_SBML.xml','Model_Example_odefun_2.m');

figure;
[T, Y] = RCGAsimulate(@Model_Example_odefun_2);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% IQMmodel object --> C source code
% iqmmodel1 --> Model_Example_odefun_2.m
IQMmakeMEXmodel(iqmmodel1,'Model_Example_C',1);
% RCGAmakeMEXmodel(iqmmodel1,'Model_Example_C',1);

figure;
[T, Y] = RCGAsimulate('Model_Example_C.c'); % Model_Example_C.mex* automatically generated.
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');

%% C source code --> MEX file
% Model_Example_C.c + Model_Example_C.h --> Model_Example_C.mex*
mexcompileIQM('Model_Example_C');

figure;
[T, Y] = RCGAsimulate(@Model_Example_C);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% IQMmodel --> MEX file
% iqmmodel1 --> Model_Example_mex_1.mex*
RCGAmakeMEXmodel(iqmmodel1,'Model_Example_mex_1');
% IQMmakeMEXmodel(iqmmodel1,'Model_Example_mex_1');

figure;
[T, Y] = RCGAsimulate(@Model_Example_mex);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');


%% SBML file --> MEX file
% Model_Example_SBML.xml --> Model_Example_mex_2.mex*
RCGAmakeMEXmodel('Model_Example_SBML.xml','Model_Example_mex_2');

figure;
[T, Y] = RCGAsimulate(@Model_Example_mex_2);
plot(T, Y);
xlabel('Time');
ylabel('States');
legend('X1','X2');
