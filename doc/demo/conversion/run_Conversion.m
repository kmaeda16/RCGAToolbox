% This script demonstrates different format convresions. IQM Tools are
% required for the demonstration.

clear all;

%% Simulate the model in Model_Example_conciseOdefun.m
figure;
tspan = 0 : 0.1 : 10;
y0(1) = 0.000000e+00;
y0(2) = 0.000000e+00;
[T, Y] = ode23s(@Model_Example_conciseOdefun,tspan,y0);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('Original');


%% Concise ODE file (RCGAToolbox format) --> IQMmodel object
% Model_Example_conciseOdefun.m --> iqmmodel1
iqmmodel1 = RCGAreadConciseODEfile('Model_Example_conciseOdefun.m');

figure;
[T, Y] = RCGAsimulate(iqmmodel1); % Model_Example_odefun.m automatically created.
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('Concise ODE file (RCGAToolbox format) --> IQMmodel object');


%% IQMmodel object --> Concise ODE file (RCGAToolbox format)
% iqmmodel1 --> Model_Example_conciseOdefun_1.m
RCGAcreateConciseODEfile(iqmmodel1,'Model_Example_conciseOdefun_1.m');

figure;
tspan = 0 : 0.1 : 10;
y0(1) = 0.000000e+00;
y0(2) = 0.000000e+00;
[T, Y] = ode23s(@Model_Example_conciseOdefun_1,tspan,y0);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('IQMmodel object --> Concise ODE file (RCGAToolbox format)');


%% IQMmodel object --> SBML file
% iqmmodel1 -> Model_Example_SBML.xml
IQMexportSBML(iqmmodel1,'Model_Example_SBML.xml');

figure;
[T, Y] = RCGAsimulate('Model_Example_SBML.xml'); % Model_Example_SBML_odefun.m automatically created.
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('IQMmodel object --> SBML file');


%% SBML file --> Concise ODE file (RCGAToolbox format)
% Model_Example_SBML.xml --> Model_Example_conciseOdefun_2.m
RCGAcreateConciseODEfile('Model_Example_SBML.xml','Model_Example_conciseOdefun_2.m');

figure;
tspan = 0 : 0.1 : 10;
y0(1) = 0.000000e+00;
y0(2) = 0.000000e+00;
[T, Y] = ode23s(@Model_Example_conciseOdefun_2,tspan,y0);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('SBML file --> Concise ODE file (RCGAToolbox format)');


%% SBML file --> IQMmodel object
% Model_Example_SBML.xml --> iqmmodel2
iqmmodel2 = IQMmodel('Model_Example_SBML.xml');

figure;
[T, Y] = RCGAsimulate(iqmmodel2); % Model_Example_SBML_odefun.m automatically created.
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('SBML file --> IQMmodel object');


%% IQMmodel object --> ODE file (IQM Tools format)
% iqmmodel1 --> Model_Example_odefun_1.m
RCGAcreateODEfile(iqmmodel1,'Model_Example_odefun_1.m');
% IQMcreateODEfile(iqmmodel1,'Model_Example_odefun_1');

figure;
[T, Y] = RCGAsimulate(@Model_Example_odefun_1);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('IQMmodel object --> ODE file (IQM Tools format)');


%% SBML file --> ODE file (IQM Tools format)
% Model_Example_SBML.xml --> Model_Example_odefun_2.m
RCGAcreateODEfile('Model_Example_SBML.xml','Model_Example_odefun_2.m');

figure;
[T, Y] = RCGAsimulate(@Model_Example_odefun_2);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('SBML file --> ODE file (IQM Tools format)');


%% IQMmodel object --> C source code
% iqmmodel1 --> Model_Example_C.c + Model_Example_C.h
RCGAmakeMEXmodel(iqmmodel1,'Model_Example_C',1);
% IQMmakeMEXmodel(iqmmodel1,'Model_Example_C',1);

figure;
[T, Y] = RCGAsimulate('Model_Example_C.c'); % Model_Example_C.mex* automatically generated.
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('IQMmodel object --> C source code');

%% C source code --> MEX file
% Model_Example_C.c + Model_Example_C.h --> Model_Example_C.mex*
mexcompileIQM('Model_Example_C');

figure;
[T, Y] = RCGAsimulate(@Model_Example_C);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('C source code --> MEX file');


%% IQMmodel --> MEX file
% iqmmodel1 --> Model_Example_mex_1.mex*
RCGAmakeMEXmodel(iqmmodel1,'Model_Example_mex_1');
% IQMmakeMEXmodel(iqmmodel1,'Model_Example_mex_1');

figure;
[T, Y] = RCGAsimulate(@Model_Example_mex_1);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('IQMmodel --> MEX file');


%% SBML file --> MEX file
% Model_Example_SBML.xml --> Model_Example_mex_2.mex*
RCGAmakeMEXmodel('Model_Example_SBML.xml','Model_Example_mex_2');

figure;
[T, Y] = RCGAsimulate(@Model_Example_mex_2);
plot(T,Y,'-','LineWidth',2);
xlabel('Time');
ylabel('Concentration');
legend('X_1','X_2','Location','best');
title('SBML file --> MEX file');
