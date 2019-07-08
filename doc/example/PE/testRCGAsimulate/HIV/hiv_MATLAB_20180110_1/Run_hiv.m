tic

T0 = 0;
Tend = 3600;
Intvl = 0.4;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%% HIV model
drivsfun = @hiv;

kmd  = 0.1;   % Fixed
kdm  = 0.001; % Fixed
kon  = 100;   % Fixed
ks   = 46.349292; % (known optimum) % 10^-6 -- 10^+6
kcat = 5.491365;  % (known optimum) % 10^-6 -- 10^+6
kp   = 269.804443;% (known optimum) % 10^-6 -- 10^+6
ki   = 0.000177;  % (known optimum) % 10^-6 -- 10^+6
kde  = 0.000582;  % (known optimum) % 10^-6 -- 10^+6

ModelInput = [ 0      24.637840 0.005387 -0.004763 ]; % Experiment 1
% ModelInput = [ 0.0015 23.456802 0.005183 -0.004950 ]; % Experiment 2
% ModelInput = [ 0.003  27.159763 0.006000 -0.017078 ]; % Experiment 3
% ModelInput = [ 0.004  16.190568 0.004119 -0.007473 ]; % Experiment 4
% ModelInput = [ 0.004  24.672660 0.003051  0.002483 ]; % Experiment 5

M  = 0;
P  = 0;
S  = ModelInput(2);
I  = ModelInput(1);
ES = 0;
EP = 0;
E  = ModelInput(3);
EI = 0;
EJ = 0;

y0 = [ M; P; S; I; ES; EP; E; EI; EJ ];

%%
param = [ kmd, kdm, kon, ks, kcat, kp, ki, kde ]; % 1 - 8

f = @(t,y)drivsfun(t,y,param);
[T, Y] = ode15s(f,T0:Intvl:Tend,y0,options);

figure;
plot(T,Y);
legend('M','P','S','I','ES','E','EI','EJ');
xlabel('Time (s^{-1})');
ylabel('Concentration (\muM)');
% ylim([0 1]);

toc
