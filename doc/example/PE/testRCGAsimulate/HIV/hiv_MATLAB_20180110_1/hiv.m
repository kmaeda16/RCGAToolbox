function ydot = hiv(t,y,param)

M  = y(1);
P  = y(2);
S  = y(3);
I  = y(4);
ES = y(5);
EP = y(6);
E  = y(7);
EI = y(8);
EJ = y(9);

kmd  = param(1); % 0.1   % Fixed
kdm  = param(2); % 0.001 % Fixed
kon  = param(3); % 100   % Fixed
ks   = param(4); % 46.349292 (known optimum)  % 10^-6 -- 10^+6
kcat = param(5); % 5.491365 (known optimum)   % 10^-6 -- 10^+6
kp   = param(6); % 269.804443 (known optimum) % 10^-6 -- 10^+6
ki   = param(7); % 0.000177 (known optimum)   % 10^-6 -- 10^+6
kde  = param(8); % 0.000582 (known optimum)   % 10^-6 -- 10^+6

M_dot  = - 2 * kmd * M * M + 2 * kdm * E;
P_dot  =   kcat * ES - kon * P * E + kp * EP;
S_dot  = - kon * S * E + ks * ES;
I_dot  = - kon * I * E + ki * EI;
ES_dot = kon * S * E - ks * ES - kcat * ES;
EP_dot = kon * P * E - kp * EP;
E_dot  = kmd * M * M - kdm * E - kon * S * E + ks * ES ...
    + kcat * ES - kon * P * E + kp * EP - kon * I * E + ki * EI;
EI_dot = kon * I * E - ki * EI - kde * EI;
EJ_dot = kde * EI;

ydot = zeros(9,1);
ydot(1) = M_dot;
ydot(2) = P_dot;
ydot(3) = S_dot;
ydot(4) = I_dot;
ydot(5) = ES_dot;
ydot(6) = EP_dot;
ydot(7) = E_dot;
ydot(8) = EI_dot;
ydot(9) = EJ_dot;
