function dSdt = hill(t,S,p)
% function [dSdt, flag, new_data] = hill(t,S,p)

% hill

S1 = S(1);
S2 = S(2);
S3 = S(3);

S4 = 0; % {isSpecie:compart:concentration}
S0 = 5; % {isSpecie:compart:concentration}
J1_Vmax = 5.5; % {isParameter}
J1_n = 4; % {isParameter}
J1_K = 0.5; % {isParameter}
J2_J2_k = 0.10000000000000001; % {isParameter}
J3_J3_k = 0.10000000000000001; % {isParameter}
J0_J0_k = 0.01; % {isParameter}
compart = 1; % {isCompartment:}

S4 = p(1); % {isSpecie:compart:concentration}
S0 = p(2); % {isSpecie:compart:concentration}
J1_Vmax = p(3); % {isParameter}
J1_n = p(4); % {isParameter}
J1_K = p(5); % {isParameter}
J2_J2_k = p(6); % {isParameter}
J3_J3_k = p(7); % {isParameter}
J0_J0_k = p(8); % {isParameter}
compart = p(9); % {isCompartment:}

J1 = J1_Vmax * power(S1, J1_n) / (power(J1_K, J1_n) + power(S1, J1_n));
J2 = J2_J2_k * S2;
J3 = J3_J3_k * S3;
J0 = J0_J0_k * S0;

dSdt = zeros(3,1);
dSdt(1) = (-J1+J0)/compart; % {isSpecie:compart:concentration}
dSdt(2) = (+J1-J2)/compart; % {isSpecie:compart:concentration}
dSdt(3) = (+J2-J3)/compart; % {isSpecie:compart:concentration}

% flag = 0;
% new_data = [];

% fprintf('%e\n',t);