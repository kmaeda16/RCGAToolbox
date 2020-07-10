function [output] = Threestep_mod_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NoName
% Generated: 08-Jul-2020 12:04:14
% 
% [output] = Threestep_mod_odefun() => output = initial conditions in column vector
% [output] = Threestep_mod_odefun('states') => output = state names in cell-array
% [output] = Threestep_mod_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = Threestep_mod_odefun('parameters') => output = parameter names in cell-array
% [output] = Threestep_mod_odefun('parametervalues') => output = parameter values in column vector
% [output] = Threestep_mod_odefun('variablenames') => output = variable names in cell-array
% [output] = Threestep_mod_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = Threestep_mod_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): E1
% statevector(2): E2
% statevector(3): E3
% statevector(4): G1
% statevector(5): G2
% statevector(6): G3
% statevector(7): M1
% statevector(8): M2
% statevector(9): P
% statevector(10): S
% statevector(11): VOID
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0.4, 0.36409, 0.29457, 0.66667, 0.57254, 0.41758, 1.419, 0.93464, 0.05, 0.1, ...
		0];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'E1', 'E2', 'E3', 'G1', 'G2', 'G3', 'M1', 'M2', 'P', 'S', ...
			'VOID'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'K4', 'K5', 'K6', 'Ka1', 'Ka2', 'Ka3', 'Ki1', 'Ki2', 'Ki3', 'Km1', ...
			'Km2', 'Km3', 'Km4', 'Km5', 'Km6', 'V1', 'V2', 'V3', 'V4', 'V5', ...
			'V6', 'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'kcat1', 'kcat2', 'kcat3', ...
			'na1', 'na2', 'na3', 'ni1', 'ni2', 'ni3', 'default0'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
			1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, ...
			0.1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, ...
			2, 2, 2, 2, 2, 2, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {};
	else
		error('Wrong input arguments! Please read the help text to the ODE file.');
	end
	output = output(:);
	return
elseif nargin == 2,
	time = varargin{1};
	statevector = varargin{2};
elseif nargin == 3,
	time = varargin{1};
	statevector = varargin{2};
	parameterValuesNew = varargin{3};
	if length(parameterValuesNew) ~= 37,
		parameterValuesNew = [];
	end
elseif nargin == 4,
	time = varargin{1};
	statevector = varargin{2};
	parameterValuesNew = varargin{4};
else
	error('Wrong input arguments! Please read the help text to the ODE file.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1 = statevector(1);
E2 = statevector(2);
E3 = statevector(3);
G1 = statevector(4);
G2 = statevector(5);
G3 = statevector(6);
M1 = statevector(7);
M2 = statevector(8);
P = statevector(9);
S = statevector(10);
VOID = statevector(11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	K4 = 1;
	K5 = 1;
	K6 = 1;
	Ka1 = 1;
	Ka2 = 1;
	Ka3 = 1;
	Ki1 = 1;
	Ki2 = 1;
	Ki3 = 1;
	Km1 = 1;
	Km2 = 1;
	Km3 = 1;
	Km4 = 1;
	Km5 = 1;
	Km6 = 1;
	V1 = 1;
	V2 = 1;
	V3 = 1;
	V4 = 0.1;
	V5 = 0.1;
	V6 = 0.1;
	k1 = 1;
	k2 = 1;
	k3 = 1;
	k4 = 0.1;
	k5 = 0.1;
	k6 = 0.1;
	kcat1 = 1;
	kcat2 = 1;
	kcat3 = 1;
	na1 = 2;
	na2 = 2;
	na3 = 2;
	ni1 = 2;
	ni2 = 2;
	ni3 = 2;
	default0 = 1;
else
	K4 = parameterValuesNew(1);
	K5 = parameterValuesNew(2);
	K6 = parameterValuesNew(3);
	Ka1 = parameterValuesNew(4);
	Ka2 = parameterValuesNew(5);
	Ka3 = parameterValuesNew(6);
	Ki1 = parameterValuesNew(7);
	Ki2 = parameterValuesNew(8);
	Ki3 = parameterValuesNew(9);
	Km1 = parameterValuesNew(10);
	Km2 = parameterValuesNew(11);
	Km3 = parameterValuesNew(12);
	Km4 = parameterValuesNew(13);
	Km5 = parameterValuesNew(14);
	Km6 = parameterValuesNew(15);
	V1 = parameterValuesNew(16);
	V2 = parameterValuesNew(17);
	V3 = parameterValuesNew(18);
	V4 = parameterValuesNew(19);
	V5 = parameterValuesNew(20);
	V6 = parameterValuesNew(21);
	k1 = parameterValuesNew(22);
	k2 = parameterValuesNew(23);
	k3 = parameterValuesNew(24);
	k4 = parameterValuesNew(25);
	k5 = parameterValuesNew(26);
	k6 = parameterValuesNew(27);
	kcat1 = parameterValuesNew(28);
	kcat2 = parameterValuesNew(29);
	kcat3 = parameterValuesNew(30);
	na1 = parameterValuesNew(31);
	na2 = parameterValuesNew(32);
	na3 = parameterValuesNew(33);
	ni1 = parameterValuesNew(34);
	ni2 = parameterValuesNew(35);
	ni3 = parameterValuesNew(36);
	default0 = parameterValuesNew(37);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vdeg1 = default0 * Function_for_vdeg1(G1, default0, k1);
vdeg2 = default0 * Function_for_vdeg2(G2, default0, k2);
vdeg3 = default0 * Function_for_vdeg3(G3, default0, k3);
vdeg4 = default0 * Function_for_vdeg4(E1, default0, k4);
vdeg5 = default0 * Function_for_vdeg5(E2, default0, k5);
vdeg6 = default0 * Function_for_vdeg6(E3, default0, k6);
vmet1 = default0 * Function_for_vmet1(E1, Km1, Km2, M1, S, default0, kcat1);
vmet2 = default0 * Function_for_vmet2(E2, Km3, Km4, M1, M2, default0, kcat2);
vmet3 = default0 * Function_for_vmet3(E3, Km5, Km6, M2, P, default0, kcat3);
vtl1 = default0 * Function_for_vtl1(G1, K4, V4, default0);
vtl2 = default0 * Function_for_vtl2(G2, K5, V5, default0);
vtl3 = default0 * Function_for_vtl3(G3, K6, V6, default0);
vts1 = default0 * Function_for_vts1(Ka1, Ki1, P, S, V1, default0, na1, ni1);
vts2 = default0 * Function_for_vts2(Ka2, Ki2, M1, P, V2, default0, na2, ni2);
vts3 = default0 * Function_for_vts3(Ka3, Ki3, M2, P, V3, default0, na3, ni3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1_dot = (-vdeg4+vtl1)/default0;
E2_dot = (-vdeg5+vtl2)/default0;
E3_dot = (-vdeg6+vtl3)/default0;
G1_dot = (-vdeg1+vts1)/default0;
G2_dot = (-vdeg2+vts2)/default0;
G3_dot = (-vdeg3+vts3)/default0;
M1_dot = (+vmet1-vmet2)/default0;
M2_dot = (+vmet2-vmet3)/default0;
P_dot = 0;
S_dot = 0;
VOID_dot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = E1_dot;
output(2) = E2_dot;
output(3) = E3_dot;
output(4) = G1_dot;
output(5) = G2_dot;
output(6) = G3_dot;
output(7) = M1_dot;
output(8) = M2_dot;
output(9) = P_dot;
output(10) = S_dot;
output(11) = VOID_dot;
% return a column vector 
output = output(:);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = Function_for_vdeg1(G1,default0,k1)
global time
result = k1*G1/default0;
return

function [result] = Function_for_vdeg5(E2,default0,k5)
global time
result = k5*E2/default0;
return

function [result] = Function_for_vdeg3(G3,default0,k3)
global time
result = k3*G3/default0;
return

function [result] = Function_for_vmet1(E1,Km1,Km2,M1,S,default0,kcat1)
global time
result = kcat1*E1*(1/Km1)*(S-M1)/(1+S/Km1+M1/Km2)/default0;
return

function [result] = Function_for_vmet2(E2,Km3,Km4,M1,M2,default0,kcat2)
global time
result = kcat2*E2*(1/Km3)*(M1-M2)/(1+M1/Km3+M2/Km4)/default0;
return

function [result] = Function_for_vmet3(E3,Km5,Km6,M2,P,default0,kcat3)
global time
result = kcat3*E3*(1/Km5)*(M2-P)/(1+M2/Km5+P/Km6)/default0;
return

function [result] = Function_for_vtl1(G1,K4,V4,default0)
global time
result = V4*G1/(K4+G1)/default0;
return

function [result] = Function_for_vtl2(G2,K5,V5,default0)
global time
result = V5*G2/(K5+G2)/default0;
return

function [result] = Function_for_vdeg2(G2,default0,k2)
global time
result = k2*G2/default0;
return

function [result] = Function_for_vtl3(G3,K6,V6,default0)
global time
result = V6*G3/(K6+G3)/default0;
return

function [result] = Function_for_vdeg4(E1,default0,k4)
global time
result = k4*E1/default0;
return

function [result] = Function_for_vts1(Ka1,Ki1,P,S,V1,default0,na1,ni1)
global time
result = V1/(1+(P/Ki1)^(ni1)+(Ka1/S)^(na1))/default0;
return

function [result] = Function_for_vdeg6(E3,default0,k6)
global time
result = k6*E3/default0;
return

function [result] = Function_for_vts3(Ka3,Ki3,M2,P,V3,default0,na3,ni3)
global time
result = V3/(1+(P/Ki3)^(ni3)+(Ka3/M2)^(na3))/default0;
return

function [result] = Function_for_vts2(Ka2,Ki2,M1,P,V2,default0,na2,ni2)
global time
result = V2/(1+(P/Ki2)^(ni2)+(Ka2/M1)^(na2))/default0;
return


