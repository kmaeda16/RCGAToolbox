function [output] = BIOMD0000000005_url_mod_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tyson1991 - Cell Cycle 6 var
% Generated: 05-Jul-2019 16:25:42
% 
% [output] = BIOMD0000000005_url_mod_odefun() => output = initial conditions in column vector
% [output] = BIOMD0000000005_url_mod_odefun('states') => output = state names in cell-array
% [output] = BIOMD0000000005_url_mod_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = BIOMD0000000005_url_mod_odefun('parameters') => output = parameter names in cell-array
% [output] = BIOMD0000000005_url_mod_odefun('parametervalues') => output = parameter values in column vector
% [output] = BIOMD0000000005_url_mod_odefun('variablenames') => output = variable names in cell-array
% [output] = BIOMD0000000005_url_mod_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = BIOMD0000000005_url_mod_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): C2
% statevector(2): CP
% statevector(3): M
% statevector(4): pM
% statevector(5): Y
% statevector(6): YP
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0, 0.75, 0, 0.25, 0, 0];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'C2', 'CP', 'M', 'pM', 'Y', 'YP'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'EmptySet', 'Reaction1_k6', 'Reaction2_k8notP', 'Reaction3_k9', 'Reaction4_k3', 'Reaction5_k5notP', 'Reaction6_k1aa', 'Reaction7_k2', 'Reaction8_k7', 'Reaction9_k4', ...
			'Reaction9_k4prime', 'cell'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [0, 1, 1e+06, 1000, 200, 0, 0.015, 0, 0.6, 180, ...
			0.018, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'YT', 'CT'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'(Y/cell)+(YP/cell)+(M/cell)+(pM/cell)', '(C2/cell)+(CP/cell)+(M/cell)+(pM/cell)'};
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
	if length(parameterValuesNew) ~= 12,
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
C2 = statevector(1);
CP = statevector(2);
M = statevector(3);
pM = statevector(4);
Y = statevector(5);
YP = statevector(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	EmptySet = 0;
	Reaction1_k6 = 1;
	Reaction2_k8notP = 1e+06;
	Reaction3_k9 = 1000;
	Reaction4_k3 = 200;
	Reaction5_k5notP = 0;
	Reaction6_k1aa = 0.015;
	Reaction7_k2 = 0;
	Reaction8_k7 = 0.6;
	Reaction9_k4 = 180;
	Reaction9_k4prime = 0.018;
	cell = 1;
else
	EmptySet = parameterValuesNew(1);
	Reaction1_k6 = parameterValuesNew(2);
	Reaction2_k8notP = parameterValuesNew(3);
	Reaction3_k9 = parameterValuesNew(4);
	Reaction4_k3 = parameterValuesNew(5);
	Reaction5_k5notP = parameterValuesNew(6);
	Reaction6_k1aa = parameterValuesNew(7);
	Reaction7_k2 = parameterValuesNew(8);
	Reaction8_k7 = parameterValuesNew(9);
	Reaction9_k4 = parameterValuesNew(10);
	Reaction9_k4prime = parameterValuesNew(11);
	cell = parameterValuesNew(12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YT = (Y/cell)+(YP/cell)+(M/cell)+(pM/cell);
CT = (C2/cell)+(CP/cell)+(M/cell)+(pM/cell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reaction1 = cell * Reaction1_k6 * (M/cell);
Reaction2 = cell * (C2/cell) * Reaction2_k8notP;
Reaction3 = cell * (CP/cell) * Reaction3_k9;
Reaction4 = cell * (CP/cell) * Reaction4_k3 * (Y/cell);
Reaction5 = cell * Reaction5_k5notP * (M/cell);
Reaction6 = cell * Reaction6_k1aa;
Reaction7 = cell * Reaction7_k2 * (Y/cell);
Reaction8 = cell * Reaction8_k7 * (YP/cell);
Reaction9 = cell * (pM/cell) * (Reaction9_k4prime + Reaction9_k4 * power((M/cell) / CT, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C2_dot = +Reaction1-Reaction2+Reaction3;
CP_dot = +Reaction2-Reaction3-Reaction4;
M_dot = -Reaction1-Reaction5+Reaction9;
pM_dot = +Reaction4+Reaction5-Reaction9;
Y_dot = -Reaction4+Reaction6-Reaction7;
YP_dot = +Reaction1-Reaction8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = C2_dot;
output(2) = CP_dot;
output(3) = M_dot;
output(4) = pM_dot;
output(5) = Y_dot;
output(6) = YP_dot;
% return a column vector 
output = output(:);
return


