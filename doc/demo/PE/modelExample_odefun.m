function [output] = modelExample_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelExample
% Generated: 06-Apr-2020 16:08:47
% 
% [output] = modelExample_odefun() => output = initial conditions in column vector
% [output] = modelExample_odefun('states') => output = state names in cell-array
% [output] = modelExample_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = modelExample_odefun('parameters') => output = parameter names in cell-array
% [output] = modelExample_odefun('parametervalues') => output = parameter values in column vector
% [output] = modelExample_odefun('variablenames') => output = variable names in cell-array
% [output] = modelExample_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = modelExample_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): S1
% statevector(2): S2
% statevector(3): S3
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0, 0, 0];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'S1', 'S2', 'S3'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'S4', 'S0', 'J1_Vmax', 'J1_n', 'J1_K', 'J2_J2_k', 'J3_J3_k', 'J0_J0_k', 'compart'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [0, 5, 5.5, 4, 0.5, 0.1, 0.1, 0.01, 1];
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
	if length(parameterValuesNew) ~= 9,
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
S1 = statevector(1);
S2 = statevector(2);
S3 = statevector(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	S4 = 0;
	S0 = 5;
	J1_Vmax = 5.5;
	J1_n = 4;
	J1_K = 0.5;
	J2_J2_k = 0.1;
	J3_J3_k = 0.1;
	J0_J0_k = 0.01;
	compart = 1;
else
	S4 = parameterValuesNew(1);
	S0 = parameterValuesNew(2);
	J1_Vmax = parameterValuesNew(3);
	J1_n = parameterValuesNew(4);
	J1_K = parameterValuesNew(5);
	J2_J2_k = parameterValuesNew(6);
	J3_J3_k = parameterValuesNew(7);
	J0_J0_k = parameterValuesNew(8);
	compart = parameterValuesNew(9);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J1 = J1_Vmax * power(S1, J1_n) / (power(J1_K, J1_n) + power(S1, J1_n));
J2 = J2_J2_k * S2;
J3 = J3_J3_k * S3;
J0 = J0_J0_k * S0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S1_dot = (-J1+J0)/compart;
S2_dot = (+J1-J2)/compart;
S3_dot = (+J2-J3)/compart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = S1_dot;
output(2) = S2_dot;
output(3) = S3_dot;
% return a column vector 
output = output(:);
return


