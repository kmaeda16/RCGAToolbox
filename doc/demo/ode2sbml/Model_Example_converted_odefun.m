function [output] = Model_Example_converted_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model_Example
% Generated: 22-Oct-2020 10:21:16
% 
% [output] = Model_Example_converted_odefun() => output = initial conditions in column vector
% [output] = Model_Example_converted_odefun('states') => output = state names in cell-array
% [output] = Model_Example_converted_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = Model_Example_converted_odefun('parameters') => output = parameter names in cell-array
% [output] = Model_Example_converted_odefun('parametervalues') => output = parameter values in column vector
% [output] = Model_Example_converted_odefun('variablenames') => output = variable names in cell-array
% [output] = Model_Example_converted_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = Model_Example_converted_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): X1
% statevector(2): X2
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0, 0];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'X1', 'X2'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'X0', 'k1', 'k2', 'k3', 'K2', 'K3', 'rootCompartment'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [0.1, 1, 1, 1, 1, 1, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'X12'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'(X1/rootCompartment)+(X2/rootCompartment)'};
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
	if length(parameterValuesNew) ~= 7,
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
X1 = statevector(1);
X2 = statevector(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	X0 = 0.1;
	k1 = 1;
	k2 = 1;
	k3 = 1;
	K2 = 1;
	K3 = 1;
	rootCompartment = 1;
else
	X0 = parameterValuesNew(1);
	k1 = parameterValuesNew(2);
	k2 = parameterValuesNew(3);
	k3 = parameterValuesNew(4);
	K2 = parameterValuesNew(5);
	K3 = parameterValuesNew(6);
	rootCompartment = parameterValuesNew(7);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X12 = (X1/rootCompartment)+(X2/rootCompartment);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1 = k1 * X0;
v2 = k2 * (X1/rootCompartment) / (K2 + (X1/rootCompartment));
v3 = k3 * (X2/rootCompartment) / (K3 + (X2/rootCompartment));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1_dot = +v1-v2;
X2_dot = +v2-v3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = X1_dot;
output(2) = X2_dot;
% return a column vector 
output = output(:);
return


