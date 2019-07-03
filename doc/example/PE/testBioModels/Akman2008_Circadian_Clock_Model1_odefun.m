function [output] = Akman2008_Circadian_Clock_Model1_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Akman2008_Circadian_Clock_Model1
% Generated: 03-Jul-2019 18:29:11
% 
% [output] = Akman2008_Circadian_Clock_Model1_odefun() => output = initial conditions in column vector
% [output] = Akman2008_Circadian_Clock_Model1_odefun('states') => output = state names in cell-array
% [output] = Akman2008_Circadian_Clock_Model1_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = Akman2008_Circadian_Clock_Model1_odefun('parameters') => output = parameter names in cell-array
% [output] = Akman2008_Circadian_Clock_Model1_odefun('parametervalues') => output = parameter values in column vector
% [output] = Akman2008_Circadian_Clock_Model1_odefun('variablenames') => output = variable names in cell-array
% [output] = Akman2008_Circadian_Clock_Model1_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = Akman2008_Circadian_Clock_Model1_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): FC
% statevector(2): FCp
% statevector(3): FN
% statevector(4): FNp
% statevector(5): MF
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [2.46246, 2.71231, 1.844, 2.74225, 0.725579];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'FC', 'FCp', 'FN', 'FNp', 'MF'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'vs', 'ki', 'n', 'vm', 'km', 'ks', 'vd', 'k1n', 'k2n', 'ksp', ...
			'vdp', 'k1np', 'k2np', 'amp', 'dawn', 'dusk', 'nucleus', 'cytoplasm'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [1.22363, 5.04543, 6.3958, 0.885376, 0.0846004, 0.313846, 0.161111, 0.222637, 0.331485, 0.29484, ...
			0.13975, 0.272306, 0.295421, 0, 6, 18, 1, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'Tot_FRQ'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'FC+FCp+FN+FNp'};
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
	if length(parameterValuesNew) ~= 18,
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
FC = statevector(1);
FCp = statevector(2);
FN = statevector(3);
FNp = statevector(4);
MF = statevector(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	vs = 1.22363;
	ki = 5.04543;
	n = 6.3958;
	vm = 0.885376;
	km = 0.0846004;
	ks = 0.313846;
	vd = 0.161111;
	k1n = 0.222637;
	k2n = 0.331485;
	ksp = 0.29484;
	vdp = 0.13975;
	k1np = 0.272306;
	k2np = 0.295421;
	amp = 0;
	dawn = 6;
	dusk = 18;
	nucleus = 1;
	cytoplasm = 1;
else
	vs = parameterValuesNew(1);
	ki = parameterValuesNew(2);
	n = parameterValuesNew(3);
	vm = parameterValuesNew(4);
	km = parameterValuesNew(5);
	ks = parameterValuesNew(6);
	vd = parameterValuesNew(7);
	k1n = parameterValuesNew(8);
	k2n = parameterValuesNew(9);
	ksp = parameterValuesNew(10);
	vdp = parameterValuesNew(11);
	k1np = parameterValuesNew(12);
	k2np = parameterValuesNew(13);
	amp = parameterValuesNew(14);
	dawn = parameterValuesNew(15);
	dusk = parameterValuesNew(16);
	nucleus = parameterValuesNew(17);
	cytoplasm = parameterValuesNew(18);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tot_FRQ = FC+FCp+FN+FNp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MFtrn = (vs + amp * ((1 + tanh(2 * (time - 24 * floor(time / 24) - dawn))) * (1 - tanh(2 * (time - 24 * floor(time / 24) - dusk))) / 4)) * power(ki, n) / (power(ki, n) + power(FN + FNp, n));
MFdeg = vm * MF / (km + MF);
FCtrl = ks * MF;
FCdeg = vd * FC;
FCtrs = k1n * FC - k2n * FN;
FCptrl = ksp * MF;
FCpdeg = vdp * FCp;
FCptrs = k1np * FCp - k2np * FNp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FC_dot = (+FCtrl-FCdeg-FCtrs)/cytoplasm;
FCp_dot = (+FCptrl-FCpdeg-FCptrs)/cytoplasm;
FN_dot = (+FCtrs)/nucleus;
FNp_dot = (+FCptrs)/nucleus;
MF_dot = (+MFtrn-MFdeg)/nucleus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = FC_dot;
output(2) = FCp_dot;
output(3) = FN_dot;
output(4) = FNp_dot;
output(5) = MF_dot;
% return a column vector 
output = output(:);
return


