function [output] = BIOMD0000000016_url_mod_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goldbeter1995_CircClock
% Generated: 05-Jul-2019 16:25:40
% 
% [output] = BIOMD0000000016_url_mod_odefun() => output = initial conditions in column vector
% [output] = BIOMD0000000016_url_mod_odefun('states') => output = state names in cell-array
% [output] = BIOMD0000000016_url_mod_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = BIOMD0000000016_url_mod_odefun('parameters') => output = parameter names in cell-array
% [output] = BIOMD0000000016_url_mod_odefun('parametervalues') => output = parameter values in column vector
% [output] = BIOMD0000000016_url_mod_odefun('variablenames') => output = variable names in cell-array
% [output] = BIOMD0000000016_url_mod_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = BIOMD0000000016_url_mod_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): EmptySet
% statevector(2): M
% statevector(3): P0
% statevector(4): P1
% statevector(5): P2
% statevector(6): Pn
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0, 0.1, 0.25, 0.25, 0.25, 0.25];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'EmptySet', 'M', 'P0', 'P1', 'P2', 'Pn'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'rM_Vs', 'rM_KI', 'rM_n', 'rTL_ks', 'rP01_V1', 'rP01_K1', 'rP10_V2', 'rP10_K2', 'rP12_V3', 'rP12_K3', ...
			'rP21_V4', 'rP21_K4', 'rP2n_k1', 'rPn2_k2', 'rmRNAd_Km', 'rmRNAd_Vm', 'rVd_Vd', 'rVd_Kd', 'default0', 'CYTOPLASM', ...
			'compartment_0000004'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [0.76, 1, 4, 0.38, 3.2, 2, 1.58, 2, 5, 2, ...
			2.5, 2, 1.9, 1.3, 0.5, 0.65, 0.95, 0.2, 1e-15, 1e-15, ...
			1e-15];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'Pt'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'P0+P1+P2+Pn'};
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
	if length(parameterValuesNew) ~= 21,
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
EmptySet = statevector(1);
M = statevector(2);
P0 = statevector(3);
P1 = statevector(4);
P2 = statevector(5);
Pn = statevector(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	rM_Vs = 0.76;
	rM_KI = 1;
	rM_n = 4;
	rTL_ks = 0.38;
	rP01_V1 = 3.2;
	rP01_K1 = 2;
	rP10_V2 = 1.58;
	rP10_K2 = 2;
	rP12_V3 = 5;
	rP12_K3 = 2;
	rP21_V4 = 2.5;
	rP21_K4 = 2;
	rP2n_k1 = 1.9;
	rPn2_k2 = 1.3;
	rmRNAd_Km = 0.5;
	rmRNAd_Vm = 0.65;
	rVd_Vd = 0.95;
	rVd_Kd = 0.2;
	default0 = 1e-15;
	CYTOPLASM = 1e-15;
	compartment_0000004 = 1e-15;
else
	rM_Vs = parameterValuesNew(1);
	rM_KI = parameterValuesNew(2);
	rM_n = parameterValuesNew(3);
	rTL_ks = parameterValuesNew(4);
	rP01_V1 = parameterValuesNew(5);
	rP01_K1 = parameterValuesNew(6);
	rP10_V2 = parameterValuesNew(7);
	rP10_K2 = parameterValuesNew(8);
	rP12_V3 = parameterValuesNew(9);
	rP12_K3 = parameterValuesNew(10);
	rP21_V4 = parameterValuesNew(11);
	rP21_K4 = parameterValuesNew(12);
	rP2n_k1 = parameterValuesNew(13);
	rPn2_k2 = parameterValuesNew(14);
	rmRNAd_Km = parameterValuesNew(15);
	rmRNAd_Vm = parameterValuesNew(16);
	rVd_Vd = parameterValuesNew(17);
	rVd_Kd = parameterValuesNew(18);
	default0 = parameterValuesNew(19);
	CYTOPLASM = parameterValuesNew(20);
	compartment_0000004 = parameterValuesNew(21);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pt = P0+P1+P2+Pn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rM = default0 * rM_Vs * power(rM_KI, rM_n) / (power(rM_KI, rM_n) + power(Pn, rM_n));
rTL = rTL_ks * M * default0;
rP01 = CYTOPLASM * rP01_V1 * P0 / (rP01_K1 + P0);
rP10 = CYTOPLASM * rP10_V2 * P1 / (rP10_K2 + P1);
rP12 = CYTOPLASM * rP12_V3 * P1 / (rP12_K3 + P1);
rP21 = CYTOPLASM * rP21_V4 * P2 / (rP21_K4 + P2);
rP2n = rP2n_k1 * P2 * CYTOPLASM;
rPn2 = rPn2_k2 * Pn * compartment_0000004;
rmRNAd = rmRNAd_Vm * M * CYTOPLASM / (rmRNAd_Km + M);
rVd = CYTOPLASM * rVd_Vd * P2 / (rVd_Kd + P2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EmptySet_dot = 0;
M_dot = (+rM-rmRNAd)/CYTOPLASM;
P0_dot = (+rTL-rP01+rP10)/CYTOPLASM;
P1_dot = (+rP01-rP10-rP12+rP21)/CYTOPLASM;
P2_dot = (+rP12-rP21-rP2n+rPn2-rVd)/CYTOPLASM;
Pn_dot = (+rP2n-rPn2)/compartment_0000004;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = EmptySet_dot;
output(2) = M_dot;
output(3) = P0_dot;
output(4) = P1_dot;
output(5) = P2_dot;
output(6) = Pn_dot;
% return a column vector 
output = output(:);
return


