function [output] = HIV_mod_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NoName
% Generated: 24-Mar-2020 15:19:25
% 
% [output] = HIV_mod_odefun() => output = initial conditions in column vector
% [output] = HIV_mod_odefun('states') => output = state names in cell-array
% [output] = HIV_mod_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = HIV_mod_odefun('parameters') => output = parameter names in cell-array
% [output] = HIV_mod_odefun('parametervalues') => output = parameter values in column vector
% [output] = HIV_mod_odefun('variablenames') => output = variable names in cell-array
% [output] = HIV_mod_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = HIV_mod_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): E
% statevector(2): EI
% statevector(3): EJ
% statevector(4): EP
% statevector(5): ES
% statevector(6): I
% statevector(7): M
% statevector(8): P
% statevector(9): S
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0.005387, 0, 0, 0, 0, 0, 0, 0, 24.6378];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'E', 'EI', 'EJ', 'EP', 'ES', 'I', 'M', 'P', 'S'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'kcat', 'kde', 'kdm', 'ki', 'kmd', 'kon', 'kp', 'ks', 'default0'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [5.49137, 0.000582, 0, 0.000177, 0.1, 100, 269.804, 46.3493, 1];
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
E = statevector(1);
EI = statevector(2);
EJ = statevector(3);
EP = statevector(4);
ES = statevector(5);
I = statevector(6);
M = statevector(7);
P = statevector(8);
S = statevector(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	kcat = 5.49137;
	kde = 0.000582;
	kdm = 0;
	ki = 0.000177;
	kmd = 0.1;
	kon = 100;
	kp = 269.804;
	ks = 46.3493;
	default0 = 1;
else
	kcat = parameterValuesNew(1);
	kde = parameterValuesNew(2);
	kdm = parameterValuesNew(3);
	ki = parameterValuesNew(4);
	kmd = parameterValuesNew(5);
	kon = parameterValuesNew(6);
	kp = parameterValuesNew(7);
	ks = parameterValuesNew(8);
	default0 = parameterValuesNew(9);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ved = default0 * Function_for_ved(E, default0, kdm);
vef = default0 * Function_for_vef(M, default0, kmd);
veid = default0 * Function_for_veid(EI, default0, ki);
veif = default0 * Function_for_veif(E, I, default0, kon);
vejf = default0 * Function_for_vejf(EI, default0, kde);
vepd = default0 * Function_for_vepd(EP, default0, kp);
vepf = default0 * Function_for_vepf(E, P, default0, kon);
vesd = default0 * Function_for_vesd(ES, default0, ks);
vesf = default0 * Function_for_vesf(E, S, default0, kon);
vpf = default0 * Function_for_vpf(ES, default0, kcat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_dot = (-ved+vef+veid-veif+vepd-vepf+vesd-vesf+vpf)/default0;
EI_dot = (-veid+veif-vejf)/default0;
EJ_dot = (+vejf)/default0;
EP_dot = (-vepd+vepf)/default0;
ES_dot = (-vesd+vesf-vpf)/default0;
I_dot = (+veid-veif)/default0;
M_dot = (+2*ved-2*vef)/default0;
P_dot = (+vepd-vepf+vpf)/default0;
S_dot = (+vesd-vesf)/default0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = E_dot;
output(2) = EI_dot;
output(3) = EJ_dot;
output(4) = EP_dot;
output(5) = ES_dot;
output(6) = I_dot;
output(7) = M_dot;
output(8) = P_dot;
output(9) = S_dot;
% return a column vector 
output = output(:);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = Function_for_vesd(ES,default0,ks)
global time
result = ks*ES/default0;
return

function [result] = Function_for_veif(E,I,default0,kon)
global time
result = kon*I*E/default0;
return

function [result] = Function_for_ved(E,default0,kdm)
global time
result = kdm*E/default0;
return

function [result] = Function_for_vepf(E,P,default0,kon)
global time
result = kon*P*E/default0;
return

function [result] = Function_for_vepd(EP,default0,kp)
global time
result = kp*EP/default0;
return

function [result] = Function_for_vejf(EI,default0,kde)
global time
result = kde*EI/default0;
return

function [result] = Function_for_vesf(E,S,default0,kon)
global time
result = kon*S*E/default0;
return

function [result] = Function_for_veid(EI,default0,ki)
global time
result = ki*EI/default0;
return

function [result] = Function_for_vpf(ES,default0,kcat)
global time
result = kcat*ES/default0;
return

function [result] = Function_for_vef(M,default0,kmd)
global time
result = kmd*M*M/default0;
return


