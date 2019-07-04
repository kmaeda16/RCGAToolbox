function [output] = Bruggeman2005_AmmoniumAssimilation_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bruggeman2005_AmmoniumAssimilation
% Generated: 04-Jul-2019 17:15:15
% 
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun() => output = initial conditions in column vector
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun('states') => output = state names in cell-array
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun('parameters') => output = parameter names in cell-array
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun('parametervalues') => output = parameter values in column vector
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun('variablenames') => output = variable names in cell-array
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = Bruggeman2005_AmmoniumAssimilation_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): PII
% statevector(2): UTP
% statevector(3): PIIUMP
% statevector(4): PPi
% statevector(5): GLN
% statevector(6): PIIUMP2
% statevector(7): PIIUMP3
% statevector(8): UMP
% statevector(9): GS
% statevector(10): AMP
% statevector(11): NH4
% statevector(12): KG
% statevector(13): NADPH
% statevector(14): GLU
% statevector(15): NADP
% statevector(16): AZGLU
% statevector(17): ATP
% statevector(18): ADP
% statevector(19): AZglu
% statevector(20): AZGLN
% statevector(21): AZgln
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0.003, 0.5, 0, 0.05, 1, 0, 0, 0.01, 0.014, 0, ...
		0.05, 0.2, 0.15, 1, 0.05, 1, 2.685, 2.685, 0.1, 1, ...
		0.1];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'PII', 'UTP', 'PIIUMP', 'PPi', 'GLN', 'PIIUMP2', 'PIIUMP3', 'UMP', 'GS', 'AMP', ...
			'NH4', 'KG', 'NADPH', 'GLU', 'NADP', 'AZGLU', 'ATP', 'ADP', 'AZglu', 'AZGLN', ...
			'AZgln'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'P_i', 'UT', 'kcatut', 'Kglnut', 'Kutipii', 'Kutpii', 'Kutpiiump', 'Kututp', 'Kutippi', 'UR', ...
			'kcatur', 'Kurpiiump', 'Kurump', 'Kglnur', 'a1', 'b1', 'c1', 'd1', 'Vad', 'Kadpiikg', ...
			'Kadgln', 'Kadgs', 'e1', 'f1', 'g1', 'h1', 'i1', 'j1', 'k1', 'l1', ...
			'm1', 'n1', 'o1', 'Vdead', 'Kdeadpiikg', 'Kdeadgln', 'Kdeadpiiu', 'Kdeadgsa', 'Vgdh', 'Kgdhkg', ...
			'Kgdhnh', 'Kgdhglu', 'Kgdhnadph', 'Kgdhnadp', 'Keqgdh', 'Kgdhazglu', 'Vgog', 'Kgoggln', 'Kgogkg', 'Kgognadph', ...
			'Kgogglu', 'Kgognadp', 'Kgogaz', 'Vgs', 'aamp', 'bamp', 'camp', 'damp', 'n1amp', 'n2amp', ...
			'Kgseq', 'Kgsatp', 'Kgsglu', 'Kgsnh', 'Kgsadp', 'Kgspi', 'Kgsgln', 'Keq', 'Vgludem', 'Kgludemglu', ...
			'Kgludemeq', 'Kgludemazglu', 'Vglndem', 'Kglndemgln', 'Kglndemeq', 'Kglndemazgln', 'Vazglndem', 'Kazglndemazgln', 'Kazglndemeq', 'Kazglndemazinter', ...
			'Vazgludem', 'Kazgludemazglu', 'Kazgludemeq', 'Kazgludemazinter', 'Vadp', 'Kadp', 'ATPtot', 'GStot', 'PIItot', 'Kd1', ...
			'Kd2', 'Kd3', 'Kd1piiump', 'Kd2piiump', 'Kd3piiump', 'compartment'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [10, 0.0006, 137, 0.07, 0.0018, 0.003, 0.0035, 0.04, 0.1135, 0.0006, ...
			5.5, 0.0023, 8.4, 0.07, 1e-22, 0.5166, 0.5974, 0.0387, 0.5, 1.052e-05, ...
			0.9714, 0.001703, 1e-22, 2.766, 3.323, 0.2148, 1e-22, 1e-22, 1e-22, 0.02316, ...
			0.8821, 8.491, 0.8791, 0.5, 2.274e-06, 0.04444, 1.805e-05, 0.0002015, 360, 0.32, ...
			1.1, 10, 0.04, 0.042, 1290, 2.5, 85, 0.175, 0.007, 0.0015, ...
			11, 0.0037, 0.65, 600, 10, 2.3667, 0.1012, 10.8688, 1.1456, 19.2166, ...
			460, 0.35, 4.1, 0.1, 0.0585, 3.7, 5.65, 460, 120, 8, ...
			1e+10, 0.5, 70, 2, 1e+10, 0.25, 20, 1, 1e+10, 0.5, ...
			30, 0.3, 1e+10, 0.5, 100, 0.5, 5.37, 0.014, 0.003, 0.005, ...
			0.15, 0.15, 0.025, 0.15, 0.15, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'vAPP_GS', 'nAMP', 'PIIKG1', 'PIIUMP3KG3'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'aamp*camp/((1+power(12,n1amp)*power(AMP/(bamp*GStot),n1amp))*(1+power(12,n2amp)*power(AMP/(damp*GStot),n2amp)))*Vgs', '12*AMP/GStot', '3*PII*KG/Kd1/(1+3*KG/Kd1+3*power(KG,2)/(Kd1*Kd2)+power(KG,3)/(Kd1*Kd2*Kd3))', 'PIIUMP3*power(KG,3)/(Kd1piiump*Kd2piiump*Kd3piiump)/(1+3*KG/Kd1piiump+3*power(KG,2)/(Kd1piiump*Kd2piiump)+power(KG,3)/(Kd1piiump*Kd2piiump*Kd3piiump))'};
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
	if length(parameterValuesNew) ~= 96,
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
PII = statevector(1);
UTP = statevector(2);
PIIUMP = statevector(3);
PPi = statevector(4);
GLN = statevector(5);
PIIUMP2 = statevector(6);
PIIUMP3 = statevector(7);
UMP = statevector(8);
GS = statevector(9);
AMP = statevector(10);
NH4 = statevector(11);
KG = statevector(12);
NADPH = statevector(13);
GLU = statevector(14);
NADP = statevector(15);
AZGLU = statevector(16);
ATP = statevector(17);
ADP = statevector(18);
AZglu = statevector(19);
AZGLN = statevector(20);
AZgln = statevector(21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	P_i = 10;
	UT = 0.0006;
	kcatut = 137;
	Kglnut = 0.07;
	Kutipii = 0.0018;
	Kutpii = 0.003;
	Kutpiiump = 0.0035;
	Kututp = 0.04;
	Kutippi = 0.1135;
	UR = 0.0006;
	kcatur = 5.5;
	Kurpiiump = 0.0023;
	Kurump = 8.4;
	Kglnur = 0.07;
	a1 = 1e-22;
	b1 = 0.5166;
	c1 = 0.5974;
	d1 = 0.0387;
	Vad = 0.5;
	Kadpiikg = 1.052e-05;
	Kadgln = 0.9714;
	Kadgs = 0.001703;
	e1 = 1e-22;
	f1 = 2.766;
	g1 = 3.323;
	h1 = 0.2148;
	i1 = 1e-22;
	j1 = 1e-22;
	k1 = 1e-22;
	l1 = 0.02316;
	m1 = 0.8821;
	n1 = 8.491;
	o1 = 0.8791;
	Vdead = 0.5;
	Kdeadpiikg = 2.274e-06;
	Kdeadgln = 0.04444;
	Kdeadpiiu = 1.805e-05;
	Kdeadgsa = 0.0002015;
	Vgdh = 360;
	Kgdhkg = 0.32;
	Kgdhnh = 1.1;
	Kgdhglu = 10;
	Kgdhnadph = 0.04;
	Kgdhnadp = 0.042;
	Keqgdh = 1290;
	Kgdhazglu = 2.5;
	Vgog = 85;
	Kgoggln = 0.175;
	Kgogkg = 0.007;
	Kgognadph = 0.0015;
	Kgogglu = 11;
	Kgognadp = 0.0037;
	Kgogaz = 0.65;
	Vgs = 600;
	aamp = 10;
	bamp = 2.3667;
	camp = 0.1012;
	damp = 10.8688;
	n1amp = 1.1456;
	n2amp = 19.2166;
	Kgseq = 460;
	Kgsatp = 0.35;
	Kgsglu = 4.1;
	Kgsnh = 0.1;
	Kgsadp = 0.0585;
	Kgspi = 3.7;
	Kgsgln = 5.65;
	Keq = 460;
	Vgludem = 120;
	Kgludemglu = 8;
	Kgludemeq = 1e+10;
	Kgludemazglu = 0.5;
	Vglndem = 70;
	Kglndemgln = 2;
	Kglndemeq = 1e+10;
	Kglndemazgln = 0.25;
	Vazglndem = 20;
	Kazglndemazgln = 1;
	Kazglndemeq = 1e+10;
	Kazglndemazinter = 0.5;
	Vazgludem = 30;
	Kazgludemazglu = 0.3;
	Kazgludemeq = 1e+10;
	Kazgludemazinter = 0.5;
	Vadp = 100;
	Kadp = 0.5;
	ATPtot = 5.37;
	GStot = 0.014;
	PIItot = 0.003;
	Kd1 = 0.005;
	Kd2 = 0.15;
	Kd3 = 0.15;
	Kd1piiump = 0.025;
	Kd2piiump = 0.15;
	Kd3piiump = 0.15;
	compartment = 1;
else
	P_i = parameterValuesNew(1);
	UT = parameterValuesNew(2);
	kcatut = parameterValuesNew(3);
	Kglnut = parameterValuesNew(4);
	Kutipii = parameterValuesNew(5);
	Kutpii = parameterValuesNew(6);
	Kutpiiump = parameterValuesNew(7);
	Kututp = parameterValuesNew(8);
	Kutippi = parameterValuesNew(9);
	UR = parameterValuesNew(10);
	kcatur = parameterValuesNew(11);
	Kurpiiump = parameterValuesNew(12);
	Kurump = parameterValuesNew(13);
	Kglnur = parameterValuesNew(14);
	a1 = parameterValuesNew(15);
	b1 = parameterValuesNew(16);
	c1 = parameterValuesNew(17);
	d1 = parameterValuesNew(18);
	Vad = parameterValuesNew(19);
	Kadpiikg = parameterValuesNew(20);
	Kadgln = parameterValuesNew(21);
	Kadgs = parameterValuesNew(22);
	e1 = parameterValuesNew(23);
	f1 = parameterValuesNew(24);
	g1 = parameterValuesNew(25);
	h1 = parameterValuesNew(26);
	i1 = parameterValuesNew(27);
	j1 = parameterValuesNew(28);
	k1 = parameterValuesNew(29);
	l1 = parameterValuesNew(30);
	m1 = parameterValuesNew(31);
	n1 = parameterValuesNew(32);
	o1 = parameterValuesNew(33);
	Vdead = parameterValuesNew(34);
	Kdeadpiikg = parameterValuesNew(35);
	Kdeadgln = parameterValuesNew(36);
	Kdeadpiiu = parameterValuesNew(37);
	Kdeadgsa = parameterValuesNew(38);
	Vgdh = parameterValuesNew(39);
	Kgdhkg = parameterValuesNew(40);
	Kgdhnh = parameterValuesNew(41);
	Kgdhglu = parameterValuesNew(42);
	Kgdhnadph = parameterValuesNew(43);
	Kgdhnadp = parameterValuesNew(44);
	Keqgdh = parameterValuesNew(45);
	Kgdhazglu = parameterValuesNew(46);
	Vgog = parameterValuesNew(47);
	Kgoggln = parameterValuesNew(48);
	Kgogkg = parameterValuesNew(49);
	Kgognadph = parameterValuesNew(50);
	Kgogglu = parameterValuesNew(51);
	Kgognadp = parameterValuesNew(52);
	Kgogaz = parameterValuesNew(53);
	Vgs = parameterValuesNew(54);
	aamp = parameterValuesNew(55);
	bamp = parameterValuesNew(56);
	camp = parameterValuesNew(57);
	damp = parameterValuesNew(58);
	n1amp = parameterValuesNew(59);
	n2amp = parameterValuesNew(60);
	Kgseq = parameterValuesNew(61);
	Kgsatp = parameterValuesNew(62);
	Kgsglu = parameterValuesNew(63);
	Kgsnh = parameterValuesNew(64);
	Kgsadp = parameterValuesNew(65);
	Kgspi = parameterValuesNew(66);
	Kgsgln = parameterValuesNew(67);
	Keq = parameterValuesNew(68);
	Vgludem = parameterValuesNew(69);
	Kgludemglu = parameterValuesNew(70);
	Kgludemeq = parameterValuesNew(71);
	Kgludemazglu = parameterValuesNew(72);
	Vglndem = parameterValuesNew(73);
	Kglndemgln = parameterValuesNew(74);
	Kglndemeq = parameterValuesNew(75);
	Kglndemazgln = parameterValuesNew(76);
	Vazglndem = parameterValuesNew(77);
	Kazglndemazgln = parameterValuesNew(78);
	Kazglndemeq = parameterValuesNew(79);
	Kazglndemazinter = parameterValuesNew(80);
	Vazgludem = parameterValuesNew(81);
	Kazgludemazglu = parameterValuesNew(82);
	Kazgludemeq = parameterValuesNew(83);
	Kazgludemazinter = parameterValuesNew(84);
	Vadp = parameterValuesNew(85);
	Kadp = parameterValuesNew(86);
	ATPtot = parameterValuesNew(87);
	GStot = parameterValuesNew(88);
	PIItot = parameterValuesNew(89);
	Kd1 = parameterValuesNew(90);
	Kd2 = parameterValuesNew(91);
	Kd3 = parameterValuesNew(92);
	Kd1piiump = parameterValuesNew(93);
	Kd2piiump = parameterValuesNew(94);
	Kd3piiump = parameterValuesNew(95);
	compartment = parameterValuesNew(96);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vAPP_GS = aamp*camp/((1+power(12,n1amp)*power(AMP/(bamp*GStot),n1amp))*(1+power(12,n2amp)*power(AMP/(damp*GStot),n2amp)))*Vgs;
nAMP = 12*AMP/GStot;
PIIKG1 = 3*PII*KG/Kd1/(1+3*KG/Kd1+3*power(KG,2)/(Kd1*Kd2)+power(KG,3)/(Kd1*Kd2*Kd3));
PIIUMP3KG3 = PIIUMP3*power(KG,3)/(Kd1piiump*Kd2piiump*Kd3piiump)/(1+3*KG/Kd1piiump+3*power(KG,2)/(Kd1piiump*Kd2piiump)+power(KG,3)/(Kd1piiump*Kd2piiump*Kd3piiump));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vut1 = compartment * (kcatut * UT * UTP * PII / (Kutipii * Kututp * (1 + GLN / Kglnut) * (1 + UTP / Kututp + (PII + PIIUMP + PIIUMP2) / Kutipii + UTP * (PII + PIIUMP + PIIUMP2) / (Kutipii * Kututp) + PPi * UTP * (PII + PIIUMP + PIIUMP2) / (Kutipii * Kutippi * Kututp) + Kutpii * (PIIUMP + PIIUMP2 + PIIUMP3) / (Kutipii * Kutpiiump))));
vur1 = compartment * (kcatur * UR * PIIUMP / (Kurpiiump * (1 + Kglnur / GLN) * (1 + (1 + UMP / Kurump) * (PIIUMP + PIIUMP2 + PIIUMP3) / Kurpiiump)));
vut2 = compartment * (kcatut * UT * UTP * PIIUMP / (Kutipii * Kututp * (1 + GLN / Kglnut) * (1 + UTP / Kututp + (PII + PIIUMP + PIIUMP2) / Kutipii + UTP * (PII + PIIUMP + PIIUMP2) / (Kutipii * Kututp) + PPi * UTP * (PII + PIIUMP + PIIUMP2) / (Kutipii * Kutippi * Kututp) + Kutpii * (PIIUMP + PIIUMP2 + PIIUMP3) / (Kutipii * Kutpiiump))));
vur2 = compartment * (kcatur * UR * PIIUMP2 / (Kurpiiump * (1 + Kglnur / GLN) * (1 + (1 + UMP / Kurump) * (PIIUMP + PIIUMP2 + PIIUMP3) / Kurpiiump)));
vut3 = compartment * (kcatut * UT * UTP * PIIUMP2 / (Kutipii * Kututp * (1 + GLN / Kglnut) * (1 + UTP / Kututp + (PII + PIIUMP + PIIUMP2) / Kutipii + UTP * (PII + PIIUMP + PIIUMP2) / (Kutipii * Kututp) + PPi * UTP * (PII + PIIUMP + PIIUMP2) / (Kutipii * Kutippi * Kututp) + Kutpii * (PIIUMP + PIIUMP2 + PIIUMP3) / (Kutipii * Kutpiiump))));
vur3 = compartment * (kcatur * UR * PIIUMP3 / (Kurpiiump * (1 + Kglnur / GLN) * (1 + (1 + UMP / Kurump) * (PIIUMP + PIIUMP2 + PIIUMP3) / Kurpiiump)));
vad = compartment * (Vad * GS * (b1 * GLN / Kadgln + 3 * a1 * KG * PII / (Kadpiikg * Kd1 * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3))) + 3 * c1 * KG * GLN * PII / (Kadgln * Kadpiikg * Kd1 * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3)))) / ((Kadgs + GS) * (1 + GLN / Kadgln + 3 * KG * PII / (Kadpiikg * Kd1 * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3))) + 3 * KG * GLN * PII / (d1 * Kadgln * Kadpiikg * Kd1 * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3))))));
vdead = compartment * (Vdead * AMP * (f1 * GLN / Kdeadgln + 3 * e1 * KG * PII / (Kd1 * Kdeadpiikg * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3))) + 3 * h1 * KG * GLN * PII / (Kd1 * Kdeadgln * Kdeadpiikg * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3))) + g1 * power(KG, 3) * PIIUMP3 / (Kd1piiump * Kd2piiump * Kd3piiump * Kdeadpiiu * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump))) + j1 * power(KG, 3) * GLN * PIIUMP3 / (Kd1piiump * Kd2piiump * Kd3piiump * Kdeadgln * Kdeadpiiu * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump))) + 3 * i1 * power(KG, 4) * PII * PIIUMP3 / (Kd1 * Kd1piiump * Kd2piiump * Kd3piiump * Kdeadpiikg * Kdeadpiiu * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3)) * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump))) + 3 * k1 * power(KG, 4) * GLN * PII * PIIUMP3 / (Kd1 * Kd1piiump * Kd2piiump * Kd3piiump * Kdeadgln * Kdeadpiikg * Kdeadpiiu * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3)) * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump)))) / ((Kdeadgsa + AMP) * (1 + GLN / Kdeadgln + 3 * KG * PII / (Kd1 * Kdeadpiikg * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3))) + 3 * KG * GLN * PII / (Kd1 * Kdeadgln * Kdeadpiikg * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3)) * l1) + power(KG, 3) * PIIUMP3 / (Kd1piiump * Kd2piiump * Kd3piiump * Kdeadpiiu * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump))) + power(KG, 3) * GLN * PIIUMP3 / (Kd1piiump * Kd2piiump * Kd3piiump * Kdeadgln * Kdeadpiiu * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump)) * n1) + 3 * power(KG, 4) * PII * PIIUMP3 / (Kd1 * Kd1piiump * Kd2piiump * Kd3piiump * Kdeadpiikg * Kdeadpiiu * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3)) * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump)) * m1) + 3 * power(KG, 4) * GLN * PII * PIIUMP3 / (Kd1 * Kd1piiump * Kd2piiump * Kd3piiump * Kdeadgln * Kdeadpiikg * Kdeadpiiu * (1 + 3 * KG / Kd1 + 3 * power(KG, 2) / (Kd1 * Kd2) + power(KG, 3) / (Kd1 * Kd2 * Kd3)) * (1 + 3 * KG / Kd1piiump + 3 * power(KG, 2) / (Kd1piiump * Kd2piiump) + power(KG, 3) / (Kd1piiump * Kd2piiump * Kd3piiump)) * o1))));
vgdh = compartment * (Vgdh * (KG * NADPH * NH4 - NADP * GLU / Keqgdh) / (Kgdhkg * Kgdhnadph * Kgdhnh * (1 + NADP / Kgdhnadp + NADPH / Kgdhnadph) * (1 + NH4 / Kgdhnh) * (1 + KG / Kgdhkg + GLU / Kgdhglu)));
vgog = compartment * (KG * NADPH * Vgog * GLN / (Kgoggln * Kgogkg * Kgognadph * (1 + NADP / Kgognadp + NADPH / Kgognadph) * (1 + AZGLU / Kgogaz) * (1 + KG / Kgogkg + GLU / Kgogglu) * (1 + GLN / Kgoggln + GLU / Kgogglu)));
vgs = compartment * (aamp * camp * Vgs * (-(P_i * ADP * GLN / Keq) + NH4 * ATP * GLU) / (Kgsatp * Kgsglu * Kgsnh * (1 + P_i / Kgspi + ADP / Kgsadp + P_i * ADP / (Kgsadp * Kgspi) + ATP / Kgsatp) * (1 + NH4 / Kgsnh + GLN / Kgsgln + NH4 * GLN / (Kgsgln * Kgsnh) + GLU / Kgsglu + NH4 * GLU / (Kgsglu * Kgsnh)) * (1 + power(12, n1amp) * power(AMP / (bamp * GStot), n1amp)) * (1 + power(12, n2amp) * power(AMP / (damp * GStot), n2amp))));
vgludem = compartment * (Vgludem * (-(AZGLU / Kgludemeq) + GLU) / (Kgludemglu * (1 + AZGLU / Kgludemazglu + GLU / Kgludemglu)));
vazgludem = compartment * (Vazgludem * (-(AZglu / Kazgludemeq) + AZGLU) / (Kazgludemazglu * (1 + AZglu / Kazgludemazinter + AZGLU / Kazgludemazglu)));
vglndem = compartment * (Vglndem * (-(AZGLN / Kglndemeq) + GLN) / (Kglndemgln * (1 + AZGLN / Kglndemazgln + GLN / Kglndemgln)));
vazglndem = compartment * (Vazglndem * (-(AZgln / Kazglndemeq) + AZGLN) / (Kazglndemazgln * (1 + AZgln / Kazglndemazinter + AZGLN / Kazglndemazgln)));
vatpase = compartment * (Vadp * ADP / (Kadp + ADP));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PII_dot = (-vut1+vur1)/compartment;
UTP_dot = 0;
PIIUMP_dot = (+vut1-vur1-vut2+vur2)/compartment;
PPi_dot = 0;
GLN_dot = (-vgog+vgs-vglndem)/compartment;
PIIUMP2_dot = (+vut2-vur2-vut3+vur3)/compartment;
PIIUMP3_dot = (+vut3-vur3)/compartment;
UMP_dot = 0;
GS_dot = (-vad+vdead)/compartment;
AMP_dot = (+vad-vdead)/compartment;
NH4_dot = 0;
KG_dot = 0;
NADPH_dot = 0;
GLU_dot = (+vgdh+2*vgog-vgs-vgludem)/compartment;
NADP_dot = 0;
AZGLU_dot = (+vgludem-vazgludem)/compartment;
ATP_dot = (-vgs+vatpase)/compartment;
ADP_dot = (+vgs-vatpase)/compartment;
AZglu_dot = 0;
AZGLN_dot = (+vglndem-vazglndem)/compartment;
AZgln_dot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = PII_dot;
output(2) = UTP_dot;
output(3) = PIIUMP_dot;
output(4) = PPi_dot;
output(5) = GLN_dot;
output(6) = PIIUMP2_dot;
output(7) = PIIUMP3_dot;
output(8) = UMP_dot;
output(9) = GS_dot;
output(10) = AMP_dot;
output(11) = NH4_dot;
output(12) = KG_dot;
output(13) = NADPH_dot;
output(14) = GLU_dot;
output(15) = NADP_dot;
output(16) = AZGLU_dot;
output(17) = ATP_dot;
output(18) = ADP_dot;
output(19) = AZglu_dot;
output(20) = AZGLN_dot;
output(21) = AZgln_dot;
% return a column vector 
output = output(:);
return


