function [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maeda2019_AmmoniumTransportAssimilation
% Generated: 08-Jul-2020 12:04:09
% 
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun() => output = initial conditions in column vector
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun('states') => output = state names in cell-array
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun('parameters') => output = parameter names in cell-array
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun('parametervalues') => output = parameter values in column vector
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun('variablenames') => output = variable names in cell-array
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = Maeda2019_RefinedActive_Kim_20190413_1_mod_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): ADP
% statevector(2): ATP
% statevector(3): GLN
% statevector(4): GLU
% statevector(5): GS
% statevector(6): GSAMP
% statevector(7): GlnB
% statevector(8): GlnBUMP
% statevector(9): GlnBUMP2
% statevector(10): GlnBUMP3
% statevector(11): GlnK
% statevector(12): GlnKUMP
% statevector(13): GlnKUMP2
% statevector(14): GlnKUMP3
% statevector(15): NADP
% statevector(16): NADPH
% statevector(17): NHxint
% statevector(18): PPi
% statevector(19): Pi
% statevector(20): UMP
% statevector(21): UTP
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0.56, 9.6, 2, 80, 0.0108376, 0, 0.00065, 0, 0, 0, ...
		0.00203338, 0, 0, 0, 0.076, 0.12, 0.02, 0.05, 10, 0.01, ...
		8.3];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'ADP', 'ATP', 'GLN', 'GLU', 'GS', 'GSAMP', 'GlnB', 'GlnBUMP', 'GlnBUMP2', 'GlnBUMP3', ...
			'GlnK', 'GlnKUMP', 'GlnKUMP2', 'GlnKUMP3', 'NADP', 'NADPH', 'NHxint', 'PPi', 'Pi', 'UMP', ...
			'UTP'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'Acell', 'AmtB', 'Dpsi', 'EXTERNAL', 'F', 'GLNdemf', 'GLNdemn', 'GLUdemf', 'GLUdemn', 'Kadgln', ...
			'Kadglnbog', 'Kadgs', 'Kamtbnh', 'Kdeadgln', 'Kdeadglnbog', 'Kdeadglnbump', 'Kdeadgsamp', 'Kgdheq', 'Kgdhglu', 'Kgdhnadp', ...
			'Kgdhnadph', 'Kgdhnh', 'Kgdhog', 'Kglnbog1', 'Kglnbog2', 'Kglnbog3', 'Kglnbumpog1', 'Kglnbumpog2', 'Kglnbumpog3', 'Kglnkamtb', ...
			'Kglnkog1', 'Kglnkog2', 'Kglnkog3', 'Kgoggln', 'Kgogglu', 'Kgognadp', 'Kgognadph', 'Kgogog', 'Kgrowthgln', 'Kgrowthglu', ...
			'Kgsadp', 'Kgsatp', 'Kgseq', 'Kgsgln', 'Kgsglu', 'Kgsnh', 'Kgspi', 'Kurgln', 'Kurglnbump', 'Kurglnkump', ...
			'Kurump', 'Kutgln', 'Kutglnb', 'Kutglnk', 'Kutiglnb', 'Kutiglnbump', 'Kutiglnk', 'Kutiglnkump', 'Kutippi', 'Kututp', ...
			'NHxext', 'Nintstar', 'OGbasal', 'Pcm', 'R', 'T', 'UTase', 'Vad', 'Vcell', 'Vdead', ...
			'Vgdh', 'Vgog', 'a1', 'aamp', 'b1', 'bamp', 'c1', 'camp', 'd1', 'damp', ...
			'e1', 'f1', 'g1', 'h1', 'i1', 'j1', 'k1', 'kappa', 'kcatamtb', 'kcatgs', ...
			'kcaturglnb', 'kcaturglnk', 'kcatutglnb', 'kcatutglnk', 'kdb', 'l1', 'm1', 'n1', 'n1amp', 'n2amp', ...
			'o1', 'pHext', 'pHint', 'pKa', 'tau0', 'default0'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [9.18e-12, 0.00168395, -0.15, 0, 96485, 84.6574, 645.163, 316.076, 2116.81, 0.952511, ...
			1.19929e-05, 7.48531e-05, 0.005, 0.202364, 2.58892e-06, 7.07454e-06, 0.00044786, 1290, 6.27836, 0.0348261, ...
			0.0485857, 1.1, 0.518879, 0.0049519, 0.137803, 0.148798, 0.0241288, 0.151613, 0.132972, 6.26171e-08, ...
			9.5247, 5.41511, 5.19229, 0.287442, 6.92079, 0.00339109, 0.00164105, 0.006879, 1.50715, 31.1697, ...
			0.0730688, 0.264057, 460, 5.81152, 4.12683, 0.1, 4.51385, 0.063601, 0.0046321, 0.00197953, ...
			8.33975, 0.0640159, 0.00292567, 0.00456417, 0.00185694, 0.00363233, 0.00176846, 0.0024576, 0.107936, 0.0417302, ...
			0.00411274, 0.033, 0.550038, 0.0733905, 8.314, 310, 0.0006, 0.418331, 2.15e-18, 0.659368, ...
			352.64, 83.1933, 1e-22, 10, 0.753263, 2.3667, 0.295923, 0.1012, 0.0158453, 10.8688, ...
			1e-22, 0.662946, 14.4267, 0.20749, 1e-22, 1e-22, 1e-22, 7.88981, 795355, 53049.4, ...
			2.81579, 17.3128, 132.248, 74.2685, 13.9379, 0.017405, 0.87943, 9.96306, 1.1456, 19.2166, ...
			1.29171, 7.4, 7.6, 8.95, 45.8312, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'GStotal', 'Vgs', 'Ka', 'Hint', 'Hext', 'NH4int', 'OG', 'GlnBOG1', 'GlnB_OGfree', 'GlnKOG2', ...
			'GlnBUMP3OG3', 'GlnK_OGfree', 'AmtB_GlnKfree', 'GlnKOG3', 'GlnKOG1', 'GlnKAmtB', 'GlnK_AmtBfree', 'NHxsurf', 'NH4surf', 'NH3int', ...
			'NH4ext', 'Vamtb_app', 'Vamtb', 'NH3ext', 'NH3surf', 'theta_ad', 'Vad_app', 'theta_dead', 'Vdead_app', 'nAMP', ...
			'theta_gs', 'Vgs_app', 'phi', 'kdiff', 'tau', 'mu'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'GS+GSAMP', 'kcatgs*GStotal', 'power(10,-pKa)*1000', 'power(10,-pHint)*1000', 'power(10,-pHext)*1000', 'NHxint*Hint/(Ka+Hint)', 'kappa*(1-NH4int/Nintstar)+OGbasal', '3*GlnB*OG/Kglnbog1/(1+3*OG/Kglnbog1+3*power(OG,2)/(Kglnbog1*Kglnbog2)+power(OG,3)/(Kglnbog1*Kglnbog2*Kglnbog3))', 'GlnB/(1+3*OG/Kglnbog1+3*power(OG,2)/(Kglnbog1*Kglnbog2)+power(OG,3)/(Kglnbog1*Kglnbog2*Kglnbog3))', '3*GlnK*power(OG,2)/(Kglnkog1*Kglnkog2)/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3))', ...
			'GlnBUMP3*power(OG,3)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3)/(1+3*OG/Kglnbumpog1+3*power(OG,2)/(Kglnbumpog1*Kglnbumpog2)+power(OG,3)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3))', 'GlnK/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3))', '0.5*(-GlnK_OGfree+AmtB-Kglnkamtb+power((GlnK_OGfree-AmtB+Kglnkamtb)*(GlnK_OGfree-AmtB+Kglnkamtb)+4*Kglnkamtb*AmtB,0.5))', 'GlnK*power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3)/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3))', '3*GlnK*OG/Kglnkog1/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3))', 'AmtB-AmtB_GlnKfree', 'GlnK-GlnKAmtB', 'NHxext', 'NHxsurf*Hext/(Ka+Hext)', 'NHxint*Ka/(Ka+Hint)', ...
			'NHxext*Hext/(Ka+Hext)', 'kcatamtb*AmtB_GlnKfree', 'kcatamtb*AmtB', 'NHxext*Ka/(Ka+Hext)', 'NHxsurf*Ka/(Ka+Hext)', '(a1*GlnBOG1/Kadglnbog+b1*GLN/Kadgln+c1*GlnBOG1*GLN/(Kadglnbog*Kadgln))/(1+GlnBOG1/Kadglnbog+GLN/Kadgln+GlnBOG1*GLN/(d1*Kadglnbog*Kadgln))', 'theta_ad*Vad', '(e1*GlnBOG1/Kdeadglnbog+f1*GLN/Kdeadgln+g1*GlnBUMP3OG3/Kdeadglnbump+h1*GlnBOG1*GLN/(Kdeadglnbog*Kdeadgln)+i1*GlnBOG1*GlnBUMP3OG3/(Kdeadglnbog*Kdeadglnbump)+j1*GLN*GlnBUMP3OG3/(Kdeadgln*Kdeadglnbump)+k1*GlnBOG1*GLN*GlnBUMP3OG3/(Kdeadglnbog*Kdeadgln*Kdeadglnbump))/(1+GlnBOG1/Kdeadglnbog+GLN/Kdeadgln+GlnBUMP3OG3/Kdeadglnbump+GlnBOG1*GLN/(l1*Kdeadglnbog*Kdeadgln)+GlnBOG1*GlnBUMP3OG3/(m1*Kdeadglnbog*Kdeadglnbump)+GLN*GlnBUMP3OG3/(n1*Kdeadgln*Kdeadglnbump)+GlnBOG1*GLN*GlnBUMP3OG3/(o1*Kdeadglnbog*Kdeadgln*Kdeadglnbump))', 'theta_dead*Vdead', '12*GSAMP/GStotal', ...
			'aamp/(1+power(nAMP/bamp,n1amp))*camp/(1+power(nAMP/damp,n2amp))', 'theta_gs*Vgs', 'exp(-F*Dpsi/(R*T))', 'Pcm*Acell/Vcell', 'tau0*(1+power(Kgrowthglu/GLU,2)+power(Kgrowthgln/GLN,2))', 'log(2)/tau'};
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
	if length(parameterValuesNew) ~= 106,
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
ADP = statevector(1);
ATP = statevector(2);
GLN = statevector(3);
GLU = statevector(4);
GS = statevector(5);
GSAMP = statevector(6);
GlnB = statevector(7);
GlnBUMP = statevector(8);
GlnBUMP2 = statevector(9);
GlnBUMP3 = statevector(10);
GlnK = statevector(11);
GlnKUMP = statevector(12);
GlnKUMP2 = statevector(13);
GlnKUMP3 = statevector(14);
NADP = statevector(15);
NADPH = statevector(16);
NHxint = statevector(17);
PPi = statevector(18);
Pi = statevector(19);
UMP = statevector(20);
UTP = statevector(21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	Acell = 9.18e-12;
	AmtB = 0.00168395;
	Dpsi = -0.15;
	EXTERNAL = 0;
	F = 96485;
	GLNdemf = 84.6574;
	GLNdemn = 645.163;
	GLUdemf = 316.076;
	GLUdemn = 2116.81;
	Kadgln = 0.952511;
	Kadglnbog = 1.19929e-05;
	Kadgs = 7.48531e-05;
	Kamtbnh = 0.005;
	Kdeadgln = 0.202364;
	Kdeadglnbog = 2.58892e-06;
	Kdeadglnbump = 7.07454e-06;
	Kdeadgsamp = 0.00044786;
	Kgdheq = 1290;
	Kgdhglu = 6.27836;
	Kgdhnadp = 0.0348261;
	Kgdhnadph = 0.0485857;
	Kgdhnh = 1.1;
	Kgdhog = 0.518879;
	Kglnbog1 = 0.0049519;
	Kglnbog2 = 0.137803;
	Kglnbog3 = 0.148798;
	Kglnbumpog1 = 0.0241288;
	Kglnbumpog2 = 0.151613;
	Kglnbumpog3 = 0.132972;
	Kglnkamtb = 6.26171e-08;
	Kglnkog1 = 9.5247;
	Kglnkog2 = 5.41511;
	Kglnkog3 = 5.19229;
	Kgoggln = 0.287442;
	Kgogglu = 6.92079;
	Kgognadp = 0.00339109;
	Kgognadph = 0.00164105;
	Kgogog = 0.006879;
	Kgrowthgln = 1.50715;
	Kgrowthglu = 31.1697;
	Kgsadp = 0.0730688;
	Kgsatp = 0.264057;
	Kgseq = 460;
	Kgsgln = 5.81152;
	Kgsglu = 4.12683;
	Kgsnh = 0.1;
	Kgspi = 4.51385;
	Kurgln = 0.063601;
	Kurglnbump = 0.0046321;
	Kurglnkump = 0.00197953;
	Kurump = 8.33975;
	Kutgln = 0.0640159;
	Kutglnb = 0.00292567;
	Kutglnk = 0.00456417;
	Kutiglnb = 0.00185694;
	Kutiglnbump = 0.00363233;
	Kutiglnk = 0.00176846;
	Kutiglnkump = 0.0024576;
	Kutippi = 0.107936;
	Kututp = 0.0417302;
	NHxext = 0.00411274;
	Nintstar = 0.033;
	OGbasal = 0.550038;
	Pcm = 0.0733905;
	R = 8.314;
	T = 310;
	UTase = 0.0006;
	Vad = 0.418331;
	Vcell = 2.15e-18;
	Vdead = 0.659368;
	Vgdh = 352.64;
	Vgog = 83.1933;
	a1 = 1e-22;
	aamp = 10;
	b1 = 0.753263;
	bamp = 2.3667;
	c1 = 0.295923;
	camp = 0.1012;
	d1 = 0.0158453;
	damp = 10.8688;
	e1 = 1e-22;
	f1 = 0.662946;
	g1 = 14.4267;
	h1 = 0.20749;
	i1 = 1e-22;
	j1 = 1e-22;
	k1 = 1e-22;
	kappa = 7.88981;
	kcatamtb = 795355;
	kcatgs = 53049.4;
	kcaturglnb = 2.81579;
	kcaturglnk = 17.3128;
	kcatutglnb = 132.248;
	kcatutglnk = 74.2685;
	kdb = 13.9379;
	l1 = 0.017405;
	m1 = 0.87943;
	n1 = 9.96306;
	n1amp = 1.1456;
	n2amp = 19.2166;
	o1 = 1.29171;
	pHext = 7.4;
	pHint = 7.6;
	pKa = 8.95;
	tau0 = 45.8312;
	default0 = 1;
else
	Acell = parameterValuesNew(1);
	AmtB = parameterValuesNew(2);
	Dpsi = parameterValuesNew(3);
	EXTERNAL = parameterValuesNew(4);
	F = parameterValuesNew(5);
	GLNdemf = parameterValuesNew(6);
	GLNdemn = parameterValuesNew(7);
	GLUdemf = parameterValuesNew(8);
	GLUdemn = parameterValuesNew(9);
	Kadgln = parameterValuesNew(10);
	Kadglnbog = parameterValuesNew(11);
	Kadgs = parameterValuesNew(12);
	Kamtbnh = parameterValuesNew(13);
	Kdeadgln = parameterValuesNew(14);
	Kdeadglnbog = parameterValuesNew(15);
	Kdeadglnbump = parameterValuesNew(16);
	Kdeadgsamp = parameterValuesNew(17);
	Kgdheq = parameterValuesNew(18);
	Kgdhglu = parameterValuesNew(19);
	Kgdhnadp = parameterValuesNew(20);
	Kgdhnadph = parameterValuesNew(21);
	Kgdhnh = parameterValuesNew(22);
	Kgdhog = parameterValuesNew(23);
	Kglnbog1 = parameterValuesNew(24);
	Kglnbog2 = parameterValuesNew(25);
	Kglnbog3 = parameterValuesNew(26);
	Kglnbumpog1 = parameterValuesNew(27);
	Kglnbumpog2 = parameterValuesNew(28);
	Kglnbumpog3 = parameterValuesNew(29);
	Kglnkamtb = parameterValuesNew(30);
	Kglnkog1 = parameterValuesNew(31);
	Kglnkog2 = parameterValuesNew(32);
	Kglnkog3 = parameterValuesNew(33);
	Kgoggln = parameterValuesNew(34);
	Kgogglu = parameterValuesNew(35);
	Kgognadp = parameterValuesNew(36);
	Kgognadph = parameterValuesNew(37);
	Kgogog = parameterValuesNew(38);
	Kgrowthgln = parameterValuesNew(39);
	Kgrowthglu = parameterValuesNew(40);
	Kgsadp = parameterValuesNew(41);
	Kgsatp = parameterValuesNew(42);
	Kgseq = parameterValuesNew(43);
	Kgsgln = parameterValuesNew(44);
	Kgsglu = parameterValuesNew(45);
	Kgsnh = parameterValuesNew(46);
	Kgspi = parameterValuesNew(47);
	Kurgln = parameterValuesNew(48);
	Kurglnbump = parameterValuesNew(49);
	Kurglnkump = parameterValuesNew(50);
	Kurump = parameterValuesNew(51);
	Kutgln = parameterValuesNew(52);
	Kutglnb = parameterValuesNew(53);
	Kutglnk = parameterValuesNew(54);
	Kutiglnb = parameterValuesNew(55);
	Kutiglnbump = parameterValuesNew(56);
	Kutiglnk = parameterValuesNew(57);
	Kutiglnkump = parameterValuesNew(58);
	Kutippi = parameterValuesNew(59);
	Kututp = parameterValuesNew(60);
	NHxext = parameterValuesNew(61);
	Nintstar = parameterValuesNew(62);
	OGbasal = parameterValuesNew(63);
	Pcm = parameterValuesNew(64);
	R = parameterValuesNew(65);
	T = parameterValuesNew(66);
	UTase = parameterValuesNew(67);
	Vad = parameterValuesNew(68);
	Vcell = parameterValuesNew(69);
	Vdead = parameterValuesNew(70);
	Vgdh = parameterValuesNew(71);
	Vgog = parameterValuesNew(72);
	a1 = parameterValuesNew(73);
	aamp = parameterValuesNew(74);
	b1 = parameterValuesNew(75);
	bamp = parameterValuesNew(76);
	c1 = parameterValuesNew(77);
	camp = parameterValuesNew(78);
	d1 = parameterValuesNew(79);
	damp = parameterValuesNew(80);
	e1 = parameterValuesNew(81);
	f1 = parameterValuesNew(82);
	g1 = parameterValuesNew(83);
	h1 = parameterValuesNew(84);
	i1 = parameterValuesNew(85);
	j1 = parameterValuesNew(86);
	k1 = parameterValuesNew(87);
	kappa = parameterValuesNew(88);
	kcatamtb = parameterValuesNew(89);
	kcatgs = parameterValuesNew(90);
	kcaturglnb = parameterValuesNew(91);
	kcaturglnk = parameterValuesNew(92);
	kcatutglnb = parameterValuesNew(93);
	kcatutglnk = parameterValuesNew(94);
	kdb = parameterValuesNew(95);
	l1 = parameterValuesNew(96);
	m1 = parameterValuesNew(97);
	n1 = parameterValuesNew(98);
	n1amp = parameterValuesNew(99);
	n2amp = parameterValuesNew(100);
	o1 = parameterValuesNew(101);
	pHext = parameterValuesNew(102);
	pHint = parameterValuesNew(103);
	pKa = parameterValuesNew(104);
	tau0 = parameterValuesNew(105);
	default0 = parameterValuesNew(106);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GStotal = GS+GSAMP;
Vgs = kcatgs*GStotal;
Ka = power(10,-pKa)*1000;
Hint = power(10,-pHint)*1000;
Hext = power(10,-pHext)*1000;
NH4int = NHxint*Hint/(Ka+Hint);
OG = kappa*(1-NH4int/Nintstar)+OGbasal;
GlnBOG1 = 3*GlnB*OG/Kglnbog1/(1+3*OG/Kglnbog1+3*power(OG,2)/(Kglnbog1*Kglnbog2)+power(OG,3)/(Kglnbog1*Kglnbog2*Kglnbog3));
GlnB_OGfree = GlnB/(1+3*OG/Kglnbog1+3*power(OG,2)/(Kglnbog1*Kglnbog2)+power(OG,3)/(Kglnbog1*Kglnbog2*Kglnbog3));
GlnKOG2 = 3*GlnK*power(OG,2)/(Kglnkog1*Kglnkog2)/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3));
GlnBUMP3OG3 = GlnBUMP3*power(OG,3)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3)/(1+3*OG/Kglnbumpog1+3*power(OG,2)/(Kglnbumpog1*Kglnbumpog2)+power(OG,3)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3));
GlnK_OGfree = GlnK/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3));
AmtB_GlnKfree = 0.5*(-GlnK_OGfree+AmtB-Kglnkamtb+power((GlnK_OGfree-AmtB+Kglnkamtb)*(GlnK_OGfree-AmtB+Kglnkamtb)+4*Kglnkamtb*AmtB,0.5));
GlnKOG3 = GlnK*power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3)/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3));
GlnKOG1 = 3*GlnK*OG/Kglnkog1/(1+3*OG/Kglnkog1+3*power(OG,2)/(Kglnkog1*Kglnkog2)+power(OG,3)/(Kglnkog1*Kglnkog2*Kglnkog3));
GlnKAmtB = AmtB-AmtB_GlnKfree;
GlnK_AmtBfree = GlnK-GlnKAmtB;
NHxsurf = NHxext;
NH4surf = NHxsurf*Hext/(Ka+Hext);
NH3int = NHxint*Ka/(Ka+Hint);
NH4ext = NHxext*Hext/(Ka+Hext);
Vamtb_app = kcatamtb*AmtB_GlnKfree;
Vamtb = kcatamtb*AmtB;
NH3ext = NHxext*Ka/(Ka+Hext);
NH3surf = NHxsurf*Ka/(Ka+Hext);
theta_ad = (a1*GlnBOG1/Kadglnbog+b1*GLN/Kadgln+c1*GlnBOG1*GLN/(Kadglnbog*Kadgln))/(1+GlnBOG1/Kadglnbog+GLN/Kadgln+GlnBOG1*GLN/(d1*Kadglnbog*Kadgln));
Vad_app = theta_ad*Vad;
theta_dead = (e1*GlnBOG1/Kdeadglnbog+f1*GLN/Kdeadgln+g1*GlnBUMP3OG3/Kdeadglnbump+h1*GlnBOG1*GLN/(Kdeadglnbog*Kdeadgln)+i1*GlnBOG1*GlnBUMP3OG3/(Kdeadglnbog*Kdeadglnbump)+j1*GLN*GlnBUMP3OG3/(Kdeadgln*Kdeadglnbump)+k1*GlnBOG1*GLN*GlnBUMP3OG3/(Kdeadglnbog*Kdeadgln*Kdeadglnbump))/(1+GlnBOG1/Kdeadglnbog+GLN/Kdeadgln+GlnBUMP3OG3/Kdeadglnbump+GlnBOG1*GLN/(l1*Kdeadglnbog*Kdeadgln)+GlnBOG1*GlnBUMP3OG3/(m1*Kdeadglnbog*Kdeadglnbump)+GLN*GlnBUMP3OG3/(n1*Kdeadgln*Kdeadglnbump)+GlnBOG1*GLN*GlnBUMP3OG3/(o1*Kdeadglnbog*Kdeadgln*Kdeadglnbump));
Vdead_app = theta_dead*Vdead;
nAMP = 12*GSAMP/GStotal;
theta_gs = aamp/(1+power(nAMP/bamp,n1amp))*camp/(1+power(nAMP/damp,n2amp));
Vgs_app = theta_gs*Vgs;
phi = exp(-F*Dpsi/(R*T));
kdiff = Pcm*Acell/Vcell;
tau = tau0*(1+power(Kgrowthglu/GLU,2)+power(Kgrowthgln/GLN,2));
mu = log(2)/tau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vad = default0 * Function_for_vad(GS, Kadgs, Vad_app, default0);
vamtb = default0 * Function_for_vamtb(Kamtbnh, NH4int, NH4surf, Vamtb_app, default0, phi);
vdead = default0 * Function_for_vdead(GSAMP, Kdeadgsamp, Vdead_app, default0);
vdiff = default0 * Function_for_vdiff(NH3int, NH3surf, default0, kdiff);
vgdh = default0 * Function_for_vgdh(GLU, Kgdheq, Kgdhglu, Kgdhnadp, Kgdhnadph, Kgdhnh, Kgdhog, NADP, NADPH, NH4int, OG, Vgdh, default0);
vglndemf = default0 * Function_for_vglndemf(GLNdemf, default0, mu);
vglndemn = default0 * Function_for_vglndemn(GLNdemn, default0, mu);
vgludemf = default0 * Function_for_vgludemf(GLUdemf, default0, mu);
vgludemn = default0 * Function_for_vgludemn(GLUdemn, default0, mu);
vgog = default0 * Function_for_vgog(GLN, GLU, Kgoggln, Kgogglu, Kgognadp, Kgognadph, Kgogog, NADP, NADPH, OG, Vgog, default0);
vgs = default0 * Function_for_vgs(ADP, ATP, GLN, GLU, Kgsadp, Kgsatp, Kgseq, Kgsgln, Kgsglu, Kgsnh, Kgspi, NH4int, Pi, Vgs_app, default0);
vurglnb1 = default0 * Function_for_vurglnb1(GLN, GlnBUMP, GlnBUMP2, GlnBUMP3, Kurgln, Kurglnbump, Kurump, UMP, UTase, default0, kcaturglnb);
vurglnb2 = default0 * Function_for_vurglnb2(GLN, GlnBUMP, GlnBUMP2, GlnBUMP3, Kurgln, Kurglnbump, Kurump, UMP, UTase, default0, kcaturglnb);
vurglnb3 = default0 * Function_for_vurglnb3(GLN, GlnBUMP, GlnBUMP2, GlnBUMP3, Kurgln, Kurglnbump, Kurump, UMP, UTase, default0, kcaturglnb);
vurglnk1 = default0 * Function_for_vurglnk1(GLN, GlnKUMP, GlnKUMP2, GlnKUMP3, Kurgln, Kurglnkump, Kurump, UMP, UTase, default0, kcaturglnk);
vurglnk2 = default0 * Function_for_vurglnk2(GLN, GlnKUMP, GlnKUMP2, GlnKUMP3, Kurgln, Kurglnkump, Kurump, UMP, UTase, default0, kcaturglnk);
vurglnk3 = default0 * Function_for_vurglnk3(GLN, GlnKUMP, GlnKUMP2, GlnKUMP3, Kurgln, Kurglnkump, Kurump, UMP, UTase, default0, kcaturglnk);
vutglnb1 = default0 * Function_for_vutglnb1(GLN, GlnB, GlnBUMP, GlnBUMP2, GlnBUMP3, Kutgln, Kutglnb, Kutiglnb, Kutiglnbump, Kutippi, Kututp, PPi, UTP, UTase, default0, kcatutglnb);
vutglnb2 = default0 * Function_for_vutglnb2(GLN, GlnB, GlnBUMP, GlnBUMP2, GlnBUMP3, Kutgln, Kutglnb, Kutiglnb, Kutiglnbump, Kutippi, Kututp, PPi, UTP, UTase, default0, kcatutglnb);
vutglnb3 = default0 * Function_for_vutglnb3(GLN, GlnB, GlnBUMP, GlnBUMP2, GlnBUMP3, Kutgln, Kutglnb, Kutiglnb, Kutiglnbump, Kutippi, Kututp, PPi, UTP, UTase, default0, kcatutglnb);
vutglnk1 = default0 * Function_for_vutglnk1(GLN, GlnKUMP, GlnKUMP2, GlnKUMP3, GlnK_AmtBfree, Kutgln, Kutglnk, Kutiglnk, Kutiglnkump, Kutippi, Kututp, PPi, UTP, UTase, default0, kcatutglnk);
vutglnk2 = default0 * Function_for_vutglnk2(GLN, GlnKUMP, GlnKUMP2, GlnKUMP3, GlnK_AmtBfree, Kutgln, Kutglnk, Kutiglnk, Kutiglnkump, Kutippi, Kututp, PPi, UTP, UTase, default0, kcatutglnk);
vutglnk3 = default0 * Function_for_vutglnk3(GLN, GlnKUMP, GlnKUMP2, GlnKUMP3, GlnK_AmtBfree, Kutgln, Kutglnk, Kutiglnk, Kutiglnkump, Kutippi, Kututp, PPi, UTP, UTase, default0, kcatutglnk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADP_dot = 0;
ATP_dot = 0;
GLN_dot = (-vglndemf-vglndemn-vgog+vgs)/default0;
GLU_dot = (+vgdh+vglndemn-vgludemf-vgludemn+2*vgog-vgs)/default0;
GS_dot = (-vad+vdead)/default0;
GSAMP_dot = (+vad-vdead)/default0;
GlnB_dot = (+vurglnb1-vutglnb1)/default0;
GlnBUMP_dot = (-vurglnb1+vurglnb2+vutglnb1-vutglnb2)/default0;
GlnBUMP2_dot = (-vurglnb2+vurglnb3+vutglnb2-vutglnb3)/default0;
GlnBUMP3_dot = (-vurglnb3+vutglnb3)/default0;
GlnK_dot = (+vurglnk1-vutglnk1)/default0;
GlnKUMP_dot = (-vurglnk1+vurglnk2+vutglnk1-vutglnk2)/default0;
GlnKUMP2_dot = (-vurglnk2+vurglnk3+vutglnk2-vutglnk3)/default0;
GlnKUMP3_dot = (-vurglnk3+vutglnk3)/default0;
NADP_dot = 0;
NADPH_dot = 0;
NHxint_dot = (+vamtb+vdiff-vgdh-vgs)/default0;
PPi_dot = 0;
Pi_dot = 0;
UMP_dot = 0;
UTP_dot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = ADP_dot;
output(2) = ATP_dot;
output(3) = GLN_dot;
output(4) = GLU_dot;
output(5) = GS_dot;
output(6) = GSAMP_dot;
output(7) = GlnB_dot;
output(8) = GlnBUMP_dot;
output(9) = GlnBUMP2_dot;
output(10) = GlnBUMP3_dot;
output(11) = GlnK_dot;
output(12) = GlnKUMP_dot;
output(13) = GlnKUMP2_dot;
output(14) = GlnKUMP3_dot;
output(15) = NADP_dot;
output(16) = NADPH_dot;
output(17) = NHxint_dot;
output(18) = PPi_dot;
output(19) = Pi_dot;
output(20) = UMP_dot;
output(21) = UTP_dot;
% return a column vector 
output = output(:);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = Function_for_vamtb(Kamtbnh,NH4int,NH4surf,Vamtb_app,default0,phi)
global time
result = Vamtb_app*(NH4surf-NH4int/phi)/(Kamtbnh+NH4surf)/default0;
return

function [result] = Function_for_vdiff(NH3int,NH3surf,default0,kdiff)
global time
result = kdiff*(NH3surf-NH3int)/default0;
return

function [result] = Function_for_vdead(GSAMP,Kdeadgsamp,Vdead_app,default0)
global time
result = Vdead_app*GSAMP/(Kdeadgsamp+GSAMP)/default0;
return

function [result] = Function_for_vglndemf(GLNdemf,default0,mu)
global time
result = mu*GLNdemf/default0;
return

function [result] = Function_for_vglndemn(GLNdemn,default0,mu)
global time
result = mu*GLNdemn/default0;
return

function [result] = Function_for_vgludemf(GLUdemf,default0,mu)
global time
result = mu*GLUdemf/default0;
return

function [result] = Function_for_vad(GS,Kadgs,Vad_app,default0)
global time
result = Vad_app*GS/(Kadgs+GS)/default0;
return

function [result] = Function_for_vgdh(GLU,Kgdheq,Kgdhglu,Kgdhnadp,Kgdhnadph,Kgdhnh,Kgdhog,NADP,NADPH,NH4int,OG,Vgdh,default0)
global time
result = Vgdh/(Kgdhog*Kgdhnh*Kgdhnadph)*(OG*NH4int*NADPH-GLU*NADP/Kgdheq)/((1+NH4int/Kgdhnh)*(1+OG/Kgdhog+GLU/Kgdhglu)*(1+NADPH/Kgdhnadph+NADP/Kgdhnadp))/default0;
return

function [result] = Function_for_vgludemn(GLUdemn,default0,mu)
global time
result = mu*GLUdemn/default0;
return

function [result] = Function_for_vutglnk2(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,GlnK_AmtBfree,Kutgln,Kutglnk,Kutiglnk,Kutiglnkump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnk)
global time
result = kcatutglnk*UTase*GlnKUMP*UTP/((1+GLN/Kutgln)*(Kutiglnk*Kututp+Kututp*(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)+Kutglnk*UTP+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP+Kutglnk*UTP*(GlnKUMP+GlnKUMP2+GlnKUMP3)/Kutiglnkump+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP*PPi/Kutippi))/default0;
return

function [result] = Function_for_vutglnk3(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,GlnK_AmtBfree,Kutgln,Kutglnk,Kutiglnk,Kutiglnkump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnk)
global time
result = kcatutglnk*UTase*GlnKUMP2*UTP/((1+GLN/Kutgln)*(Kutiglnk*Kututp+Kututp*(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)+Kutglnk*UTP+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP+Kutglnk*UTP*(GlnKUMP+GlnKUMP2+GlnKUMP3)/Kutiglnkump+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP*PPi/Kutippi))/default0;
return

function [result] = Function_for_vurglnb1(GLN,GlnBUMP,GlnBUMP2,GlnBUMP3,Kurgln,Kurglnbump,Kurump,UMP,UTase,default0,kcaturglnb)
global time
result = kcaturglnb*UTase*GlnBUMP/((1+Kurgln/GLN)*(Kurglnbump+(1+UMP/Kurump)*(GlnBUMP+GlnBUMP2+GlnBUMP3)))/default0;
return

function [result] = Function_for_vurglnb2(GLN,GlnBUMP,GlnBUMP2,GlnBUMP3,Kurgln,Kurglnbump,Kurump,UMP,UTase,default0,kcaturglnb)
global time
result = kcaturglnb*UTase*GlnBUMP2/((1+Kurgln/GLN)*(Kurglnbump+(1+UMP/Kurump)*(GlnBUMP+GlnBUMP2+GlnBUMP3)))/default0;
return

function [result] = Function_for_vutglnb3(GLN,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,Kutgln,Kutglnb,Kutiglnb,Kutiglnbump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnb)
global time
result = kcatutglnb*UTase*GlnBUMP2*UTP/((1+GLN/Kutgln)*(Kutiglnb*Kututp+Kututp*(GlnB+GlnBUMP+GlnBUMP2)+Kutglnb*UTP+(GlnB+GlnBUMP+GlnBUMP2)*UTP+Kutglnb*UTP*(GlnBUMP+GlnBUMP2+GlnBUMP3)/Kutiglnbump+(GlnB+GlnBUMP+GlnBUMP2)*UTP*PPi/Kutippi))/default0;
return

function [result] = Function_for_vurglnk3(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,Kurgln,Kurglnkump,Kurump,UMP,UTase,default0,kcaturglnk)
global time
result = kcaturglnk*UTase*GlnKUMP3/((1+Kurgln/GLN)*(Kurglnkump+(1+UMP/Kurump)*(GlnKUMP+GlnKUMP2+GlnKUMP3)))/default0;
return

function [result] = Function_for_vutglnk1(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,GlnK_AmtBfree,Kutgln,Kutglnk,Kutiglnk,Kutiglnkump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnk)
global time
result = kcatutglnk*UTase*GlnK_AmtBfree*UTP/((1+GLN/Kutgln)*(Kutiglnk*Kututp+Kututp*(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)+Kutglnk*UTP+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP+Kutglnk*UTP*(GlnKUMP+GlnKUMP2+GlnKUMP3)/Kutiglnkump+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP*PPi/Kutippi))/default0;
return

function [result] = Function_for_vutglnb1(GLN,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,Kutgln,Kutglnb,Kutiglnb,Kutiglnbump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnb)
global time
result = kcatutglnb*UTase*GlnB*UTP/((1+GLN/Kutgln)*(Kutiglnb*Kututp+Kututp*(GlnB+GlnBUMP+GlnBUMP2)+Kutglnb*UTP+(GlnB+GlnBUMP+GlnBUMP2)*UTP+Kutglnb*UTP*(GlnBUMP+GlnBUMP2+GlnBUMP3)/Kutiglnbump+(GlnB+GlnBUMP+GlnBUMP2)*UTP*PPi/Kutippi))/default0;
return

function [result] = Function_for_vutglnb2(GLN,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,Kutgln,Kutglnb,Kutiglnb,Kutiglnbump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnb)
global time
result = kcatutglnb*UTase*GlnBUMP*UTP/((1+GLN/Kutgln)*(Kutiglnb*Kututp+Kututp*(GlnB+GlnBUMP+GlnBUMP2)+Kutglnb*UTP+(GlnB+GlnBUMP+GlnBUMP2)*UTP+Kutglnb*UTP*(GlnBUMP+GlnBUMP2+GlnBUMP3)/Kutiglnbump+(GlnB+GlnBUMP+GlnBUMP2)*UTP*PPi/Kutippi))/default0;
return

function [result] = Function_for_vurglnb3(GLN,GlnBUMP,GlnBUMP2,GlnBUMP3,Kurgln,Kurglnbump,Kurump,UMP,UTase,default0,kcaturglnb)
global time
result = kcaturglnb*UTase*GlnBUMP3/((1+Kurgln/GLN)*(Kurglnbump+(1+UMP/Kurump)*(GlnBUMP+GlnBUMP2+GlnBUMP3)))/default0;
return

function [result] = Function_for_vgog(GLN,GLU,Kgoggln,Kgogglu,Kgognadp,Kgognadph,Kgogog,NADP,NADPH,OG,Vgog,default0)
global time
result = Vgog*GLN*OG*NADPH/(Kgoggln*Kgogog*Kgognadph)/((1+GLN/Kgoggln+GLU/Kgogglu)*(1+OG/Kgogog+GLU/Kgogglu)*(1+NADPH/Kgognadph+NADP/Kgognadp))/default0;
return

function [result] = Function_for_vurglnk1(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,Kurgln,Kurglnkump,Kurump,UMP,UTase,default0,kcaturglnk)
global time
result = kcaturglnk*UTase*GlnKUMP/((1+Kurgln/GLN)*(Kurglnkump+(1+UMP/Kurump)*(GlnKUMP+GlnKUMP2+GlnKUMP3)))/default0;
return

function [result] = Function_for_vgs(ADP,ATP,GLN,GLU,Kgsadp,Kgsatp,Kgseq,Kgsgln,Kgsglu,Kgsnh,Kgspi,NH4int,Pi,Vgs_app,default0)
global time
result = Vgs_app/(Kgsatp*Kgsnh*Kgsglu)*(ATP*NH4int*GLU-ADP*GLN*Pi/Kgseq)/((1+ATP/Kgsatp+ADP/Kgsadp+Pi/Kgspi+ADP*Pi/(Kgsadp*Kgspi))*(1+NH4int/Kgsnh+GLN/Kgsgln+GLU/Kgsglu+GLN*NH4int/(Kgsgln*Kgsnh)+GLU*NH4int/(Kgsglu*Kgsnh)))/default0;
return

function [result] = Function_for_vurglnk2(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,Kurgln,Kurglnkump,Kurump,UMP,UTase,default0,kcaturglnk)
global time
result = kcaturglnk*UTase*GlnKUMP2/((1+Kurgln/GLN)*(Kurglnkump+(1+UMP/Kurump)*(GlnKUMP+GlnKUMP2+GlnKUMP3)))/default0;
return


