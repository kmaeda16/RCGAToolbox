function [output] = BIOMD0000000051_url_mod_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chassagnole2002_Carbon_Metabolism
% Generated: 08-Jul-2019 12:18:02
% 
% [output] = BIOMD0000000051_url_mod_odefun() => output = initial conditions in column vector
% [output] = BIOMD0000000051_url_mod_odefun('states') => output = state names in cell-array
% [output] = BIOMD0000000051_url_mod_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = BIOMD0000000051_url_mod_odefun('parameters') => output = parameter names in cell-array
% [output] = BIOMD0000000051_url_mod_odefun('parametervalues') => output = parameter values in column vector
% [output] = BIOMD0000000051_url_mod_odefun('variablenames') => output = variable names in cell-array
% [output] = BIOMD0000000051_url_mod_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = BIOMD0000000051_url_mod_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): cpep
% statevector(2): cglcex
% statevector(3): cg6p
% statevector(4): cpyr
% statevector(5): cf6p
% statevector(6): cg1p
% statevector(7): cpg
% statevector(8): cfdp
% statevector(9): csed7p
% statevector(10): cgap
% statevector(11): ce4p
% statevector(12): cxyl5p
% statevector(13): crib5p
% statevector(14): cdhap
% statevector(15): cpgp
% statevector(16): cpg3
% statevector(17): cpg2
% statevector(18): cribu5p
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [2.67, 2, 3.48, 2.67, 0.6, 0.653, 0.808, 0.272, 0.276, 0.218, ...
		0.098, 0.138, 0.398, 0.167, 0.008, 2.13, 0.399, 0.111];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'cpep', 'cglcex', 'cg6p', 'cpyr', 'cf6p', 'cg1p', 'cpg', 'cfdp', 'csed7p', 'cgap', ...
			'ce4p', 'cxyl5p', 'crib5p', 'cdhap', 'cpgp', 'cpg3', 'cpg2', 'cribu5p'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'vPTS_rmaxPTS', 'vPTS_KPTSa1', 'vPTS_KPTSa2', 'vPTS_KPTSa3', 'vPTS_nPTSg6p', 'vPTS_KPTSg6p', 'vPGI_rmaxPGI', 'vPGI_KPGIeq', 'vPGI_KPGIg6p', 'vPGI_KPGIf6p', ...
			'vPGI_KPGIf6ppginh', 'vPGI_KPGIg6ppginh', 'vPGM_rmaxPGM', 'vPGM_KPGMeq', 'vPGM_KPGMg6p', 'vPGM_KPGMg1p', 'vG6PDH_rmaxG6PDH', 'vG6PDH_KG6PDHg6p', 'vG6PDH_KG6PDHnadphg6pinh', 'vG6PDH_KG6PDHnadp', ...
			'vG6PDH_KG6PDHnadphnadpinh', 'vPFK_rmaxPFK', 'vPFK_KPFKatps', 'vPFK_KPFKadpc', 'vPFK_KPFKf6ps', 'vPFK_KPFKpep', 'vPFK_KPFKadpb', 'vPFK_KPFKampb', 'vPFK_KPFKadpa', 'vPFK_KPFKampa', ...
			'vPFK_LPFK', 'vPFK_nPFK', 'vTA_rmaxTA', 'vTA_KTAeq', 'vTKA_rmaxTKa', 'vTKA_KTKaeq', 'vTKB_rmaxTKb', 'vTKB_KTKbeq', 'vMURSyNTH_rmaxMurSynth', 'vALDO_rmaxALDO', ...
			'vALDO_kALDOeq', 'vALDO_kALDOfdp', 'vALDO_kALDOgap', 'vALDO_VALDOblf', 'vALDO_kALDOdhap', 'vALDO_kALDOgapinh', 'vGAPDH_rmaxGAPDH', 'vGAPDH_KGAPDHeq', 'vGAPDH_KGAPDHgap', 'vGAPDH_KGAPDHpgp', ...
			'vGAPDH_KGAPDHnad', 'vGAPDH_KGAPDHnadh', 'vTIS_rmaxTIS', 'vTIS_kTISeq', 'vTIS_kTISdhap', 'vTIS_kTISgap', 'vTRPSYNTH_rmaxTrpSynth', 'vG3PDH_rmaxG3PDH', 'vG3PDH_KG3PDHdhap', 'vPGK_rmaxPGK', ...
			'vPGK_KPGKeq', 'vPGK_KPGKadp', 'vPGK_KPGKatp', 'vPGK_KPGKpgp', 'vPGK_KPGKpg3', 'vsersynth_rmaxSerSynth', 'vsersynth_KSerSynthpg3', 'vrpGluMu_rmaxPGluMu', 'vrpGluMu_KPGluMueq', 'vrpGluMu_KPGluMupg3', ...
			'vrpGluMu_KPGluMupg2', 'vENO_rmaxENO', 'vENO_KENOeq', 'vENO_KENOpg2', 'vENO_KENOpep', 'vPK_rmaxPK', 'vPK_KPKpep', 'vPK_nPK', 'vPK_LPK', 'vPK_KPKatp', ...
			'vPK_KPKfdp', 'vPK_KPKamp', 'vPK_KPKadp', 'vpepCxylase_rmaxpepCxylase', 'vpepCxylase_KpepCxylasefdp', 'vpepCxylase_npepCxylasefdp', 'vpepCxylase_KpepCxylasepep', 'vSynth1_rmaxSynth1', 'vSynth1_KSynth1pep', 'vSynth2_rmaxSynth2', ...
			'vSynth2_KSynth2pyr', 'vDAHPS_rmaxDAHPS', 'vDAHPS_nDAHPSe4p', 'vDAHPS_nDAHPSpep', 'vDAHPS_KDAHPSe4p', 'vDAHPS_KDAHPSpep', 'vPDH_rmaxPDH', 'vPDH_nPDH', 'vPDH_KPDHpyr', 'vMethSynth_rmaxMetSynth', ...
			'vPGDH_rmaxPGDH', 'vPGDH_KPGDHpg', 'vPGDH_KPGDHnadp', 'vPGDH_KPGDHnadphinh', 'vPGDH_KPGDHatpinh', 'vR5PI_rmaxR5PI', 'vR5PI_KR5PIeq', 'vRu5P_rmaxRu5P', 'vRu5P_KRu5Peq', 'vPPK_rmaxRPPK', ...
			'vPPK_KRPPKrib5p', 'vG1PAT_rmaxG1PAT', 'vG1PAT_KG1PATfdp', 'vG1PAT_nG1PATfdp', 'vG1PAT_KG1PATatp', 'vG1PAT_KG1PATg1p', 'vG6P_mu', 'vf6P_mu', 'vfdP_mu', 'vGAP_mu', ...
			'vDHAP_mu', 'vPGP_mu', 'vPG3_mu', 'vpg2_mu', 'vPEP_mu', 'vRibu5p_mu', 'vRIB5P_mu', 'vXYL5P_mu', 'vSED7P_mu', 'vpyr_mu', ...
			'vPG_mu', 'vE4P_mu', 'vGLP_mu', 'vEXTER_Dil', 'vEXTER_cfeed', 'extracellular', 'cytosol'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [7829.78, 3082.3, 0.01, 245.3, 3.66, 2.15, 650.988, 0.1725, 2.9, 0.266, ...
			0.2, 0.2, 0.839824, 0.196, 1.038, 0.0136, 1.3802, 14.4, 6.43, 0.0246, ...
			0.01, 1840.58, 0.123, 4.14, 0.325, 3.26, 3.89, 3.2, 128, 19.1, ...
			5.62907e+06, 11.1, 10.8716, 1.05, 9.47338, 1.2, 86.5586, 10, 0.00043711, 17.4146, ...
			0.144, 1.75, 0.088, 2, 0.088, 0.6, 921.594, 0.63, 0.683, 1.04e-05, ...
			0.252, 1.09, 68.6747, 1.39, 2.8, 0.3, 0.001037, 0.0116204, 1, 3021.77, ...
			1934.4, 0.185, 0.653, 0.0468, 0.473, 0.0257121, 1, 89.0497, 0.188, 0.2, ...
			0.369, 330.448, 6.73, 0.1, 0.135, 0.0611315, 0.31, 4, 1000, 22.5, ...
			0.19, 0.2, 0.26, 0.107021, 0.7, 4.21, 4.07, 0.019539, 1, 0.0736186, ...
			1, 0.107953, 2.6, 2.2, 0.035, 0.0053, 6.05953, 3.68, 1159, 0.0022627, ...
			16.2324, 37.5, 0.0506, 0.0138, 208, 4.83841, 4, 6.73903, 1.4, 0.0129005, ...
			0.1, 0.00752546, 0.119, 1.2, 4.42, 3.2, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, ...
			2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, ...
			2.78e-05, 2.78e-05, 2.78e-05, 2.78e-05, 110.96, 1, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'catp', 'cadp', 'camp', 'cnadph', 'cnadp', 'cnadh', 'cnad'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'4.27-4.163*(time/(0.657+1.43*time+0.0364*power(time,2)))', '0.582+1.73*power(2.731,-0.15*time)*(0.12*time+0.000214*power(time,3))', '0.123+7.25*(time/(7.25+1.47*time+0.17*power(time,2)))+1.073/(1.29+8.05*time)', '0.062+0.332*power(2.718,-0.464*time)*(0.0166*power(time,1.58)+0.000166*power(time,4.73)+0.1312*power(10,-9)*power(time,7.89)+0.1362*power(10,-12)*power(time,11)+0.1233*power(10,-15)*power(time,14.2))', '0.159-0.00554*(time/(2.8-0.271*time+0.01*power(time,2)))+0.182/(4.82+0.526*time)', '0.0934+0.00111*power(2.371,-0.123*time)*(0.844*time+0.104*power(time,3))', '1.314+1.314*power(2.73,-0.0435*time-0.342)-(time+7.871)*(power(2.73,-0.0218*time-0.171)/(8.481+time))'};
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
	if length(parameterValuesNew) ~= 137,
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
cpep = statevector(1);
cglcex = statevector(2);
cg6p = statevector(3);
cpyr = statevector(4);
cf6p = statevector(5);
cg1p = statevector(6);
cpg = statevector(7);
cfdp = statevector(8);
csed7p = statevector(9);
cgap = statevector(10);
ce4p = statevector(11);
cxyl5p = statevector(12);
crib5p = statevector(13);
cdhap = statevector(14);
cpgp = statevector(15);
cpg3 = statevector(16);
cpg2 = statevector(17);
cribu5p = statevector(18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	vPTS_rmaxPTS = 7829.78;
	vPTS_KPTSa1 = 3082.3;
	vPTS_KPTSa2 = 0.01;
	vPTS_KPTSa3 = 245.3;
	vPTS_nPTSg6p = 3.66;
	vPTS_KPTSg6p = 2.15;
	vPGI_rmaxPGI = 650.988;
	vPGI_KPGIeq = 0.1725;
	vPGI_KPGIg6p = 2.9;
	vPGI_KPGIf6p = 0.266;
	vPGI_KPGIf6ppginh = 0.2;
	vPGI_KPGIg6ppginh = 0.2;
	vPGM_rmaxPGM = 0.839824;
	vPGM_KPGMeq = 0.196;
	vPGM_KPGMg6p = 1.038;
	vPGM_KPGMg1p = 0.0136;
	vG6PDH_rmaxG6PDH = 1.3802;
	vG6PDH_KG6PDHg6p = 14.4;
	vG6PDH_KG6PDHnadphg6pinh = 6.43;
	vG6PDH_KG6PDHnadp = 0.0246;
	vG6PDH_KG6PDHnadphnadpinh = 0.01;
	vPFK_rmaxPFK = 1840.58;
	vPFK_KPFKatps = 0.123;
	vPFK_KPFKadpc = 4.14;
	vPFK_KPFKf6ps = 0.325;
	vPFK_KPFKpep = 3.26;
	vPFK_KPFKadpb = 3.89;
	vPFK_KPFKampb = 3.2;
	vPFK_KPFKadpa = 128;
	vPFK_KPFKampa = 19.1;
	vPFK_LPFK = 5.62907e+06;
	vPFK_nPFK = 11.1;
	vTA_rmaxTA = 10.8716;
	vTA_KTAeq = 1.05;
	vTKA_rmaxTKa = 9.47338;
	vTKA_KTKaeq = 1.2;
	vTKB_rmaxTKb = 86.5586;
	vTKB_KTKbeq = 10;
	vMURSyNTH_rmaxMurSynth = 0.00043711;
	vALDO_rmaxALDO = 17.4146;
	vALDO_kALDOeq = 0.144;
	vALDO_kALDOfdp = 1.75;
	vALDO_kALDOgap = 0.088;
	vALDO_VALDOblf = 2;
	vALDO_kALDOdhap = 0.088;
	vALDO_kALDOgapinh = 0.6;
	vGAPDH_rmaxGAPDH = 921.594;
	vGAPDH_KGAPDHeq = 0.63;
	vGAPDH_KGAPDHgap = 0.683;
	vGAPDH_KGAPDHpgp = 1.04e-05;
	vGAPDH_KGAPDHnad = 0.252;
	vGAPDH_KGAPDHnadh = 1.09;
	vTIS_rmaxTIS = 68.6747;
	vTIS_kTISeq = 1.39;
	vTIS_kTISdhap = 2.8;
	vTIS_kTISgap = 0.3;
	vTRPSYNTH_rmaxTrpSynth = 0.001037;
	vG3PDH_rmaxG3PDH = 0.0116204;
	vG3PDH_KG3PDHdhap = 1;
	vPGK_rmaxPGK = 3021.77;
	vPGK_KPGKeq = 1934.4;
	vPGK_KPGKadp = 0.185;
	vPGK_KPGKatp = 0.653;
	vPGK_KPGKpgp = 0.0468;
	vPGK_KPGKpg3 = 0.473;
	vsersynth_rmaxSerSynth = 0.0257121;
	vsersynth_KSerSynthpg3 = 1;
	vrpGluMu_rmaxPGluMu = 89.0497;
	vrpGluMu_KPGluMueq = 0.188;
	vrpGluMu_KPGluMupg3 = 0.2;
	vrpGluMu_KPGluMupg2 = 0.369;
	vENO_rmaxENO = 330.448;
	vENO_KENOeq = 6.73;
	vENO_KENOpg2 = 0.1;
	vENO_KENOpep = 0.135;
	vPK_rmaxPK = 0.0611315;
	vPK_KPKpep = 0.31;
	vPK_nPK = 4;
	vPK_LPK = 1000;
	vPK_KPKatp = 22.5;
	vPK_KPKfdp = 0.19;
	vPK_KPKamp = 0.2;
	vPK_KPKadp = 0.26;
	vpepCxylase_rmaxpepCxylase = 0.107021;
	vpepCxylase_KpepCxylasefdp = 0.7;
	vpepCxylase_npepCxylasefdp = 4.21;
	vpepCxylase_KpepCxylasepep = 4.07;
	vSynth1_rmaxSynth1 = 0.019539;
	vSynth1_KSynth1pep = 1;
	vSynth2_rmaxSynth2 = 0.0736186;
	vSynth2_KSynth2pyr = 1;
	vDAHPS_rmaxDAHPS = 0.107953;
	vDAHPS_nDAHPSe4p = 2.6;
	vDAHPS_nDAHPSpep = 2.2;
	vDAHPS_KDAHPSe4p = 0.035;
	vDAHPS_KDAHPSpep = 0.0053;
	vPDH_rmaxPDH = 6.05953;
	vPDH_nPDH = 3.68;
	vPDH_KPDHpyr = 1159;
	vMethSynth_rmaxMetSynth = 0.0022627;
	vPGDH_rmaxPGDH = 16.2324;
	vPGDH_KPGDHpg = 37.5;
	vPGDH_KPGDHnadp = 0.0506;
	vPGDH_KPGDHnadphinh = 0.0138;
	vPGDH_KPGDHatpinh = 208;
	vR5PI_rmaxR5PI = 4.83841;
	vR5PI_KR5PIeq = 4;
	vRu5P_rmaxRu5P = 6.73903;
	vRu5P_KRu5Peq = 1.4;
	vPPK_rmaxRPPK = 0.0129005;
	vPPK_KRPPKrib5p = 0.1;
	vG1PAT_rmaxG1PAT = 0.00752546;
	vG1PAT_KG1PATfdp = 0.119;
	vG1PAT_nG1PATfdp = 1.2;
	vG1PAT_KG1PATatp = 4.42;
	vG1PAT_KG1PATg1p = 3.2;
	vG6P_mu = 2.78e-05;
	vf6P_mu = 2.78e-05;
	vfdP_mu = 2.78e-05;
	vGAP_mu = 2.78e-05;
	vDHAP_mu = 2.78e-05;
	vPGP_mu = 2.78e-05;
	vPG3_mu = 2.78e-05;
	vpg2_mu = 2.78e-05;
	vPEP_mu = 2.78e-05;
	vRibu5p_mu = 2.78e-05;
	vRIB5P_mu = 2.78e-05;
	vXYL5P_mu = 2.78e-05;
	vSED7P_mu = 2.78e-05;
	vpyr_mu = 2.78e-05;
	vPG_mu = 2.78e-05;
	vE4P_mu = 2.78e-05;
	vGLP_mu = 2.78e-05;
	vEXTER_Dil = 2.78e-05;
	vEXTER_cfeed = 110.96;
	extracellular = 1;
	cytosol = 1;
else
	vPTS_rmaxPTS = parameterValuesNew(1);
	vPTS_KPTSa1 = parameterValuesNew(2);
	vPTS_KPTSa2 = parameterValuesNew(3);
	vPTS_KPTSa3 = parameterValuesNew(4);
	vPTS_nPTSg6p = parameterValuesNew(5);
	vPTS_KPTSg6p = parameterValuesNew(6);
	vPGI_rmaxPGI = parameterValuesNew(7);
	vPGI_KPGIeq = parameterValuesNew(8);
	vPGI_KPGIg6p = parameterValuesNew(9);
	vPGI_KPGIf6p = parameterValuesNew(10);
	vPGI_KPGIf6ppginh = parameterValuesNew(11);
	vPGI_KPGIg6ppginh = parameterValuesNew(12);
	vPGM_rmaxPGM = parameterValuesNew(13);
	vPGM_KPGMeq = parameterValuesNew(14);
	vPGM_KPGMg6p = parameterValuesNew(15);
	vPGM_KPGMg1p = parameterValuesNew(16);
	vG6PDH_rmaxG6PDH = parameterValuesNew(17);
	vG6PDH_KG6PDHg6p = parameterValuesNew(18);
	vG6PDH_KG6PDHnadphg6pinh = parameterValuesNew(19);
	vG6PDH_KG6PDHnadp = parameterValuesNew(20);
	vG6PDH_KG6PDHnadphnadpinh = parameterValuesNew(21);
	vPFK_rmaxPFK = parameterValuesNew(22);
	vPFK_KPFKatps = parameterValuesNew(23);
	vPFK_KPFKadpc = parameterValuesNew(24);
	vPFK_KPFKf6ps = parameterValuesNew(25);
	vPFK_KPFKpep = parameterValuesNew(26);
	vPFK_KPFKadpb = parameterValuesNew(27);
	vPFK_KPFKampb = parameterValuesNew(28);
	vPFK_KPFKadpa = parameterValuesNew(29);
	vPFK_KPFKampa = parameterValuesNew(30);
	vPFK_LPFK = parameterValuesNew(31);
	vPFK_nPFK = parameterValuesNew(32);
	vTA_rmaxTA = parameterValuesNew(33);
	vTA_KTAeq = parameterValuesNew(34);
	vTKA_rmaxTKa = parameterValuesNew(35);
	vTKA_KTKaeq = parameterValuesNew(36);
	vTKB_rmaxTKb = parameterValuesNew(37);
	vTKB_KTKbeq = parameterValuesNew(38);
	vMURSyNTH_rmaxMurSynth = parameterValuesNew(39);
	vALDO_rmaxALDO = parameterValuesNew(40);
	vALDO_kALDOeq = parameterValuesNew(41);
	vALDO_kALDOfdp = parameterValuesNew(42);
	vALDO_kALDOgap = parameterValuesNew(43);
	vALDO_VALDOblf = parameterValuesNew(44);
	vALDO_kALDOdhap = parameterValuesNew(45);
	vALDO_kALDOgapinh = parameterValuesNew(46);
	vGAPDH_rmaxGAPDH = parameterValuesNew(47);
	vGAPDH_KGAPDHeq = parameterValuesNew(48);
	vGAPDH_KGAPDHgap = parameterValuesNew(49);
	vGAPDH_KGAPDHpgp = parameterValuesNew(50);
	vGAPDH_KGAPDHnad = parameterValuesNew(51);
	vGAPDH_KGAPDHnadh = parameterValuesNew(52);
	vTIS_rmaxTIS = parameterValuesNew(53);
	vTIS_kTISeq = parameterValuesNew(54);
	vTIS_kTISdhap = parameterValuesNew(55);
	vTIS_kTISgap = parameterValuesNew(56);
	vTRPSYNTH_rmaxTrpSynth = parameterValuesNew(57);
	vG3PDH_rmaxG3PDH = parameterValuesNew(58);
	vG3PDH_KG3PDHdhap = parameterValuesNew(59);
	vPGK_rmaxPGK = parameterValuesNew(60);
	vPGK_KPGKeq = parameterValuesNew(61);
	vPGK_KPGKadp = parameterValuesNew(62);
	vPGK_KPGKatp = parameterValuesNew(63);
	vPGK_KPGKpgp = parameterValuesNew(64);
	vPGK_KPGKpg3 = parameterValuesNew(65);
	vsersynth_rmaxSerSynth = parameterValuesNew(66);
	vsersynth_KSerSynthpg3 = parameterValuesNew(67);
	vrpGluMu_rmaxPGluMu = parameterValuesNew(68);
	vrpGluMu_KPGluMueq = parameterValuesNew(69);
	vrpGluMu_KPGluMupg3 = parameterValuesNew(70);
	vrpGluMu_KPGluMupg2 = parameterValuesNew(71);
	vENO_rmaxENO = parameterValuesNew(72);
	vENO_KENOeq = parameterValuesNew(73);
	vENO_KENOpg2 = parameterValuesNew(74);
	vENO_KENOpep = parameterValuesNew(75);
	vPK_rmaxPK = parameterValuesNew(76);
	vPK_KPKpep = parameterValuesNew(77);
	vPK_nPK = parameterValuesNew(78);
	vPK_LPK = parameterValuesNew(79);
	vPK_KPKatp = parameterValuesNew(80);
	vPK_KPKfdp = parameterValuesNew(81);
	vPK_KPKamp = parameterValuesNew(82);
	vPK_KPKadp = parameterValuesNew(83);
	vpepCxylase_rmaxpepCxylase = parameterValuesNew(84);
	vpepCxylase_KpepCxylasefdp = parameterValuesNew(85);
	vpepCxylase_npepCxylasefdp = parameterValuesNew(86);
	vpepCxylase_KpepCxylasepep = parameterValuesNew(87);
	vSynth1_rmaxSynth1 = parameterValuesNew(88);
	vSynth1_KSynth1pep = parameterValuesNew(89);
	vSynth2_rmaxSynth2 = parameterValuesNew(90);
	vSynth2_KSynth2pyr = parameterValuesNew(91);
	vDAHPS_rmaxDAHPS = parameterValuesNew(92);
	vDAHPS_nDAHPSe4p = parameterValuesNew(93);
	vDAHPS_nDAHPSpep = parameterValuesNew(94);
	vDAHPS_KDAHPSe4p = parameterValuesNew(95);
	vDAHPS_KDAHPSpep = parameterValuesNew(96);
	vPDH_rmaxPDH = parameterValuesNew(97);
	vPDH_nPDH = parameterValuesNew(98);
	vPDH_KPDHpyr = parameterValuesNew(99);
	vMethSynth_rmaxMetSynth = parameterValuesNew(100);
	vPGDH_rmaxPGDH = parameterValuesNew(101);
	vPGDH_KPGDHpg = parameterValuesNew(102);
	vPGDH_KPGDHnadp = parameterValuesNew(103);
	vPGDH_KPGDHnadphinh = parameterValuesNew(104);
	vPGDH_KPGDHatpinh = parameterValuesNew(105);
	vR5PI_rmaxR5PI = parameterValuesNew(106);
	vR5PI_KR5PIeq = parameterValuesNew(107);
	vRu5P_rmaxRu5P = parameterValuesNew(108);
	vRu5P_KRu5Peq = parameterValuesNew(109);
	vPPK_rmaxRPPK = parameterValuesNew(110);
	vPPK_KRPPKrib5p = parameterValuesNew(111);
	vG1PAT_rmaxG1PAT = parameterValuesNew(112);
	vG1PAT_KG1PATfdp = parameterValuesNew(113);
	vG1PAT_nG1PATfdp = parameterValuesNew(114);
	vG1PAT_KG1PATatp = parameterValuesNew(115);
	vG1PAT_KG1PATg1p = parameterValuesNew(116);
	vG6P_mu = parameterValuesNew(117);
	vf6P_mu = parameterValuesNew(118);
	vfdP_mu = parameterValuesNew(119);
	vGAP_mu = parameterValuesNew(120);
	vDHAP_mu = parameterValuesNew(121);
	vPGP_mu = parameterValuesNew(122);
	vPG3_mu = parameterValuesNew(123);
	vpg2_mu = parameterValuesNew(124);
	vPEP_mu = parameterValuesNew(125);
	vRibu5p_mu = parameterValuesNew(126);
	vRIB5P_mu = parameterValuesNew(127);
	vXYL5P_mu = parameterValuesNew(128);
	vSED7P_mu = parameterValuesNew(129);
	vpyr_mu = parameterValuesNew(130);
	vPG_mu = parameterValuesNew(131);
	vE4P_mu = parameterValuesNew(132);
	vGLP_mu = parameterValuesNew(133);
	vEXTER_Dil = parameterValuesNew(134);
	vEXTER_cfeed = parameterValuesNew(135);
	extracellular = parameterValuesNew(136);
	cytosol = parameterValuesNew(137);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catp = 4.27-4.163*(time/(0.657+1.43*time+0.0364*power(time,2)));
cadp = 0.582+1.73*power(2.731,-0.15*time)*(0.12*time+0.000214*power(time,3));
camp = 0.123+7.25*(time/(7.25+1.47*time+0.17*power(time,2)))+1.073/(1.29+8.05*time);
cnadph = 0.062+0.332*power(2.718,-0.464*time)*(0.0166*power(time,1.58)+0.000166*power(time,4.73)+0.1312*power(10,-9)*power(time,7.89)+0.1362*power(10,-12)*power(time,11)+0.1233*power(10,-15)*power(time,14.2));
cnadp = 0.159-0.00554*(time/(2.8-0.271*time+0.01*power(time,2)))+0.182/(4.82+0.526*time);
cnadh = 0.0934+0.00111*power(2.371,-0.123*time)*(0.844*time+0.104*power(time,3));
cnad = 1.314+1.314*power(2.73,-0.0435*time-0.342)-(time+7.871)*(power(2.73,-0.0218*time-0.171)/(8.481+time));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vPTS = extracellular * vPTS_rmaxPTS * cglcex * (cpep / cpyr) / ((vPTS_KPTSa1 + vPTS_KPTSa2 * (cpep / cpyr) + vPTS_KPTSa3 * cglcex + cglcex * (cpep / cpyr)) * (1 + power(cg6p, vPTS_nPTSg6p) / vPTS_KPTSg6p));
vPGI = cytosol * vPGI_rmaxPGI * (cg6p - cf6p / vPGI_KPGIeq) / (vPGI_KPGIg6p * (1 + cf6p / (vPGI_KPGIf6p * (1 + cpg / vPGI_KPGIf6ppginh)) + cpg / vPGI_KPGIg6ppginh) + cg6p);
vPGM = cytosol * vPGM_rmaxPGM * (cg6p - cg1p / vPGM_KPGMeq) / (vPGM_KPGMg6p * (1 + cg1p / vPGM_KPGMg1p) + cg6p);
vG6PDH = cytosol * vG6PDH_rmaxG6PDH * cg6p * cnadp / ((cg6p + vG6PDH_KG6PDHg6p) * (1 + cnadph / vG6PDH_KG6PDHnadphg6pinh) * (vG6PDH_KG6PDHnadp * (1 + cnadph / vG6PDH_KG6PDHnadphnadpinh) + cnadp));
vPFK = cytosol * vPFK_rmaxPFK * catp * cf6p / ((catp + vPFK_KPFKatps * (1 + cadp / vPFK_KPFKadpc)) * (cf6p + vPFK_KPFKf6ps * (1 + cpep / vPFK_KPFKpep + cadp / vPFK_KPFKadpb + camp / vPFK_KPFKampb) / (1 + cadp / vPFK_KPFKadpa + camp / vPFK_KPFKampa)) * (1 + vPFK_LPFK / power(1 + cf6p * (1 + cadp / vPFK_KPFKadpa + camp / vPFK_KPFKampa) / (vPFK_KPFKf6ps * (1 + cpep / vPFK_KPFKpep + cadp / vPFK_KPFKadpb + camp / vPFK_KPFKampb)), vPFK_nPFK)));
vTA = cytosol * vTA_rmaxTA * (cgap * csed7p - ce4p * cf6p / vTA_KTAeq);
vTKA = cytosol * vTKA_rmaxTKa * (crib5p * cxyl5p - csed7p * cgap / vTKA_KTKaeq);
vTKB = cytosol * vTKB_rmaxTKb * (cxyl5p * ce4p - cf6p * cgap / vTKB_KTKbeq);
vMURSyNTH = cytosol * vMURSyNTH_rmaxMurSynth;
vALDO = cytosol * vALDO_rmaxALDO * (cfdp - cgap * cdhap / vALDO_kALDOeq) / (vALDO_kALDOfdp + cfdp + vALDO_kALDOgap * cdhap / (vALDO_kALDOeq * vALDO_VALDOblf) + vALDO_kALDOdhap * cgap / (vALDO_kALDOeq * vALDO_VALDOblf) + cfdp * cgap / vALDO_kALDOgapinh + cgap * cdhap / (vALDO_VALDOblf * vALDO_kALDOeq));
vGAPDH = cytosol * vGAPDH_rmaxGAPDH * (cgap * cnad - cpgp * cnadh / vGAPDH_KGAPDHeq) / ((vGAPDH_KGAPDHgap * (1 + cpgp / vGAPDH_KGAPDHpgp) + cgap) * (vGAPDH_KGAPDHnad * (1 + cnadh / vGAPDH_KGAPDHnadh) + cnad));
vTIS = cytosol * vTIS_rmaxTIS * (cdhap - cgap / vTIS_kTISeq) / (vTIS_kTISdhap * (1 + cgap / vTIS_kTISgap) + cdhap);
vTRPSYNTH = cytosol * vTRPSYNTH_rmaxTrpSynth;
vG3PDH = cytosol * vG3PDH_rmaxG3PDH * cdhap / (vG3PDH_KG3PDHdhap + cdhap);
vPGK = cytosol * vPGK_rmaxPGK * (cadp * cpgp - catp * cpg3 / vPGK_KPGKeq) / ((vPGK_KPGKadp * (1 + catp / vPGK_KPGKatp) + cadp) * (vPGK_KPGKpgp * (1 + cpg3 / vPGK_KPGKpg3) + cpgp));
vsersynth = cytosol * vsersynth_rmaxSerSynth * cpg3 / (vsersynth_KSerSynthpg3 + cpg3);
vrpGluMu = cytosol * vrpGluMu_rmaxPGluMu * (cpg3 - cpg2 / vrpGluMu_KPGluMueq) / (vrpGluMu_KPGluMupg3 * (1 + cpg2 / vrpGluMu_KPGluMupg2) + cpg3);
vENO = cytosol * vENO_rmaxENO * (cpg2 - cpep / vENO_KENOeq) / (vENO_KENOpg2 * (1 + cpep / vENO_KENOpep) + cpg2);
vPK = cytosol * vPK_rmaxPK * cpep * power(cpep / vPK_KPKpep + 1, vPK_nPK - 1) * cadp / (vPK_KPKpep * (vPK_LPK * power((1 + catp / vPK_KPKatp) / (cfdp / vPK_KPKfdp + camp / vPK_KPKamp + 1), vPK_nPK) + power(cpep / vPK_KPKpep + 1, vPK_nPK)) * (cadp + vPK_KPKadp));
vpepCxylase = cytosol * vpepCxylase_rmaxpepCxylase * cpep * (1 + power(cfdp / vpepCxylase_KpepCxylasefdp, vpepCxylase_npepCxylasefdp)) / (vpepCxylase_KpepCxylasepep + cpep);
vSynth1 = cytosol * vSynth1_rmaxSynth1 * cpep / (vSynth1_KSynth1pep + cpep);
vSynth2 = cytosol * vSynth2_rmaxSynth2 * cpyr / (vSynth2_KSynth2pyr + cpyr);
vDAHPS = cytosol * vDAHPS_rmaxDAHPS * power(ce4p, vDAHPS_nDAHPSe4p) * power(cpep, vDAHPS_nDAHPSpep) / ((vDAHPS_KDAHPSe4p + power(ce4p, vDAHPS_nDAHPSe4p)) * (vDAHPS_KDAHPSpep + power(cpep, vDAHPS_nDAHPSpep)));
vPDH = cytosol * vPDH_rmaxPDH * power(cpyr, vPDH_nPDH) / (vPDH_KPDHpyr + power(cpyr, vPDH_nPDH));
vMethSynth = cytosol * vMethSynth_rmaxMetSynth;
vPGDH = cytosol * vPGDH_rmaxPGDH * cpg * cnadp / ((cpg + vPGDH_KPGDHpg) * (cnadp + vPGDH_KPGDHnadp * (1 + cnadph / vPGDH_KPGDHnadphinh) * (1 + catp / vPGDH_KPGDHatpinh)));
vR5PI = cytosol * vR5PI_rmaxR5PI * (cribu5p - crib5p / vR5PI_KR5PIeq);
vRu5P = cytosol * vRu5P_rmaxRu5P * (cribu5p - cxyl5p / vRu5P_KRu5Peq);
vPPK = cytosol * vPPK_rmaxRPPK * crib5p / (vPPK_KRPPKrib5p + crib5p);
vG1PAT = cytosol * vG1PAT_rmaxG1PAT * cg1p * catp * (1 + power(cfdp / vG1PAT_KG1PATfdp, vG1PAT_nG1PATfdp)) / ((vG1PAT_KG1PATatp + catp) * (vG1PAT_KG1PATg1p + cg1p));
vG6P = cytosol * vG6P_mu * cg6p;
vf6P = cytosol * vf6P_mu * cf6p;
vfdP = cytosol * vfdP_mu * cfdp;
vGAP = cytosol * vGAP_mu * cgap;
vDHAP = cytosol * vDHAP_mu * cdhap;
vPGP = cytosol * vPGP_mu * cpgp;
vPG3 = cytosol * vPG3_mu * cpg3;
vpg2 = cytosol * vpg2_mu * cpg2;
vPEP = cytosol * vPEP_mu * cpep;
vRibu5p = cytosol * vRibu5p_mu * cribu5p;
vRIB5P = cytosol * vRIB5P_mu * crib5p;
vXYL5P = cytosol * vXYL5P_mu * cxyl5p;
vSED7P = cytosol * vSED7P_mu * csed7p;
vpyr = cytosol * vpyr_mu * cpyr;
vPG = cytosol * vPG_mu * cpg;
vE4P = cytosol * vE4P_mu * ce4p;
vGLP = cytosol * vGLP_mu * cg1p;
vEXTER = extracellular * vEXTER_Dil * (vEXTER_cfeed - cglcex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cpep_dot = (-65*vPTS+vENO-vPK-vpepCxylase-vSynth1-vDAHPS-vPEP)/cytosol;
cglcex_dot = (-vPTS+vEXTER)/extracellular;
cg6p_dot = (+65*vPTS-vPGI-vPGM-vG6PDH-vG6P)/cytosol;
cpyr_dot = (+65*vPTS+vTRPSYNTH+vPK-vSynth2-vPDH+vMethSynth-vpyr)/cytosol;
cf6p_dot = (+vPGI-vPFK+vTA+vTKB-2*vMURSyNTH-vf6P)/cytosol;
cg1p_dot = (+vPGM-vG1PAT-vGLP)/cytosol;
cpg_dot = (+vG6PDH-vPGDH-vPG)/cytosol;
cfdp_dot = (+vPFK-vALDO-vfdP)/cytosol;
csed7p_dot = (-vTA+vTKA-vSED7P)/cytosol;
cgap_dot = (-vTA+vTKA+vTKB+vALDO-vGAPDH+vTIS+vTRPSYNTH-vGAP)/cytosol;
ce4p_dot = (+vTA-vTKB-vDAHPS-vE4P)/cytosol;
cxyl5p_dot = (-vTKA-vTKB+vRu5P-vXYL5P)/cytosol;
crib5p_dot = (-vTKA+vR5PI-vPPK-vRIB5P)/cytosol;
cdhap_dot = (+vALDO-vTIS-vG3PDH-vDHAP)/cytosol;
cpgp_dot = (+vGAPDH-vPGK-vPGP)/cytosol;
cpg3_dot = (+vPGK-vsersynth-vrpGluMu-vPG3)/cytosol;
cpg2_dot = (+vrpGluMu-vENO-vpg2)/cytosol;
cribu5p_dot = (+vPGDH-vR5PI-vRu5P-vRibu5p)/cytosol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = cpep_dot;
output(2) = cglcex_dot;
output(3) = cg6p_dot;
output(4) = cpyr_dot;
output(5) = cf6p_dot;
output(6) = cg1p_dot;
output(7) = cpg_dot;
output(8) = cfdp_dot;
output(9) = csed7p_dot;
output(10) = cgap_dot;
output(11) = ce4p_dot;
output(12) = cxyl5p_dot;
output(13) = crib5p_dot;
output(14) = cdhap_dot;
output(15) = cpgp_dot;
output(16) = cpg3_dot;
output(17) = cpg2_dot;
output(18) = cribu5p_dot;
% return a column vector 
output = output(:);
return


