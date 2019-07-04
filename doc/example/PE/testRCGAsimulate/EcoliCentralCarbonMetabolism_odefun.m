function [output] = EcoliCentralCarbonMetabolism_odefun(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EcoliCentralCarbonMetabolism
% Generated: 04-Jul-2019 18:58:56
% 
% [output] = EcoliCentralCarbonMetabolism_odefun() => output = initial conditions in column vector
% [output] = EcoliCentralCarbonMetabolism_odefun('states') => output = state names in cell-array
% [output] = EcoliCentralCarbonMetabolism_odefun('algebraic') => output = algebraic variable names in cell-array
% [output] = EcoliCentralCarbonMetabolism_odefun('parameters') => output = parameter names in cell-array
% [output] = EcoliCentralCarbonMetabolism_odefun('parametervalues') => output = parameter values in column vector
% [output] = EcoliCentralCarbonMetabolism_odefun('variablenames') => output = variable names in cell-array
% [output] = EcoliCentralCarbonMetabolism_odefun('variableformulas') => output = variable formulas in cell-array
% [output] = EcoliCentralCarbonMetabolism_odefun(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): ACEex
% statevector(2): AcCoA
% statevector(3): AcP
% statevector(4): AceK
% statevector(5): Acs
% statevector(6): CS
% statevector(7): E4P
% statevector(8): EIIAP
% statevector(9): F6P
% statevector(10): FBP
% statevector(11): FUM
% statevector(12): Fba
% statevector(13): Fbp
% statevector(14): Fum
% statevector(15): G6P
% statevector(16): GAP
% statevector(17): GAPDH
% statevector(18): GLC
% statevector(19): GLCex
% statevector(20): GLCfeed
% statevector(21): GOX
% statevector(22): Glk
% statevector(23): ICDH
% statevector(24): ICDHP
% statevector(25): ICIT
% statevector(26): Icl
% statevector(27): KDPG
% statevector(28): MAL
% statevector(29): MDH
% statevector(30): MS
% statevector(31): Mez
% statevector(32): OAA
% statevector(33): PDH
% statevector(34): PEP
% statevector(35): PYR
% statevector(36): Pck
% statevector(37): Pfk
% statevector(38): Ppc
% statevector(39): Pps
% statevector(40): Pyk
% statevector(41): R5P
% statevector(42): RU5P
% statevector(43): S7P
% statevector(44): SDH
% statevector(45): SUC
% statevector(46): X
% statevector(47): X5P
% statevector(48): aKG
% statevector(49): aKGDH
% statevector(50): cAMP
% statevector(51): sixPG
% statevector(52): sixPGL
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [1.11524e-10, 0.422054, 0.0237175, 0.00107739, 0.000406765, 0.0200874, 0.0261385, 0.00367648, 0.650133, 0.0992524, ...
		0.0372591, 0.0651811, 0.000460524, 0.0152583, 0.913391, 0.173332, 0.0652304, 0.012, 22.2222, 0, ...
		3.77198e-10, 0.00548198, 0.0154742, 0.0128716, 0.00736962, 0.0367581, 1.38744, 0.417978, 0.00976235, 0.0187144, ...
		0.0133654, 0.00254143, 0.00183994, 1.16312, 0.144336, 0.00360495, 0.00967106, 0.00323499, 0.0173145, 0.0107065, ...
		0.156479, 0.334105, 0.141549, 0.193564, 0.116135, 0.011, 0.486039, 0.0010271, 0.00819793, 0.246959, ...
		0.391757, 0.0109159];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'ACEex', 'AcCoA', 'AcP', 'AceK', 'Acs', 'CS', 'E4P', 'EIIAP', 'F6P', 'FBP', ...
			'FUM', 'Fba', 'Fbp', 'Fum', 'G6P', 'GAP', 'GAPDH', 'GLC', 'GLCex', 'GLCfeed', ...
			'GOX', 'Glk', 'ICDH', 'ICDHP', 'ICIT', 'Icl', 'KDPG', 'MAL', 'MDH', 'MS', ...
			'Mez', 'OAA', 'PDH', 'PEP', 'PYR', 'Pck', 'Pfk', 'Ppc', 'Pps', 'Pyk', ...
			'R5P', 'RU5P', 'S7P', 'SDH', 'SUC', 'X', 'X5P', 'aKG', 'aKGDH', 'cAMP', ...
			'sixPG', 'sixPGL'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'ADP', 'AMP', 'ATP', 'CoA', 'Cratotal', 'Crptotal', 'D', 'EIIAtotal', 'EXTERNAL', 'Factor_aceB', ...
			'Factor_aceK', 'IclRtotal', 'K6PGDH_6PG', 'K6PGDH_ATP_inh', 'K6PGDH_NADP', 'K6PGDH_NADPH_inh', 'KAceK_3PG', 'KAceK_GOX', 'KAceK_ICDH', 'KAceK_ICDHP', ...
			'KAceK_ICIT', 'KAceK_OAA', 'KAceK_PEP', 'KAceK_PYR', 'KAceK_aKG', 'KAck_ACE_m', 'KAck_ADP_m', 'KAck_ATP_m', 'KAck_AcP_m', 'KAck_eq', ...
			'KAcs_ACE', 'KCS_AcCoA', 'KCS_OAA', 'KCS_OAA_AcCoA', 'KCS_aKG', 'KCraFBP', 'KCrpcAMP', 'KCya_EIIAP', 'KEda_GAP_m', 'KEda_KDPG_m', ...
			'KEda_PYR_m', 'KEda_eq', 'KEdd_6PG_m', 'KEdd_KDPG_m', 'KEdd_eq', 'KFba_DHAP', 'KFba_FBP', 'KFba_GAP', 'KFba_GAP_inh', 'KFba_eq', ...
			'KFbp_FBP', 'KFbp_PEP', 'KFum_FUM_m', 'KFum_eq', 'KG6PDH_G6P', 'KG6PDH_NADP', 'KG6PDH_NADPH_g6pinh', 'KG6PDH_NADPH_nadpinh', 'KGAPDH_GAP', 'KGAPDH_NAD', ...
			'KGAPDH_NADH', 'KGAPDH_PGP', 'KGAPDH_eq', 'KGlk_ATP_m', 'KGlk_G6P_i', 'KGlk_GLC_m', 'KICDH_ICIT', 'KICDH_PEP', 'KIcl_3PG', 'KIcl_ICIT', ...
			'KIcl_PEP', 'KIcl_aKG', 'KMDH_MAL_I', 'KMDH_MAL_m', 'KMDH_NADH_I', 'KMDH_NADH_m', 'KMDH_NAD_I', 'KMDH_NAD_II', 'KMDH_NAD_m', 'KMDH_OAA_I', ...
			'KMDH_OAA_II', 'KMDH_OAA_m', 'KMDH_eq', 'KMS_AcCoA', 'KMS_GOX', 'KMS_GOX_AcCoA', 'KMez_AcCoA', 'KMez_MAL', 'KMez_cAMP', 'KNonPTS_I', ...
			'KNonPTS_S', 'KPDH_AcCoA_m', 'KPDH_CoA_m', 'KPDH_NADH_m', 'KPDH_NAD_m', 'KPDH_PYR_m', 'KPDH_i', 'KPTS_EIIA', 'KPTS_GLC', 'KPck_ADP_i', ...
			'KPck_ATP_I', 'KPck_ATP_i', 'KPck_OAA', 'KPck_OAA_I', 'KPck_PEP', 'KPck_PEP_i', 'KPdhRPYR', 'KPfk_ADP_a', 'KPfk_ADP_b', 'KPfk_ADP_c', ...
			'KPfk_AMP_a', 'KPfk_AMP_b', 'KPfk_ATP_s', 'KPfk_F6P_s', 'KPfk_PEP', 'KPgi_F6P', 'KPgi_F6P_6pginh', 'KPgi_G6P', 'KPgi_G6P_6pginh', 'KPgi_eq', ...
			'KPgl_6PGL_m', 'KPgl_6PG_m', 'KPgl_eq', 'KPgl_h1', 'KPgl_h2', 'KPpc_FBP', 'KPpc_PEP', 'KPps_PEP', 'KPps_PYR', 'KPta_AcCoA_i', ...
			'KPta_AcP_i', 'KPta_AcP_m', 'KPta_CoA_i', 'KPta_Pi_i', 'KPta_Pi_m', 'KPta_eq', 'KPyk_ADP', 'KPyk_AMP', 'KPyk_ATP', 'KPyk_FBP', ...
			'KPyk_PEP', 'KR5PI_eq', 'KRu5P_eq', 'KSDH_SUC_m', 'KSDH_eq', 'KTal_eq', 'KTktA_eq', 'KTktB_eq', 'KaKGDH_CoA_m', 'KaKGDH_NADH_I', ...
			'KaKGDH_NAD_m', 'KaKGDH_SUC_I', 'KaKGDH_Z', 'KaKGDH_aKG_I', 'KaKGDH_aKG_m', 'KaceBAK_Cra', 'KaceBAK_Crp', 'KaceBAK_DNA', 'KaceBAK_GOX', 'KaceBAK_PYR', ...
			'KaceBAK_PYRprime', 'Kacs_Crp', 'KcAMPdegr_cAMP', 'KfbaA_Cra', 'KfbaA_Crp', 'Kfbp_Cra', 'KfumABC_Crp', 'KgapA_Cra', 'KgapA_Crp', 'Kglk_Cra', ...
			'KgltA_Crp', 'KicdA_Cra', 'Kmdh_Crp', 'KpckA_Cra', 'Kpdh_PdhR', 'KpfkA_Cra', 'KppsA_Cra', 'KpykF_Cra', 'KsdhCDAB_Crp', 'KsucAB_Crp', ...
			'LAceK', 'LFbp', 'LICDH', 'LIcl', 'LMez', 'LPfk', 'LPpc', 'LPps', 'LPyk', 'LaceBAK', ...
			'NAD', 'NADH', 'NADP', 'NADPH', 'POratio', 'POratio_prime', 'PdhRtotal', 'Pi', 'SS_Mez_ACE', 'SS_Mez_GLC', ...
			'SS_Ppc_ACE', 'SS_Ppc_GLC', 'VFba_blf', 'aceBAK_DNA', 'kATP', 'kAceKki_cat', 'kAceKph_cat', 'kAcs_cat', 'kBM_ACE_AcCoA', 'kBM_ACE_E4P', ...
			'kBM_ACE_F6P', 'kBM_ACE_FUM', 'kBM_ACE_G6P', 'kBM_ACE_GAP', 'kBM_ACE_OAA', 'kBM_ACE_PEP', 'kBM_ACE_PYR', 'kBM_ACE_R5P', 'kBM_ACE_SUC', 'kBM_ACE_aKG', ...
			'kBM_GLC_AcCoA', 'kBM_GLC_E4P', 'kBM_GLC_F6P', 'kBM_GLC_FUM', 'kBM_GLC_G6P', 'kBM_GLC_GAP', 'kBM_GLC_OAA', 'kBM_GLC_PEP', 'kBM_GLC_PYR', 'kBM_GLC_R5P', ...
			'kBM_GLC_SUC', 'kBM_GLC_aKG', 'kCS_cat', 'kFba_cat', 'kFbp_cat', 'kFum1_cat', 'kFum2_cat', 'kGAPDH_cat', 'kGlk_cat', 'kICDH_cat', ...
			'kIcl_cat', 'kMDH1_cat', 'kMDH2_cat', 'kMS_cat', 'kMez_cat', 'kPDH_cat', 'kPTS1', 'kPck_cat', 'kPfk_cat', 'kPpc_cat', ...
			'kPps_cat', 'kPyk_cat', 'kSDH1_cat', 'kSDH2_cat', 'kaKGDH_cat', 'kaceBAK_cat_IclR', 'kdegr', 'kexpr', 'kmPTS1', 'nAceK', ...
			'nCraFBP', 'nCrpcAMP', 'nFbp', 'nICDH', 'nIcl', 'nMez', 'nPdhRPYR', 'nPfk', 'nPpc', 'nPps', ...
			'nPyk', 'nacs', 'nfumABC', 'ngltA', 'nsdhCDAB', 'nsucAB', 'pH', 'pH_Eda_m', 'pH_Edd_m', 'pK_Eda', ...
			'pK_Edd', 'rho', 'v6PGDH_max', 'vAck_max', 'vCya_max', 'vEda_max', 'vEdd_max', 'vG6PDH_max', 'vNonPTS_max', 'vPTS4_max', ...
			'vPgi_max', 'vPgl_max', 'vPta_max', 'vR5PI_max', 'vRu5P_max', 'vTal_max', 'vTktA_max', 'vTktB_max', 'vaceBAK_Cra_bound', 'vaceBAK_Cra_unbound', ...
			'vaceBAK_Crp_bound', 'vaceBAK_Crp_unbound', 'vacs_Crp_bound', 'vacs_Crp_unbound', 'vcAMPdegr_max', 'vfbaA_Cra_bound', 'vfbaA_Cra_unbound', 'vfbaA_Crp_bound', 'vfbaA_Crp_unbound', 'vfbp_Cra_bound', ...
			'vfbp_Cra_unbound', 'vfumABC_Crp_bound', 'vfumABC_Crp_unbound', 'vgapA_Cra_bound', 'vgapA_Cra_unbound', 'vgapA_Crp_bound', 'vgapA_Crp_unbound', 'vglk_Cra_bound', 'vglk_Cra_unbound', 'vgltA_Crp_bound', ...
			'vgltA_Crp_unbound', 'vicdA_Cra_bound', 'vicdA_Cra_unbound', 'vmdh_Crp_bound', 'vmdh_Crp_unbound', 'vpckA_Cra_bound', 'vpckA_Cra_unbound', 'vpdh_PdhR_bound', 'vpdh_PdhR_unbound', 'vpfkA_Cra_bound', ...
			'vpfkA_Cra_unbound', 'vppsA_Cra_bound', 'vppsA_Cra_unbound', 'vpykF_Cra_bound', 'vpykF_Cra_unbound', 'vsdhCDAB_Crp_bound', 'vsdhCDAB_Crp_unbound', 'vsucAB_Crp_bound', 'vsucAB_Crp_unbound', 'default0'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [0.56, 0.28, 9.6, 1.4, 0.0003, 0.0115, 0, 0.0769, 0, 0.509128, ...
			0.0293109, 8.3e-05, 0.0999956, 3.03764, 0.00774325, 0.0180077, 1.0927, 0.415957, 0.731925, 8.19069, ...
			0.0658771, 0.0345864, 0.232436, 0.0311645, 0.723859, 6.33002, 0.296006, 0.0881512, 0.0742685, 222.087, ...
			0.0233726, 0.172622, 0.0127217, 0.0177911, 0.271462, 0.0322704, 0.368586, 0.00242536, 1.37717, 0.350113, ...
			11.6552, 0.500041, 0.347444, 0.648821, 1000.01, 0.0879915, 0.13302, 0.0879908, 0.771367, 0.33082, ...
			0.00255028, 0.225362, 0.0554577, 18.9482, 0.0699847, 0.00607658, 0.192067, 0.0251647, 0.104745, 0.449961, ...
			0.0200034, 0.0995525, 0.217878, 0.799718, 14.9942, 0.219961, 9.55731e-05, 0.209919, 0.492465, 0.0257752, ...
			0.025628, 0.202044, 1.32753, 1.33028, 0.0242784, 0.0168477, 0.309937, 0.330748, 0.099998, 0.25207, ...
			0.396105, 0.269873, 0.557468, 0.381655, 0.361087, 0.311119, 2.99104, 0.00293864, 2.40366, 0.009166, ...
			2.34579, 0.0080012, 0.00698447, 0.248943, 0.399784, 1.20561, 21.7028, 0.00660253, 0.0127905, 0.0388343, ...
			0.0396311, 0.0415227, 0.669902, 0.532385, 0.0701171, 0.059769, 0.0856307, 238.935, 0.250106, 0.359964, ...
			6.05077, 0.0207231, 0.160031, 0.0231651, 1.93645, 0.199983, 0.19997, 2.45964, 0.200129, 0.717321, ...
			0.0229252, 9.96181, 42.7993, 0.00225241, 9.71558e-06, 0.126971, 0.0385597, 0.000604094, 0.0006405, 0.088354, ...
			0.0678373, 0.523402, 0.0244299, 1.61647, 0.699177, 0.0326788, 0.213774, 0.209251, 22.4085, 0.190009, ...
			0.370704, 0.484213, 1.62946, 0.220018, 8.67636, 1.17897, 1.2001, 10, 0.00391074, 0.0180044, ...
			0.0700101, 0.533941, 1.41777, 1.78176, 0.99931, 0.00018881, 4.56257, 1.13147e-06, 0.00566376, 0.249242, ...
			0.00858996, 0.00846341, 0.0662446, 0.00373129, 0.03915, 0.000122915, 0.122193, 0.00494581, 0.0296522, 2.91919e-08, ...
			0.0312059, 5.3399e-05, 0.255211, 0.000275018, 1.24734e-05, 2.41644e-08, 0.00062288, 0.000117808, 0.045294, 0.117486, ...
			7.59785e+07, 1.11817e+07, 126.46, 33725.8, 197779, 1.26952e+06, 3.8947e+06, 4.03223e-80, 1000.15, 1243.52, ...
			2.6, 0.083, 0.0021, 0.12, 3.30385, 1.50583, 6.66e-05, 10, 0.0334723, 0.0133738, ...
			0.000852337, 0.00323703, 2.00107, 5.15e-07, 1.30699e-05, 9.19676e+15, 8.85104e+12, 116474, 164.317, 284.114, ...
			321.291, 143.731, 118.486, 734.492, 7011.8, 202.218, 12876.3, 240.315, 150.808, 324.664, ...
			2656.7, 603.434, 966.423, 1091.82, 52.0836, 375.616, 19606.2, 708.44, 707.651, 247.339, ...
			2467.94, 4673.85, 792390, 1.60675e+06, 2.41475e+06, 681983, 681983, 8.79789e+07, 1.98397e+06, 211215, ...
			1.09342e+06, 328277, 328277, 13791.5, 1.65343e+06, 6.59379e+07, 42555.9, 2.60105e+06, 1.49637e+10, 1.0471e+07, ...
			1070.78, 17530.3, 22311.9, 22311.9, 2.04051e+08, 3.36933, 0.265324, 4.85939, 12994.3, 2, ...
			2, 1, 4, 2, 4, 1.92774, 1, 4, 3, 2, ...
			4, 2.31, 0.74, 1.07, 0.74, 0.74, 7.5, 7.49978, 7.15931, 13.1433, ...
			3.42317, 564, 193585, 533045, 13.2427, 582.202, 944.737, 53780.5, 4679.47, 10119.9, ...
			3.62077e+06, 45257.5, 5915.4, 42445.8, 23336.1, 42057.1, 7656.02, 285535, 0.0596534, 7.53474e-05, ...
			1.59563e-05, 0.00262305, 0.000592461, 0, 9.91179, 0, 0.0173187, 0.0150909, 0, 0.000701946, ...
			0, 0.0540414, 0, 0, 0.0156125, 0.0242728, 0, 0.00139121, 0.184677, 0.0505062, ...
			0, 0.0159016, 0.00401129, 0.158286, 0, 0.0110092, 0, 9.45249e-06, 0.00156527, 0.00272437, ...
			0.0754715, 0.113449, 0, 8.98315e-05, 0.0038128, 0.357936, 0, 0.0282708, 0, 1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'CraFBP', 'B', 'A', 'Cra', 'CrpcAMP', 'Crp', 'H', 'EIIA', 'alphaGLC', 'alphaACE', ...
			'SS_Ppc', 'SS_Mez', 'OP_NADH', 'OP_FADH2', 'PdhRPYR', 'PdhR', 'QEda_pH', 'QEdd_pH', 'vMDH1_max', 'vFum2_max', ...
			'vSDH2_max', 'vSDH1_max', 'vATP', 'mu', 'vMDH2_max', 'vFum1_max'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'Cratotal*power(FBP,nCraFBP)/(power(FBP,nCraFBP)+power(KCraFBP,nCraFBP))', '1+ADP/KPfk_ADP_a+AMP/KPfk_AMP_a', '1+PEP/KPfk_PEP+ADP/KPfk_ADP_b+AMP/KPfk_AMP_b', 'Cratotal-CraFBP', 'Crptotal*power(cAMP,nCrpcAMP)/(power(cAMP,nCrpcAMP)+power(KCrpcAMP,nCrpcAMP))', 'Crptotal-CrpcAMP', 'power(10,-pH)*1000', 'EIIAtotal-EIIAP', 'GLCex/(GLCex+KPTS_GLC)', 'ACEex/(ACEex+KAcs_ACE)*(1-alphaGLC)', ...
			'alphaGLC*SS_Ppc_GLC+alphaACE*SS_Ppc_ACE', 'alphaGLC*SS_Mez_GLC+alphaACE*SS_Mez_ACE', '(vE_GAPDH+vE_PDH+vE_aKGDH+vE_MDH)*POratio', 'vE_SDH*POratio_prime', 'PdhRtotal*power(PYR,nPdhRPYR)/(power(PYR,nPdhRPYR)+power(KPdhRPYR,nPdhRPYR))', 'PdhRtotal-PdhRPYR', '1+2*power(10,pH_Eda_m-pK_Eda)/(1+power(10,pH-pK_Eda)+power(10,2*pH_Eda_m-pH-pK_Eda))', '1+2*power(10,pH_Edd_m-pK_Edd)/(1+power(10,pH-pK_Edd)+power(10,2*pH_Edd_m-pH-pK_Edd))', 'MDH*kMDH1_cat', 'Fum*kFum2_cat', ...
			'SDH*kSDH2_cat', 'SDH*kSDH1_cat', 'OP_NADH+OP_FADH2-vE_Glk-vE_Pfk+vE_GAPDH+vE_Pyk-vE_Pps+vE_Ack-vE_Acs+vE_aKGDH-vE_Pck-vE_AceKki-vE_Cya', 'kATP*vATP', 'MDH*kMDH2_cat', 'Fum*kFum1_cat'};
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
	if length(parameterValuesNew) ~= 340,
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
ACEex = statevector(1);
AcCoA = statevector(2);
AcP = statevector(3);
AceK = statevector(4);
Acs = statevector(5);
CS = statevector(6);
E4P = statevector(7);
EIIAP = statevector(8);
F6P = statevector(9);
FBP = statevector(10);
FUM = statevector(11);
Fba = statevector(12);
Fbp = statevector(13);
Fum = statevector(14);
G6P = statevector(15);
GAP = statevector(16);
GAPDH = statevector(17);
GLC = statevector(18);
GLCex = statevector(19);
GLCfeed = statevector(20);
GOX = statevector(21);
Glk = statevector(22);
ICDH = statevector(23);
ICDHP = statevector(24);
ICIT = statevector(25);
Icl = statevector(26);
KDPG = statevector(27);
MAL = statevector(28);
MDH = statevector(29);
MS = statevector(30);
Mez = statevector(31);
OAA = statevector(32);
PDH = statevector(33);
PEP = statevector(34);
PYR = statevector(35);
Pck = statevector(36);
Pfk = statevector(37);
Ppc = statevector(38);
Pps = statevector(39);
Pyk = statevector(40);
R5P = statevector(41);
RU5P = statevector(42);
S7P = statevector(43);
SDH = statevector(44);
SUC = statevector(45);
X = statevector(46);
X5P = statevector(47);
aKG = statevector(48);
aKGDH = statevector(49);
cAMP = statevector(50);
sixPG = statevector(51);
sixPGL = statevector(52);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	ADP = 0.56;
	AMP = 0.28;
	ATP = 9.6;
	CoA = 1.4;
	Cratotal = 0.0003;
	Crptotal = 0.0115;
	D = 0;
	EIIAtotal = 0.0769;
	EXTERNAL = 0;
	Factor_aceB = 0.509128;
	Factor_aceK = 0.0293109;
	IclRtotal = 8.3e-05;
	K6PGDH_6PG = 0.0999956;
	K6PGDH_ATP_inh = 3.03764;
	K6PGDH_NADP = 0.00774325;
	K6PGDH_NADPH_inh = 0.0180077;
	KAceK_3PG = 1.0927;
	KAceK_GOX = 0.415957;
	KAceK_ICDH = 0.731925;
	KAceK_ICDHP = 8.19069;
	KAceK_ICIT = 0.0658771;
	KAceK_OAA = 0.0345864;
	KAceK_PEP = 0.232436;
	KAceK_PYR = 0.0311645;
	KAceK_aKG = 0.723859;
	KAck_ACE_m = 6.33002;
	KAck_ADP_m = 0.296006;
	KAck_ATP_m = 0.0881512;
	KAck_AcP_m = 0.0742685;
	KAck_eq = 222.087;
	KAcs_ACE = 0.0233726;
	KCS_AcCoA = 0.172622;
	KCS_OAA = 0.0127217;
	KCS_OAA_AcCoA = 0.0177911;
	KCS_aKG = 0.271462;
	KCraFBP = 0.0322704;
	KCrpcAMP = 0.368586;
	KCya_EIIAP = 0.00242536;
	KEda_GAP_m = 1.37717;
	KEda_KDPG_m = 0.350113;
	KEda_PYR_m = 11.6552;
	KEda_eq = 0.500041;
	KEdd_6PG_m = 0.347444;
	KEdd_KDPG_m = 0.648821;
	KEdd_eq = 1000.01;
	KFba_DHAP = 0.0879915;
	KFba_FBP = 0.13302;
	KFba_GAP = 0.0879908;
	KFba_GAP_inh = 0.771367;
	KFba_eq = 0.33082;
	KFbp_FBP = 0.00255028;
	KFbp_PEP = 0.225362;
	KFum_FUM_m = 0.0554577;
	KFum_eq = 18.9482;
	KG6PDH_G6P = 0.0699847;
	KG6PDH_NADP = 0.00607658;
	KG6PDH_NADPH_g6pinh = 0.192067;
	KG6PDH_NADPH_nadpinh = 0.0251647;
	KGAPDH_GAP = 0.104745;
	KGAPDH_NAD = 0.449961;
	KGAPDH_NADH = 0.0200034;
	KGAPDH_PGP = 0.0995525;
	KGAPDH_eq = 0.217878;
	KGlk_ATP_m = 0.799718;
	KGlk_G6P_i = 14.9942;
	KGlk_GLC_m = 0.219961;
	KICDH_ICIT = 9.55731e-05;
	KICDH_PEP = 0.209919;
	KIcl_3PG = 0.492465;
	KIcl_ICIT = 0.0257752;
	KIcl_PEP = 0.025628;
	KIcl_aKG = 0.202044;
	KMDH_MAL_I = 1.32753;
	KMDH_MAL_m = 1.33028;
	KMDH_NADH_I = 0.0242784;
	KMDH_NADH_m = 0.0168477;
	KMDH_NAD_I = 0.309937;
	KMDH_NAD_II = 0.330748;
	KMDH_NAD_m = 0.099998;
	KMDH_OAA_I = 0.25207;
	KMDH_OAA_II = 0.396105;
	KMDH_OAA_m = 0.269873;
	KMDH_eq = 0.557468;
	KMS_AcCoA = 0.381655;
	KMS_GOX = 0.361087;
	KMS_GOX_AcCoA = 0.311119;
	KMez_AcCoA = 2.99104;
	KMez_MAL = 0.00293864;
	KMez_cAMP = 2.40366;
	KNonPTS_I = 0.009166;
	KNonPTS_S = 2.34579;
	KPDH_AcCoA_m = 0.0080012;
	KPDH_CoA_m = 0.00698447;
	KPDH_NADH_m = 0.248943;
	KPDH_NAD_m = 0.399784;
	KPDH_PYR_m = 1.20561;
	KPDH_i = 21.7028;
	KPTS_EIIA = 0.00660253;
	KPTS_GLC = 0.0127905;
	KPck_ADP_i = 0.0388343;
	KPck_ATP_I = 0.0396311;
	KPck_ATP_i = 0.0415227;
	KPck_OAA = 0.669902;
	KPck_OAA_I = 0.532385;
	KPck_PEP = 0.0701171;
	KPck_PEP_i = 0.059769;
	KPdhRPYR = 0.0856307;
	KPfk_ADP_a = 238.935;
	KPfk_ADP_b = 0.250106;
	KPfk_ADP_c = 0.359964;
	KPfk_AMP_a = 6.05077;
	KPfk_AMP_b = 0.0207231;
	KPfk_ATP_s = 0.160031;
	KPfk_F6P_s = 0.0231651;
	KPfk_PEP = 1.93645;
	KPgi_F6P = 0.199983;
	KPgi_F6P_6pginh = 0.19997;
	KPgi_G6P = 2.45964;
	KPgi_G6P_6pginh = 0.200129;
	KPgi_eq = 0.717321;
	KPgl_6PGL_m = 0.0229252;
	KPgl_6PG_m = 9.96181;
	KPgl_eq = 42.7993;
	KPgl_h1 = 0.00225241;
	KPgl_h2 = 9.71558e-06;
	KPpc_FBP = 0.126971;
	KPpc_PEP = 0.0385597;
	KPps_PEP = 0.000604094;
	KPps_PYR = 0.0006405;
	KPta_AcCoA_i = 0.088354;
	KPta_AcP_i = 0.0678373;
	KPta_AcP_m = 0.523402;
	KPta_CoA_i = 0.0244299;
	KPta_Pi_i = 1.61647;
	KPta_Pi_m = 0.699177;
	KPta_eq = 0.0326788;
	KPyk_ADP = 0.213774;
	KPyk_AMP = 0.209251;
	KPyk_ATP = 22.4085;
	KPyk_FBP = 0.190009;
	KPyk_PEP = 0.370704;
	KR5PI_eq = 0.484213;
	KRu5P_eq = 1.62946;
	KSDH_SUC_m = 0.220018;
	KSDH_eq = 8.67636;
	KTal_eq = 1.17897;
	KTktA_eq = 1.2001;
	KTktB_eq = 10;
	KaKGDH_CoA_m = 0.00391074;
	KaKGDH_NADH_I = 0.0180044;
	KaKGDH_NAD_m = 0.0700101;
	KaKGDH_SUC_I = 0.533941;
	KaKGDH_Z = 1.41777;
	KaKGDH_aKG_I = 1.78176;
	KaKGDH_aKG_m = 0.99931;
	KaceBAK_Cra = 0.00018881;
	KaceBAK_Crp = 4.56257;
	KaceBAK_DNA = 1.13147e-06;
	KaceBAK_GOX = 0.00566376;
	KaceBAK_PYR = 0.249242;
	KaceBAK_PYRprime = 0.00858996;
	Kacs_Crp = 0.00846341;
	KcAMPdegr_cAMP = 0.0662446;
	KfbaA_Cra = 0.00373129;
	KfbaA_Crp = 0.03915;
	Kfbp_Cra = 0.000122915;
	KfumABC_Crp = 0.122193;
	KgapA_Cra = 0.00494581;
	KgapA_Crp = 0.0296522;
	Kglk_Cra = 2.91919e-08;
	KgltA_Crp = 0.0312059;
	KicdA_Cra = 5.3399e-05;
	Kmdh_Crp = 0.255211;
	KpckA_Cra = 0.000275018;
	Kpdh_PdhR = 1.24734e-05;
	KpfkA_Cra = 2.41644e-08;
	KppsA_Cra = 0.00062288;
	KpykF_Cra = 0.000117808;
	KsdhCDAB_Crp = 0.045294;
	KsucAB_Crp = 0.117486;
	LAceK = 7.59785e+07;
	LFbp = 1.11817e+07;
	LICDH = 126.46;
	LIcl = 33725.8;
	LMez = 197779;
	LPfk = 1.26952e+06;
	LPpc = 3.8947e+06;
	LPps = 4.03223e-80;
	LPyk = 1000.15;
	LaceBAK = 1243.52;
	NAD = 2.6;
	NADH = 0.083;
	NADP = 0.0021;
	NADPH = 0.12;
	POratio = 3.30385;
	POratio_prime = 1.50583;
	PdhRtotal = 6.66e-05;
	Pi = 10;
	SS_Mez_ACE = 0.0334723;
	SS_Mez_GLC = 0.0133738;
	SS_Ppc_ACE = 0.000852337;
	SS_Ppc_GLC = 0.00323703;
	VFba_blf = 2.00107;
	aceBAK_DNA = 5.15e-07;
	kATP = 1.30699e-05;
	kAceKki_cat = 9.19676e+15;
	kAceKph_cat = 8.85104e+12;
	kAcs_cat = 116474;
	kBM_ACE_AcCoA = 164.317;
	kBM_ACE_E4P = 284.114;
	kBM_ACE_F6P = 321.291;
	kBM_ACE_FUM = 143.731;
	kBM_ACE_G6P = 118.486;
	kBM_ACE_GAP = 734.492;
	kBM_ACE_OAA = 7011.8;
	kBM_ACE_PEP = 202.218;
	kBM_ACE_PYR = 12876.3;
	kBM_ACE_R5P = 240.315;
	kBM_ACE_SUC = 150.808;
	kBM_ACE_aKG = 324.664;
	kBM_GLC_AcCoA = 2656.7;
	kBM_GLC_E4P = 603.434;
	kBM_GLC_F6P = 966.423;
	kBM_GLC_FUM = 1091.82;
	kBM_GLC_G6P = 52.0836;
	kBM_GLC_GAP = 375.616;
	kBM_GLC_OAA = 19606.2;
	kBM_GLC_PEP = 708.44;
	kBM_GLC_PYR = 707.651;
	kBM_GLC_R5P = 247.339;
	kBM_GLC_SUC = 2467.94;
	kBM_GLC_aKG = 4673.85;
	kCS_cat = 792390;
	kFba_cat = 1.60675e+06;
	kFbp_cat = 2.41475e+06;
	kFum1_cat = 681983;
	kFum2_cat = 681983;
	kGAPDH_cat = 8.79789e+07;
	kGlk_cat = 1.98397e+06;
	kICDH_cat = 211215;
	kIcl_cat = 1.09342e+06;
	kMDH1_cat = 328277;
	kMDH2_cat = 328277;
	kMS_cat = 13791.5;
	kMez_cat = 1.65343e+06;
	kPDH_cat = 6.59379e+07;
	kPTS1 = 42555.9;
	kPck_cat = 2.60105e+06;
	kPfk_cat = 1.49637e+10;
	kPpc_cat = 1.0471e+07;
	kPps_cat = 1070.78;
	kPyk_cat = 17530.3;
	kSDH1_cat = 22311.9;
	kSDH2_cat = 22311.9;
	kaKGDH_cat = 2.04051e+08;
	kaceBAK_cat_IclR = 3.36933;
	kdegr = 0.265324;
	kexpr = 4.85939;
	kmPTS1 = 12994.3;
	nAceK = 2;
	nCraFBP = 2;
	nCrpcAMP = 1;
	nFbp = 4;
	nICDH = 2;
	nIcl = 4;
	nMez = 1.92774;
	nPdhRPYR = 1;
	nPfk = 4;
	nPpc = 3;
	nPps = 2;
	nPyk = 4;
	nacs = 2.31;
	nfumABC = 0.74;
	ngltA = 1.07;
	nsdhCDAB = 0.74;
	nsucAB = 0.74;
	pH = 7.5;
	pH_Eda_m = 7.49978;
	pH_Edd_m = 7.15931;
	pK_Eda = 13.1433;
	pK_Edd = 3.42317;
	rho = 564;
	v6PGDH_max = 193585;
	vAck_max = 533045;
	vCya_max = 13.2427;
	vEda_max = 582.202;
	vEdd_max = 944.737;
	vG6PDH_max = 53780.5;
	vNonPTS_max = 4679.47;
	vPTS4_max = 10119.9;
	vPgi_max = 3.62077e+06;
	vPgl_max = 45257.5;
	vPta_max = 5915.4;
	vR5PI_max = 42445.8;
	vRu5P_max = 23336.1;
	vTal_max = 42057.1;
	vTktA_max = 7656.02;
	vTktB_max = 285535;
	vaceBAK_Cra_bound = 0.0596534;
	vaceBAK_Cra_unbound = 7.53474e-05;
	vaceBAK_Crp_bound = 1.59563e-05;
	vaceBAK_Crp_unbound = 0.00262305;
	vacs_Crp_bound = 0.000592461;
	vacs_Crp_unbound = 0;
	vcAMPdegr_max = 9.91179;
	vfbaA_Cra_bound = 0;
	vfbaA_Cra_unbound = 0.0173187;
	vfbaA_Crp_bound = 0.0150909;
	vfbaA_Crp_unbound = 0;
	vfbp_Cra_bound = 0.000701946;
	vfbp_Cra_unbound = 0;
	vfumABC_Crp_bound = 0.0540414;
	vfumABC_Crp_unbound = 0;
	vgapA_Cra_bound = 0;
	vgapA_Cra_unbound = 0.0156125;
	vgapA_Crp_bound = 0.0242728;
	vgapA_Crp_unbound = 0;
	vglk_Cra_bound = 0.00139121;
	vglk_Cra_unbound = 0.184677;
	vgltA_Crp_bound = 0.0505062;
	vgltA_Crp_unbound = 0;
	vicdA_Cra_bound = 0.0159016;
	vicdA_Cra_unbound = 0.00401129;
	vmdh_Crp_bound = 0.158286;
	vmdh_Crp_unbound = 0;
	vpckA_Cra_bound = 0.0110092;
	vpckA_Cra_unbound = 0;
	vpdh_PdhR_bound = 9.45249e-06;
	vpdh_PdhR_unbound = 0.00156527;
	vpfkA_Cra_bound = 0.00272437;
	vpfkA_Cra_unbound = 0.0754715;
	vppsA_Cra_bound = 0.113449;
	vppsA_Cra_unbound = 0;
	vpykF_Cra_bound = 8.98315e-05;
	vpykF_Cra_unbound = 0.0038128;
	vsdhCDAB_Crp_bound = 0.357936;
	vsdhCDAB_Crp_unbound = 0;
	vsucAB_Crp_bound = 0.0282708;
	vsucAB_Crp_unbound = 0;
	default0 = 1;
else
	ADP = parameterValuesNew(1);
	AMP = parameterValuesNew(2);
	ATP = parameterValuesNew(3);
	CoA = parameterValuesNew(4);
	Cratotal = parameterValuesNew(5);
	Crptotal = parameterValuesNew(6);
	D = parameterValuesNew(7);
	EIIAtotal = parameterValuesNew(8);
	EXTERNAL = parameterValuesNew(9);
	Factor_aceB = parameterValuesNew(10);
	Factor_aceK = parameterValuesNew(11);
	IclRtotal = parameterValuesNew(12);
	K6PGDH_6PG = parameterValuesNew(13);
	K6PGDH_ATP_inh = parameterValuesNew(14);
	K6PGDH_NADP = parameterValuesNew(15);
	K6PGDH_NADPH_inh = parameterValuesNew(16);
	KAceK_3PG = parameterValuesNew(17);
	KAceK_GOX = parameterValuesNew(18);
	KAceK_ICDH = parameterValuesNew(19);
	KAceK_ICDHP = parameterValuesNew(20);
	KAceK_ICIT = parameterValuesNew(21);
	KAceK_OAA = parameterValuesNew(22);
	KAceK_PEP = parameterValuesNew(23);
	KAceK_PYR = parameterValuesNew(24);
	KAceK_aKG = parameterValuesNew(25);
	KAck_ACE_m = parameterValuesNew(26);
	KAck_ADP_m = parameterValuesNew(27);
	KAck_ATP_m = parameterValuesNew(28);
	KAck_AcP_m = parameterValuesNew(29);
	KAck_eq = parameterValuesNew(30);
	KAcs_ACE = parameterValuesNew(31);
	KCS_AcCoA = parameterValuesNew(32);
	KCS_OAA = parameterValuesNew(33);
	KCS_OAA_AcCoA = parameterValuesNew(34);
	KCS_aKG = parameterValuesNew(35);
	KCraFBP = parameterValuesNew(36);
	KCrpcAMP = parameterValuesNew(37);
	KCya_EIIAP = parameterValuesNew(38);
	KEda_GAP_m = parameterValuesNew(39);
	KEda_KDPG_m = parameterValuesNew(40);
	KEda_PYR_m = parameterValuesNew(41);
	KEda_eq = parameterValuesNew(42);
	KEdd_6PG_m = parameterValuesNew(43);
	KEdd_KDPG_m = parameterValuesNew(44);
	KEdd_eq = parameterValuesNew(45);
	KFba_DHAP = parameterValuesNew(46);
	KFba_FBP = parameterValuesNew(47);
	KFba_GAP = parameterValuesNew(48);
	KFba_GAP_inh = parameterValuesNew(49);
	KFba_eq = parameterValuesNew(50);
	KFbp_FBP = parameterValuesNew(51);
	KFbp_PEP = parameterValuesNew(52);
	KFum_FUM_m = parameterValuesNew(53);
	KFum_eq = parameterValuesNew(54);
	KG6PDH_G6P = parameterValuesNew(55);
	KG6PDH_NADP = parameterValuesNew(56);
	KG6PDH_NADPH_g6pinh = parameterValuesNew(57);
	KG6PDH_NADPH_nadpinh = parameterValuesNew(58);
	KGAPDH_GAP = parameterValuesNew(59);
	KGAPDH_NAD = parameterValuesNew(60);
	KGAPDH_NADH = parameterValuesNew(61);
	KGAPDH_PGP = parameterValuesNew(62);
	KGAPDH_eq = parameterValuesNew(63);
	KGlk_ATP_m = parameterValuesNew(64);
	KGlk_G6P_i = parameterValuesNew(65);
	KGlk_GLC_m = parameterValuesNew(66);
	KICDH_ICIT = parameterValuesNew(67);
	KICDH_PEP = parameterValuesNew(68);
	KIcl_3PG = parameterValuesNew(69);
	KIcl_ICIT = parameterValuesNew(70);
	KIcl_PEP = parameterValuesNew(71);
	KIcl_aKG = parameterValuesNew(72);
	KMDH_MAL_I = parameterValuesNew(73);
	KMDH_MAL_m = parameterValuesNew(74);
	KMDH_NADH_I = parameterValuesNew(75);
	KMDH_NADH_m = parameterValuesNew(76);
	KMDH_NAD_I = parameterValuesNew(77);
	KMDH_NAD_II = parameterValuesNew(78);
	KMDH_NAD_m = parameterValuesNew(79);
	KMDH_OAA_I = parameterValuesNew(80);
	KMDH_OAA_II = parameterValuesNew(81);
	KMDH_OAA_m = parameterValuesNew(82);
	KMDH_eq = parameterValuesNew(83);
	KMS_AcCoA = parameterValuesNew(84);
	KMS_GOX = parameterValuesNew(85);
	KMS_GOX_AcCoA = parameterValuesNew(86);
	KMez_AcCoA = parameterValuesNew(87);
	KMez_MAL = parameterValuesNew(88);
	KMez_cAMP = parameterValuesNew(89);
	KNonPTS_I = parameterValuesNew(90);
	KNonPTS_S = parameterValuesNew(91);
	KPDH_AcCoA_m = parameterValuesNew(92);
	KPDH_CoA_m = parameterValuesNew(93);
	KPDH_NADH_m = parameterValuesNew(94);
	KPDH_NAD_m = parameterValuesNew(95);
	KPDH_PYR_m = parameterValuesNew(96);
	KPDH_i = parameterValuesNew(97);
	KPTS_EIIA = parameterValuesNew(98);
	KPTS_GLC = parameterValuesNew(99);
	KPck_ADP_i = parameterValuesNew(100);
	KPck_ATP_I = parameterValuesNew(101);
	KPck_ATP_i = parameterValuesNew(102);
	KPck_OAA = parameterValuesNew(103);
	KPck_OAA_I = parameterValuesNew(104);
	KPck_PEP = parameterValuesNew(105);
	KPck_PEP_i = parameterValuesNew(106);
	KPdhRPYR = parameterValuesNew(107);
	KPfk_ADP_a = parameterValuesNew(108);
	KPfk_ADP_b = parameterValuesNew(109);
	KPfk_ADP_c = parameterValuesNew(110);
	KPfk_AMP_a = parameterValuesNew(111);
	KPfk_AMP_b = parameterValuesNew(112);
	KPfk_ATP_s = parameterValuesNew(113);
	KPfk_F6P_s = parameterValuesNew(114);
	KPfk_PEP = parameterValuesNew(115);
	KPgi_F6P = parameterValuesNew(116);
	KPgi_F6P_6pginh = parameterValuesNew(117);
	KPgi_G6P = parameterValuesNew(118);
	KPgi_G6P_6pginh = parameterValuesNew(119);
	KPgi_eq = parameterValuesNew(120);
	KPgl_6PGL_m = parameterValuesNew(121);
	KPgl_6PG_m = parameterValuesNew(122);
	KPgl_eq = parameterValuesNew(123);
	KPgl_h1 = parameterValuesNew(124);
	KPgl_h2 = parameterValuesNew(125);
	KPpc_FBP = parameterValuesNew(126);
	KPpc_PEP = parameterValuesNew(127);
	KPps_PEP = parameterValuesNew(128);
	KPps_PYR = parameterValuesNew(129);
	KPta_AcCoA_i = parameterValuesNew(130);
	KPta_AcP_i = parameterValuesNew(131);
	KPta_AcP_m = parameterValuesNew(132);
	KPta_CoA_i = parameterValuesNew(133);
	KPta_Pi_i = parameterValuesNew(134);
	KPta_Pi_m = parameterValuesNew(135);
	KPta_eq = parameterValuesNew(136);
	KPyk_ADP = parameterValuesNew(137);
	KPyk_AMP = parameterValuesNew(138);
	KPyk_ATP = parameterValuesNew(139);
	KPyk_FBP = parameterValuesNew(140);
	KPyk_PEP = parameterValuesNew(141);
	KR5PI_eq = parameterValuesNew(142);
	KRu5P_eq = parameterValuesNew(143);
	KSDH_SUC_m = parameterValuesNew(144);
	KSDH_eq = parameterValuesNew(145);
	KTal_eq = parameterValuesNew(146);
	KTktA_eq = parameterValuesNew(147);
	KTktB_eq = parameterValuesNew(148);
	KaKGDH_CoA_m = parameterValuesNew(149);
	KaKGDH_NADH_I = parameterValuesNew(150);
	KaKGDH_NAD_m = parameterValuesNew(151);
	KaKGDH_SUC_I = parameterValuesNew(152);
	KaKGDH_Z = parameterValuesNew(153);
	KaKGDH_aKG_I = parameterValuesNew(154);
	KaKGDH_aKG_m = parameterValuesNew(155);
	KaceBAK_Cra = parameterValuesNew(156);
	KaceBAK_Crp = parameterValuesNew(157);
	KaceBAK_DNA = parameterValuesNew(158);
	KaceBAK_GOX = parameterValuesNew(159);
	KaceBAK_PYR = parameterValuesNew(160);
	KaceBAK_PYRprime = parameterValuesNew(161);
	Kacs_Crp = parameterValuesNew(162);
	KcAMPdegr_cAMP = parameterValuesNew(163);
	KfbaA_Cra = parameterValuesNew(164);
	KfbaA_Crp = parameterValuesNew(165);
	Kfbp_Cra = parameterValuesNew(166);
	KfumABC_Crp = parameterValuesNew(167);
	KgapA_Cra = parameterValuesNew(168);
	KgapA_Crp = parameterValuesNew(169);
	Kglk_Cra = parameterValuesNew(170);
	KgltA_Crp = parameterValuesNew(171);
	KicdA_Cra = parameterValuesNew(172);
	Kmdh_Crp = parameterValuesNew(173);
	KpckA_Cra = parameterValuesNew(174);
	Kpdh_PdhR = parameterValuesNew(175);
	KpfkA_Cra = parameterValuesNew(176);
	KppsA_Cra = parameterValuesNew(177);
	KpykF_Cra = parameterValuesNew(178);
	KsdhCDAB_Crp = parameterValuesNew(179);
	KsucAB_Crp = parameterValuesNew(180);
	LAceK = parameterValuesNew(181);
	LFbp = parameterValuesNew(182);
	LICDH = parameterValuesNew(183);
	LIcl = parameterValuesNew(184);
	LMez = parameterValuesNew(185);
	LPfk = parameterValuesNew(186);
	LPpc = parameterValuesNew(187);
	LPps = parameterValuesNew(188);
	LPyk = parameterValuesNew(189);
	LaceBAK = parameterValuesNew(190);
	NAD = parameterValuesNew(191);
	NADH = parameterValuesNew(192);
	NADP = parameterValuesNew(193);
	NADPH = parameterValuesNew(194);
	POratio = parameterValuesNew(195);
	POratio_prime = parameterValuesNew(196);
	PdhRtotal = parameterValuesNew(197);
	Pi = parameterValuesNew(198);
	SS_Mez_ACE = parameterValuesNew(199);
	SS_Mez_GLC = parameterValuesNew(200);
	SS_Ppc_ACE = parameterValuesNew(201);
	SS_Ppc_GLC = parameterValuesNew(202);
	VFba_blf = parameterValuesNew(203);
	aceBAK_DNA = parameterValuesNew(204);
	kATP = parameterValuesNew(205);
	kAceKki_cat = parameterValuesNew(206);
	kAceKph_cat = parameterValuesNew(207);
	kAcs_cat = parameterValuesNew(208);
	kBM_ACE_AcCoA = parameterValuesNew(209);
	kBM_ACE_E4P = parameterValuesNew(210);
	kBM_ACE_F6P = parameterValuesNew(211);
	kBM_ACE_FUM = parameterValuesNew(212);
	kBM_ACE_G6P = parameterValuesNew(213);
	kBM_ACE_GAP = parameterValuesNew(214);
	kBM_ACE_OAA = parameterValuesNew(215);
	kBM_ACE_PEP = parameterValuesNew(216);
	kBM_ACE_PYR = parameterValuesNew(217);
	kBM_ACE_R5P = parameterValuesNew(218);
	kBM_ACE_SUC = parameterValuesNew(219);
	kBM_ACE_aKG = parameterValuesNew(220);
	kBM_GLC_AcCoA = parameterValuesNew(221);
	kBM_GLC_E4P = parameterValuesNew(222);
	kBM_GLC_F6P = parameterValuesNew(223);
	kBM_GLC_FUM = parameterValuesNew(224);
	kBM_GLC_G6P = parameterValuesNew(225);
	kBM_GLC_GAP = parameterValuesNew(226);
	kBM_GLC_OAA = parameterValuesNew(227);
	kBM_GLC_PEP = parameterValuesNew(228);
	kBM_GLC_PYR = parameterValuesNew(229);
	kBM_GLC_R5P = parameterValuesNew(230);
	kBM_GLC_SUC = parameterValuesNew(231);
	kBM_GLC_aKG = parameterValuesNew(232);
	kCS_cat = parameterValuesNew(233);
	kFba_cat = parameterValuesNew(234);
	kFbp_cat = parameterValuesNew(235);
	kFum1_cat = parameterValuesNew(236);
	kFum2_cat = parameterValuesNew(237);
	kGAPDH_cat = parameterValuesNew(238);
	kGlk_cat = parameterValuesNew(239);
	kICDH_cat = parameterValuesNew(240);
	kIcl_cat = parameterValuesNew(241);
	kMDH1_cat = parameterValuesNew(242);
	kMDH2_cat = parameterValuesNew(243);
	kMS_cat = parameterValuesNew(244);
	kMez_cat = parameterValuesNew(245);
	kPDH_cat = parameterValuesNew(246);
	kPTS1 = parameterValuesNew(247);
	kPck_cat = parameterValuesNew(248);
	kPfk_cat = parameterValuesNew(249);
	kPpc_cat = parameterValuesNew(250);
	kPps_cat = parameterValuesNew(251);
	kPyk_cat = parameterValuesNew(252);
	kSDH1_cat = parameterValuesNew(253);
	kSDH2_cat = parameterValuesNew(254);
	kaKGDH_cat = parameterValuesNew(255);
	kaceBAK_cat_IclR = parameterValuesNew(256);
	kdegr = parameterValuesNew(257);
	kexpr = parameterValuesNew(258);
	kmPTS1 = parameterValuesNew(259);
	nAceK = parameterValuesNew(260);
	nCraFBP = parameterValuesNew(261);
	nCrpcAMP = parameterValuesNew(262);
	nFbp = parameterValuesNew(263);
	nICDH = parameterValuesNew(264);
	nIcl = parameterValuesNew(265);
	nMez = parameterValuesNew(266);
	nPdhRPYR = parameterValuesNew(267);
	nPfk = parameterValuesNew(268);
	nPpc = parameterValuesNew(269);
	nPps = parameterValuesNew(270);
	nPyk = parameterValuesNew(271);
	nacs = parameterValuesNew(272);
	nfumABC = parameterValuesNew(273);
	ngltA = parameterValuesNew(274);
	nsdhCDAB = parameterValuesNew(275);
	nsucAB = parameterValuesNew(276);
	pH = parameterValuesNew(277);
	pH_Eda_m = parameterValuesNew(278);
	pH_Edd_m = parameterValuesNew(279);
	pK_Eda = parameterValuesNew(280);
	pK_Edd = parameterValuesNew(281);
	rho = parameterValuesNew(282);
	v6PGDH_max = parameterValuesNew(283);
	vAck_max = parameterValuesNew(284);
	vCya_max = parameterValuesNew(285);
	vEda_max = parameterValuesNew(286);
	vEdd_max = parameterValuesNew(287);
	vG6PDH_max = parameterValuesNew(288);
	vNonPTS_max = parameterValuesNew(289);
	vPTS4_max = parameterValuesNew(290);
	vPgi_max = parameterValuesNew(291);
	vPgl_max = parameterValuesNew(292);
	vPta_max = parameterValuesNew(293);
	vR5PI_max = parameterValuesNew(294);
	vRu5P_max = parameterValuesNew(295);
	vTal_max = parameterValuesNew(296);
	vTktA_max = parameterValuesNew(297);
	vTktB_max = parameterValuesNew(298);
	vaceBAK_Cra_bound = parameterValuesNew(299);
	vaceBAK_Cra_unbound = parameterValuesNew(300);
	vaceBAK_Crp_bound = parameterValuesNew(301);
	vaceBAK_Crp_unbound = parameterValuesNew(302);
	vacs_Crp_bound = parameterValuesNew(303);
	vacs_Crp_unbound = parameterValuesNew(304);
	vcAMPdegr_max = parameterValuesNew(305);
	vfbaA_Cra_bound = parameterValuesNew(306);
	vfbaA_Cra_unbound = parameterValuesNew(307);
	vfbaA_Crp_bound = parameterValuesNew(308);
	vfbaA_Crp_unbound = parameterValuesNew(309);
	vfbp_Cra_bound = parameterValuesNew(310);
	vfbp_Cra_unbound = parameterValuesNew(311);
	vfumABC_Crp_bound = parameterValuesNew(312);
	vfumABC_Crp_unbound = parameterValuesNew(313);
	vgapA_Cra_bound = parameterValuesNew(314);
	vgapA_Cra_unbound = parameterValuesNew(315);
	vgapA_Crp_bound = parameterValuesNew(316);
	vgapA_Crp_unbound = parameterValuesNew(317);
	vglk_Cra_bound = parameterValuesNew(318);
	vglk_Cra_unbound = parameterValuesNew(319);
	vgltA_Crp_bound = parameterValuesNew(320);
	vgltA_Crp_unbound = parameterValuesNew(321);
	vicdA_Cra_bound = parameterValuesNew(322);
	vicdA_Cra_unbound = parameterValuesNew(323);
	vmdh_Crp_bound = parameterValuesNew(324);
	vmdh_Crp_unbound = parameterValuesNew(325);
	vpckA_Cra_bound = parameterValuesNew(326);
	vpckA_Cra_unbound = parameterValuesNew(327);
	vpdh_PdhR_bound = parameterValuesNew(328);
	vpdh_PdhR_unbound = parameterValuesNew(329);
	vpfkA_Cra_bound = parameterValuesNew(330);
	vpfkA_Cra_unbound = parameterValuesNew(331);
	vppsA_Cra_bound = parameterValuesNew(332);
	vppsA_Cra_unbound = parameterValuesNew(333);
	vpykF_Cra_bound = parameterValuesNew(334);
	vpykF_Cra_unbound = parameterValuesNew(335);
	vsdhCDAB_Crp_bound = parameterValuesNew(336);
	vsdhCDAB_Crp_unbound = parameterValuesNew(337);
	vsucAB_Crp_bound = parameterValuesNew(338);
	vsucAB_Crp_unbound = parameterValuesNew(339);
	default0 = parameterValuesNew(340);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CraFBP = Cratotal*power(FBP,nCraFBP)/(power(FBP,nCraFBP)+power(KCraFBP,nCraFBP));
B = 1+ADP/KPfk_ADP_a+AMP/KPfk_AMP_a;
A = 1+PEP/KPfk_PEP+ADP/KPfk_ADP_b+AMP/KPfk_AMP_b;
Cra = Cratotal-CraFBP;
CrpcAMP = Crptotal*power(cAMP,nCrpcAMP)/(power(cAMP,nCrpcAMP)+power(KCrpcAMP,nCrpcAMP));
Crp = Crptotal-CrpcAMP;
H = power(10,-pH)*1000;
EIIA = EIIAtotal-EIIAP;
alphaGLC = GLCex/(GLCex+KPTS_GLC);
alphaACE = ACEex/(ACEex+KAcs_ACE)*(1-alphaGLC);
SS_Ppc = alphaGLC*SS_Ppc_GLC+alphaACE*SS_Ppc_ACE;
SS_Mez = alphaGLC*SS_Mez_GLC+alphaACE*SS_Mez_ACE;
OP_NADH = (vE_GAPDH+vE_PDH+vE_aKGDH+vE_MDH)*POratio;
OP_FADH2 = vE_SDH*POratio_prime;
PdhRPYR = PdhRtotal*power(PYR,nPdhRPYR)/(power(PYR,nPdhRPYR)+power(KPdhRPYR,nPdhRPYR));
PdhR = PdhRtotal-PdhRPYR;
QEda_pH = 1+2*power(10,pH_Eda_m-pK_Eda)/(1+power(10,pH-pK_Eda)+power(10,2*pH_Eda_m-pH-pK_Eda));
QEdd_pH = 1+2*power(10,pH_Edd_m-pK_Edd)/(1+power(10,pH-pK_Edd)+power(10,2*pH_Edd_m-pH-pK_Edd));
vMDH1_max = MDH*kMDH1_cat;
vFum2_max = Fum*kFum2_cat;
vSDH2_max = SDH*kSDH2_cat;
vSDH1_max = SDH*kSDH1_cat;
vATP = OP_NADH+OP_FADH2-vE_Glk-vE_Pfk+vE_GAPDH+vE_Pyk-vE_Pps+vE_Ack-vE_Acs+vE_aKGDH-vE_Pck-vE_AceKki-vE_Cya;
mu = kATP*vATP;
vMDH2_max = MDH*kMDH2_cat;
vFum1_max = Fum*kFum1_cat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vBM_AcCoA = default0 * Function_for_vBM_AcCoA(AcCoA, alphaACE, alphaGLC, default0, kBM_ACE_AcCoA, kBM_GLC_AcCoA);
vBM_E4P = default0 * Function_for_vBM_E4P(E4P, alphaACE, alphaGLC, default0, kBM_ACE_E4P, kBM_GLC_E4P);
vBM_F6P = default0 * Function_for_vBM_F6P(F6P, alphaACE, alphaGLC, default0, kBM_ACE_F6P, kBM_GLC_F6P);
vBM_FUM = default0 * Function_for_vBM_FUM(FUM, alphaACE, alphaGLC, default0, kBM_ACE_FUM, kBM_GLC_FUM);
vBM_G6P = default0 * Function_for_vBM_G6P(G6P, alphaACE, alphaGLC, default0, kBM_ACE_G6P, kBM_GLC_G6P);
vBM_GAP = default0 * Function_for_vBM_GAP(GAP, alphaACE, alphaGLC, default0, kBM_ACE_GAP, kBM_GLC_GAP);
vBM_OAA = default0 * Function_for_vBM_OAA(OAA, alphaACE, alphaGLC, default0, kBM_ACE_OAA, kBM_GLC_OAA);
vBM_PEP = default0 * Function_for_vBM_PEP(PEP, alphaACE, alphaGLC, default0, kBM_ACE_PEP, kBM_GLC_PEP);
vBM_PYR = default0 * Function_for_vBM_PYR(PYR, alphaACE, alphaGLC, default0, kBM_ACE_PYR, kBM_GLC_PYR);
vBM_R5P = default0 * Function_for_vBM_R5P(R5P, alphaACE, alphaGLC, default0, kBM_ACE_R5P, kBM_GLC_R5P);
vBM_SUC = default0 * Function_for_vBM_SUC(SUC, alphaACE, alphaGLC, default0, kBM_ACE_SUC, kBM_GLC_SUC);
vBM_aKG = default0 * Function_for_vBM_aKG(aKG, alphaACE, alphaGLC, default0, kBM_ACE_aKG, kBM_GLC_aKG);
vD_6PG = default0 * Function_for_vD_6PG(default0, mu, sixPG);
vD_6PGL = default0 * Function_for_vD_6PGL(default0, mu, sixPGL);
vD_ACEex = default0 * Function_for_vD_ACEex(ACEex, D, default0);
vD_AcCoA = default0 * Function_for_vD_AcCoA(AcCoA, default0, mu);
vD_AcP = default0 * Function_for_vD_AcP(AcP, default0, mu);
vD_AceK = default0 * Function_for_vD_AceK(AceK, default0, kdegr, mu);
vD_Acs = default0 * Function_for_vD_Acs(Acs, default0, kdegr, mu);
vD_CS = default0 * Function_for_vD_CS(CS, default0, kdegr, mu);
vD_E4P = default0 * Function_for_vD_E4P(E4P, default0, mu);
vD_F6P = default0 * Function_for_vD_F6P(F6P, default0, mu);
vD_FBP = default0 * Function_for_vD_FBP(FBP, default0, mu);
vD_FUM = default0 * Function_for_vD_FUM(FUM, default0, mu);
vD_Fba = default0 * Function_for_vD_Fba(Fba, default0, kdegr, mu);
vD_Fbp = default0 * Function_for_vD_Fbp(Fbp, default0, kdegr, mu);
vD_Fum = default0 * Function_for_vD_Fum(Fum, default0, kdegr, mu);
vD_G6P = default0 * Function_for_vD_G6P(G6P, default0, mu);
vD_GAP = default0 * Function_for_vD_GAP(GAP, default0, mu);
vD_GAPDH = default0 * Function_for_vD_GAPDH(GAPDH, default0, kdegr, mu);
vD_GLC = default0 * Function_for_vD_GLC(GLC, default0, mu);
vD_GLCex = default0 * Function_for_vD_GLCex(D, GLCex, default0);
vD_GLCfeed = default0 * Function_for_vD_GLCfeed(D, GLCfeed, default0);
vD_GOX = default0 * Function_for_vD_GOX(GOX, default0, mu);
vD_Glk = default0 * Function_for_vD_Glk(Glk, default0, kdegr, mu);
vD_ICDH = default0 * Function_for_vD_ICDH(ICDH, default0, kdegr, mu);
vD_ICDHP = default0 * Function_for_vD_ICDHP(ICDHP, default0, kdegr, mu);
vD_ICIT = default0 * Function_for_vD_ICIT(ICIT, default0, mu);
vD_Icl = default0 * Function_for_vD_Icl(Icl, default0, kdegr, mu);
vD_KDPG = default0 * Function_for_vD_KDPG(KDPG, default0, mu);
vD_MAL = default0 * Function_for_vD_MAL(MAL, default0, mu);
vD_MDH = default0 * Function_for_vD_MDH(MDH, default0, kdegr, mu);
vD_MS = default0 * Function_for_vD_MS(MS, default0, kdegr, mu);
vD_Mez = default0 * Function_for_vD_Mez(Mez, default0, kdegr, mu);
vD_OAA = default0 * Function_for_vD_OAA(OAA, default0, mu);
vD_PDH = default0 * Function_for_vD_PDH(PDH, default0, kdegr, mu);
vD_PEP = default0 * Function_for_vD_PEP(PEP, default0, mu);
vD_PYR = default0 * Function_for_vD_PYR(PYR, default0, mu);
vD_Pck = default0 * Function_for_vD_Pck(Pck, default0, kdegr, mu);
vD_Pfk = default0 * Function_for_vD_Pfk(Pfk, default0, kdegr, mu);
vD_Ppc = default0 * Function_for_vD_Ppc(Ppc, default0, kdegr, mu);
vD_Pps = default0 * Function_for_vD_Pps(Pps, default0, kdegr, mu);
vD_Pyk = default0 * Function_for_vD_Pyk(Pyk, default0, kdegr, mu);
vD_R5P = default0 * Function_for_vD_R5P(R5P, default0, mu);
vD_RU5P = default0 * Function_for_vD_RU5P(RU5P, default0, mu);
vD_S7P = default0 * Function_for_vD_S7P(S7P, default0, mu);
vD_SDH = default0 * Function_for_vD_SDH(SDH, default0, kdegr, mu);
vD_SUC = default0 * Function_for_vD_SUC(SUC, default0, mu);
vD_X = default0 * Function_for_vD_X(D, X, default0);
vD_X5P = default0 * Function_for_vD_X5P(X5P, default0, mu);
vD_aKG = default0 * Function_for_vD_aKG(aKG, default0, mu);
vD_aKGDH = default0 * Function_for_vD_aKGDH(aKGDH, default0, kdegr, mu);
vD_cAMP = default0 * Function_for_vD_cAMP(cAMP, default0, mu);
vE_6PGDH = default0 * Function_for_vE_6PGDH(ATP, K6PGDH_6PG, K6PGDH_ATP_inh, K6PGDH_NADP, K6PGDH_NADPH_inh, NADP, NADPH, default0, sixPG, v6PGDH_max);
vE_AceKki = default0 * Function_for_vE_AceKki(AceK, GAP, GOX, ICDH, ICIT, KAceK_3PG, KAceK_GOX, KAceK_ICDH, KAceK_ICIT, KAceK_OAA, KAceK_PEP, KAceK_PYR, KAceK_aKG, LAceK, OAA, PEP, PYR, aKG, default0, kAceKki_cat, nAceK);
vE_AceKph = default0 * Function_for_vE_AceKph(AceK, GAP, ICDHP, KAceK_3PG, KAceK_ICDHP, KAceK_OAA, KAceK_PEP, KAceK_PYR, KAceK_aKG, LAceK, OAA, PEP, PYR, aKG, default0, kAceKph_cat, nAceK);
vE_Ack = default0 * Function_for_vE_Ack(ACEex, ADP, ATP, AcP, KAck_ACE_m, KAck_ADP_m, KAck_ATP_m, KAck_AcP_m, KAck_eq, default0, vAck_max);
vE_Ack_medium = default0 * Function_for_vE_Ack_medium(ACEex, ADP, ATP, AcP, KAck_ACE_m, KAck_ADP_m, KAck_ATP_m, KAck_AcP_m, KAck_eq, X, default0, rho, vAck_max);
vE_Acs = default0 * Function_for_vE_Acs(ACEex, Acs, KAcs_ACE, default0, kAcs_cat);
vE_Acs_medium = default0 * Function_for_vE_Acs_medium(ACEex, Acs, KAcs_ACE, X, default0, kAcs_cat, rho);
vE_CS = default0 * Function_for_vE_CS(AcCoA, CS, KCS_AcCoA, KCS_OAA, KCS_OAA_AcCoA, KCS_aKG, OAA, aKG, default0, kCS_cat);
vE_Cya = default0 * Function_for_vE_Cya(EIIAP, KCya_EIIAP, default0, vCya_max);
vE_Eda = default0 * Function_for_vE_Eda(GAP, KDPG, KEda_GAP_m, KEda_KDPG_m, KEda_PYR_m, KEda_eq, PYR, QEda_pH, default0, vEda_max);
vE_Edd = default0 * Function_for_vE_Edd(KDPG, KEdd_6PG_m, KEdd_KDPG_m, KEdd_eq, QEdd_pH, default0, sixPG, vEdd_max);
vE_Fba = default0 * Function_for_vE_Fba(FBP, Fba, GAP, KFba_DHAP, KFba_FBP, KFba_GAP, KFba_GAP_inh, KFba_eq, VFba_blf, default0, kFba_cat);
vE_Fbp = default0 * Function_for_vE_Fbp(FBP, Fbp, KFbp_FBP, KFbp_PEP, LFbp, PEP, default0, kFbp_cat, nFbp);
vE_Fum = default0 * Function_for_vE_Fum(FUM, KFum_FUM_m, KFum_eq, MAL, default0, vFum1_max, vFum2_max);
vE_G6PDH = default0 * Function_for_vE_G6PDH(G6P, KG6PDH_G6P, KG6PDH_NADP, KG6PDH_NADPH_g6pinh, KG6PDH_NADPH_nadpinh, NADP, NADPH, default0, vG6PDH_max);
vE_GAPDH = default0 * Function_for_vE_GAPDH(GAP, GAPDH, KGAPDH_GAP, KGAPDH_NAD, KGAPDH_NADH, KGAPDH_PGP, KGAPDH_eq, NAD, NADH, PEP, default0, kGAPDH_cat);
vE_Glk = default0 * Function_for_vE_Glk(ATP, G6P, GLC, Glk, KGlk_ATP_m, KGlk_G6P_i, KGlk_GLC_m, default0, kGlk_cat);
vE_ICDH = default0 * Function_for_vE_ICDH(ICDH, ICIT, KICDH_ICIT, KICDH_PEP, LICDH, PEP, default0, kICDH_cat, nICDH);
vE_Icl = default0 * Function_for_vE_Icl(GAP, ICIT, Icl, KIcl_3PG, KIcl_ICIT, KIcl_PEP, KIcl_aKG, LIcl, PEP, aKG, default0, kIcl_cat, nIcl);
vE_MDH = default0 * Function_for_vE_MDH(KMDH_MAL_I, KMDH_MAL_m, KMDH_NADH_I, KMDH_NADH_m, KMDH_NAD_I, KMDH_NAD_II, KMDH_NAD_m, KMDH_OAA_I, KMDH_OAA_II, KMDH_OAA_m, KMDH_eq, MAL, NAD, NADH, OAA, default0, vMDH1_max, vMDH2_max);
vE_MS = default0 * Function_for_vE_MS(AcCoA, GOX, KMS_AcCoA, KMS_GOX, KMS_GOX_AcCoA, MS, default0, kMS_cat);
vE_Mez = default0 * Function_for_vE_Mez(AcCoA, KMez_AcCoA, KMez_MAL, KMez_cAMP, LMez, MAL, Mez, cAMP, default0, kMez_cat, nMez);
vE_PDH = default0 * Function_for_vE_PDH(AcCoA, CoA, KPDH_AcCoA_m, KPDH_CoA_m, KPDH_NADH_m, KPDH_NAD_m, KPDH_PYR_m, KPDH_i, NAD, NADH, PDH, PYR, default0, kPDH_cat);
vE_Pck = default0 * Function_for_vE_Pck(ADP, ATP, KPck_ADP_i, KPck_ATP_I, KPck_ATP_i, KPck_OAA, KPck_OAA_I, KPck_PEP, KPck_PEP_i, OAA, PEP, Pck, default0, kPck_cat);
vE_Pfk = default0 * Function_for_vE_Pfk(A, ADP, ATP, B, F6P, KPfk_ADP_c, KPfk_ATP_s, KPfk_F6P_s, LPfk, Pfk, default0, kPfk_cat, nPfk);
vE_Pgi = default0 * Function_for_vE_Pgi(F6P, G6P, KPgi_F6P, KPgi_F6P_6pginh, KPgi_G6P, KPgi_G6P_6pginh, KPgi_eq, default0, sixPG, vPgi_max);
vE_Pgl = default0 * Function_for_vE_Pgl(H, KPgl_6PGL_m, KPgl_6PG_m, KPgl_eq, KPgl_h1, KPgl_h2, default0, sixPG, sixPGL, vPgl_max);
vE_Ppc = default0 * Function_for_vE_Ppc(FBP, KPpc_FBP, KPpc_PEP, LPpc, PEP, Ppc, default0, kPpc_cat, nPpc);
vE_Pps = default0 * Function_for_vE_Pps(KPps_PEP, KPps_PYR, LPps, PEP, PYR, Pps, default0, kPps_cat, nPps);
vE_Pta = default0 * Function_for_vE_Pta(AcCoA, AcP, CoA, KPta_AcCoA_i, KPta_AcP_i, KPta_AcP_m, KPta_CoA_i, KPta_Pi_i, KPta_Pi_m, KPta_eq, Pi, default0, vPta_max);
vE_Pyk = default0 * Function_for_vE_Pyk(ADP, AMP, ATP, FBP, KPyk_ADP, KPyk_AMP, KPyk_ATP, KPyk_FBP, KPyk_PEP, LPyk, PEP, Pyk, default0, kPyk_cat, nPyk);
vE_R5PI = default0 * Function_for_vE_R5PI(KR5PI_eq, R5P, RU5P, default0, vR5PI_max);
vE_Ru5P = default0 * Function_for_vE_Ru5P(KRu5P_eq, RU5P, X5P, default0, vRu5P_max);
vE_SDH = default0 * Function_for_vE_SDH(FUM, KSDH_SUC_m, KSDH_eq, SUC, default0, vSDH1_max, vSDH2_max);
vE_Tal = default0 * Function_for_vE_Tal(E4P, F6P, GAP, KTal_eq, S7P, default0, vTal_max);
vE_TktA = default0 * Function_for_vE_TktA(GAP, KTktA_eq, R5P, S7P, X5P, default0, vTktA_max);
vE_TktB = default0 * Function_for_vE_TktB(E4P, F6P, GAP, KTktB_eq, X5P, default0, vTktB_max);
vE_aKGDH = default0 * Function_for_vE_aKGDH(CoA, KaKGDH_CoA_m, KaKGDH_NADH_I, KaKGDH_NAD_m, KaKGDH_SUC_I, KaKGDH_Z, KaKGDH_aKG_I, KaKGDH_aKG_m, NAD, NADH, SUC, aKG, aKGDH, default0, kaKGDH_cat);
vE_cAMPdegr = default0 * Function_for_vE_cAMPdegr(KcAMPdegr_cAMP, cAMP, default0, vcAMPdegr_max);
vG_aceA = default0 * Function_for_vG_aceA(Cra, CrpcAMP, GOX, IclRtotal, KaceBAK_Cra, KaceBAK_Crp, KaceBAK_DNA, KaceBAK_GOX, KaceBAK_PYR, KaceBAK_PYRprime, LaceBAK, PYR, aceBAK_DNA, default0, kaceBAK_cat_IclR, kexpr, mu, vaceBAK_Cra_bound, vaceBAK_Cra_unbound, vaceBAK_Crp_bound, vaceBAK_Crp_unbound);
vG_aceB = default0 * Function_for_vG_aceB(Cra, CrpcAMP, Factor_aceB, GOX, IclRtotal, KaceBAK_Cra, KaceBAK_Crp, KaceBAK_DNA, KaceBAK_GOX, KaceBAK_PYR, KaceBAK_PYRprime, LaceBAK, PYR, aceBAK_DNA, default0, kaceBAK_cat_IclR, kexpr, mu, vaceBAK_Cra_bound, vaceBAK_Cra_unbound, vaceBAK_Crp_bound, vaceBAK_Crp_unbound);
vG_aceK = default0 * Function_for_vG_aceK(Cra, CrpcAMP, Factor_aceK, GOX, IclRtotal, KaceBAK_Cra, KaceBAK_Crp, KaceBAK_DNA, KaceBAK_GOX, KaceBAK_PYR, KaceBAK_PYRprime, LaceBAK, PYR, aceBAK_DNA, default0, kaceBAK_cat_IclR, kexpr, mu, vaceBAK_Cra_bound, vaceBAK_Cra_unbound, vaceBAK_Crp_bound, vaceBAK_Crp_unbound);
vG_acs = default0 * Function_for_vG_acs(CrpcAMP, Kacs_Crp, default0, kexpr, mu, nacs, vacs_Crp_bound, vacs_Crp_unbound);
vG_fbaA = default0 * Function_for_vG_fbaA(Cra, CrpcAMP, KfbaA_Cra, KfbaA_Crp, default0, kexpr, mu, vfbaA_Cra_bound, vfbaA_Cra_unbound, vfbaA_Crp_bound, vfbaA_Crp_unbound);
vG_fbp = default0 * Function_for_vG_fbp(Cra, Kfbp_Cra, default0, kexpr, mu, vfbp_Cra_bound, vfbp_Cra_unbound);
vG_fumABC = default0 * Function_for_vG_fumABC(CrpcAMP, KfumABC_Crp, default0, kexpr, mu, nfumABC, vfumABC_Crp_bound, vfumABC_Crp_unbound);
vG_gapA = default0 * Function_for_vG_gapA(Cra, CrpcAMP, KgapA_Cra, KgapA_Crp, default0, kexpr, mu, vgapA_Cra_bound, vgapA_Cra_unbound, vgapA_Crp_bound, vgapA_Crp_unbound);
vG_glk = default0 * Function_for_vG_glk(Cra, Kglk_Cra, default0, kexpr, mu, vglk_Cra_bound, vglk_Cra_unbound);
vG_gltA = default0 * Function_for_vG_gltA(CrpcAMP, KgltA_Crp, default0, kexpr, mu, ngltA, vgltA_Crp_bound, vgltA_Crp_unbound);
vG_icdA = default0 * Function_for_vG_icdA(Cra, KicdA_Cra, default0, kexpr, mu, vicdA_Cra_bound, vicdA_Cra_unbound);
vG_maeB = default0 * Function_for_vG_maeB(SS_Mez, default0, kdegr, mu);
vG_mdh = default0 * Function_for_vG_mdh(CrpcAMP, Kmdh_Crp, default0, kexpr, mu, vmdh_Crp_bound, vmdh_Crp_unbound);
vG_pckA = default0 * Function_for_vG_pckA(Cra, KpckA_Cra, default0, kexpr, mu, vpckA_Cra_bound, vpckA_Cra_unbound);
vG_pdh = default0 * Function_for_vG_pdh(Kpdh_PdhR, PdhR, default0, kexpr, mu, vpdh_PdhR_bound, vpdh_PdhR_unbound);
vG_pfkA = default0 * Function_for_vG_pfkA(Cra, KpfkA_Cra, default0, kexpr, mu, vpfkA_Cra_bound, vpfkA_Cra_unbound);
vG_ppc = default0 * Function_for_vG_ppc(SS_Ppc, default0, kdegr, mu);
vG_ppsA = default0 * Function_for_vG_ppsA(Cra, KppsA_Cra, default0, kexpr, mu, vppsA_Cra_bound, vppsA_Cra_unbound);
vG_pykF = default0 * Function_for_vG_pykF(Cra, KpykF_Cra, default0, kexpr, mu, vpykF_Cra_bound, vpykF_Cra_unbound);
vG_sdhCDAB = default0 * Function_for_vG_sdhCDAB(CrpcAMP, KsdhCDAB_Crp, default0, kexpr, mu, nsdhCDAB, vsdhCDAB_Crp_bound, vsdhCDAB_Crp_unbound);
vG_sucAB = default0 * Function_for_vG_sucAB(CrpcAMP, KsucAB_Crp, default0, kexpr, mu, nsucAB, vsucAB_Crp_bound, vsucAB_Crp_unbound);
vNonPTS = default0 * Function_for_vNonPTS(EIIA, GLCex, KNonPTS_I, KNonPTS_S, default0, vNonPTS_max);
vNonPTS_medium = default0 * Function_for_vNonPTS_medium(EIIA, GLCex, KNonPTS_I, KNonPTS_S, X, default0, rho, vNonPTS_max);
vPTS1 = default0 * Function_for_vPTS1(EIIA, EIIAP, PEP, PYR, default0, kPTS1, kmPTS1);
vPTS4 = default0 * Function_for_vPTS4(EIIAP, GLCex, KPTS_EIIA, KPTS_GLC, default0, vPTS4_max);
vPTS4_medium = default0 * Function_for_vPTS4_medium(EIIAP, GLCex, KPTS_EIIA, KPTS_GLC, X, default0, rho, vPTS4_max);
vgrowth = default0 * Function_for_vgrowth(X, default0, mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACEex_dot = (-vD_ACEex+vE_Ack_medium-vE_Acs_medium)/default0;
AcCoA_dot = (-vBM_AcCoA-vD_AcCoA+vE_Acs-vE_CS-vE_MS+vE_PDH-vE_Pta)/default0;
AcP_dot = (-vD_AcP-vE_Ack+vE_Pta)/default0;
AceK_dot = (-vD_AceK+vG_aceK)/default0;
Acs_dot = (-vD_Acs+vG_acs)/default0;
CS_dot = (-vD_CS+vG_gltA)/default0;
E4P_dot = (-vBM_E4P-vD_E4P+vE_Tal-vE_TktB)/default0;
EIIAP_dot = (+vPTS1-vPTS4)/default0;
F6P_dot = (-vBM_F6P-vD_F6P+vE_Fbp-vE_Pfk+vE_Pgi+vE_Tal+vE_TktB)/default0;
FBP_dot = (-vD_FBP-vE_Fba-vE_Fbp+vE_Pfk)/default0;
FUM_dot = (-vBM_FUM-vD_FUM-vE_Fum+vE_SDH)/default0;
Fba_dot = (-vD_Fba+vG_fbaA)/default0;
Fbp_dot = (-vD_Fbp+vG_fbp)/default0;
Fum_dot = (-vD_Fum+vG_fumABC)/default0;
G6P_dot = (-vBM_G6P-vD_G6P-vE_G6PDH+vE_Glk-vE_Pgi+vPTS4)/default0;
GAP_dot = (-vBM_GAP-vD_GAP+vE_Eda+2*vE_Fba-vE_GAPDH-vE_Tal+vE_TktA+vE_TktB)/default0;
GAPDH_dot = (-vD_GAPDH+vG_gapA)/default0;
GLC_dot = (-vD_GLC-vE_Glk+vNonPTS)/default0;
GLCex_dot = (-vD_GLCex-vNonPTS_medium-vPTS4_medium)/default0;
GLCfeed_dot = (-vD_GLCfeed)/default0;
GOX_dot = (-vD_GOX+vE_Icl-vE_MS)/default0;
Glk_dot = (-vD_Glk+vG_glk)/default0;
ICDH_dot = (-vD_ICDH-vE_AceKki+vE_AceKph+vG_icdA)/default0;
ICDHP_dot = (-vD_ICDHP+vE_AceKki-vE_AceKph)/default0;
ICIT_dot = (-vD_ICIT+vE_CS-vE_ICDH-vE_Icl)/default0;
Icl_dot = (-vD_Icl+vG_aceA)/default0;
KDPG_dot = (-vD_KDPG-vE_Eda+vE_Edd)/default0;
MAL_dot = (-vD_MAL+vE_Fum-vE_MDH+vE_MS-vE_Mez)/default0;
MDH_dot = (-vD_MDH+vG_mdh)/default0;
MS_dot = (-vD_MS+vG_aceB)/default0;
Mez_dot = (-vD_Mez+vG_maeB)/default0;
OAA_dot = (-vBM_OAA-vD_OAA-vE_CS+vE_MDH-vE_Pck+vE_Ppc)/default0;
PDH_dot = (-vD_PDH+vG_pdh)/default0;
PEP_dot = (-vBM_PEP-vD_PEP+vE_GAPDH+vE_Pck-vE_Ppc+vE_Pps-vE_Pyk-vPTS1)/default0;
PYR_dot = (-vBM_PYR-vD_PYR+vE_Eda+vE_Mez-vE_PDH-vE_Pps+vE_Pyk+vPTS1)/default0;
Pck_dot = (-vD_Pck+vG_pckA)/default0;
Pfk_dot = (-vD_Pfk+vG_pfkA)/default0;
Ppc_dot = (-vD_Ppc+vG_ppc)/default0;
Pps_dot = (-vD_Pps+vG_ppsA)/default0;
Pyk_dot = (-vD_Pyk+vG_pykF)/default0;
R5P_dot = (-vBM_R5P-vD_R5P+vE_R5PI-vE_TktA)/default0;
RU5P_dot = (-vD_RU5P+vE_6PGDH-vE_R5PI-vE_Ru5P)/default0;
S7P_dot = (-vD_S7P-vE_Tal+vE_TktA)/default0;
SDH_dot = (-vD_SDH+vG_sdhCDAB)/default0;
SUC_dot = (-vBM_SUC-vD_SUC+vE_Icl-vE_SDH+vE_aKGDH)/default0;
X_dot = (-vD_X+vgrowth)/default0;
X5P_dot = (-vD_X5P+vE_Ru5P-vE_TktA-vE_TktB)/default0;
aKG_dot = (-vBM_aKG-vD_aKG+vE_ICDH-vE_aKGDH)/default0;
aKGDH_dot = (-vD_aKGDH+vG_sucAB)/default0;
cAMP_dot = (-vD_cAMP+vE_Cya-vE_cAMPdegr)/default0;
sixPG_dot = (-vD_6PG-vE_6PGDH-vE_Edd+vE_Pgl)/default0;
sixPGL_dot = (-vD_6PGL+vE_G6PDH-vE_Pgl)/default0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = ACEex_dot;
output(2) = AcCoA_dot;
output(3) = AcP_dot;
output(4) = AceK_dot;
output(5) = Acs_dot;
output(6) = CS_dot;
output(7) = E4P_dot;
output(8) = EIIAP_dot;
output(9) = F6P_dot;
output(10) = FBP_dot;
output(11) = FUM_dot;
output(12) = Fba_dot;
output(13) = Fbp_dot;
output(14) = Fum_dot;
output(15) = G6P_dot;
output(16) = GAP_dot;
output(17) = GAPDH_dot;
output(18) = GLC_dot;
output(19) = GLCex_dot;
output(20) = GLCfeed_dot;
output(21) = GOX_dot;
output(22) = Glk_dot;
output(23) = ICDH_dot;
output(24) = ICDHP_dot;
output(25) = ICIT_dot;
output(26) = Icl_dot;
output(27) = KDPG_dot;
output(28) = MAL_dot;
output(29) = MDH_dot;
output(30) = MS_dot;
output(31) = Mez_dot;
output(32) = OAA_dot;
output(33) = PDH_dot;
output(34) = PEP_dot;
output(35) = PYR_dot;
output(36) = Pck_dot;
output(37) = Pfk_dot;
output(38) = Ppc_dot;
output(39) = Pps_dot;
output(40) = Pyk_dot;
output(41) = R5P_dot;
output(42) = RU5P_dot;
output(43) = S7P_dot;
output(44) = SDH_dot;
output(45) = SUC_dot;
output(46) = X_dot;
output(47) = X5P_dot;
output(48) = aKG_dot;
output(49) = aKGDH_dot;
output(50) = cAMP_dot;
output(51) = sixPG_dot;
output(52) = sixPGL_dot;
% return a column vector 
output = output(:);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = Function_for_vD_CS(CS,default0,kdegr,mu)
global time
result = (mu+kdegr)*CS/default0;
return

function [result] = Function_for_vD_Fum(Fum,default0,kdegr,mu)
global time
result = (mu+kdegr)*Fum/default0;
return

function [result] = Function_for_vE_TktA(GAP,KTktA_eq,R5P,S7P,X5P,default0,vTktA_max)
global time
result = vTktA_max*(R5P*X5P-S7P*GAP/KTktA_eq)/default0;
return

function [result] = Function_for_vBM_G6P(G6P,alphaACE,alphaGLC,default0,kBM_ACE_G6P,kBM_GLC_G6P)
global time
result = (alphaGLC*kBM_GLC_G6P+alphaACE*kBM_ACE_G6P)*G6P/default0;
return

function [result] = Function_for_vD_F6P(F6P,default0,mu)
global time
result = mu*F6P/default0;
return

function [result] = Function_for_vD_E4P(E4P,default0,mu)
global time
result = mu*E4P/default0;
return

function [result] = Function_for_vE_CS(AcCoA,CS,KCS_AcCoA,KCS_OAA,KCS_OAA_AcCoA,KCS_aKG,OAA,aKG,default0,kCS_cat)
global time
result = CS*kCS_cat*OAA*AcCoA/((1+aKG/KCS_aKG)*KCS_OAA_AcCoA*KCS_AcCoA+KCS_AcCoA*OAA+(1+aKG/KCS_aKG)*KCS_OAA*AcCoA+OAA*AcCoA)/default0;
return

function [result] = Function_for_vE_Pta(AcCoA,AcP,CoA,KPta_AcCoA_i,KPta_AcP_i,KPta_AcP_m,KPta_CoA_i,KPta_Pi_i,KPta_Pi_m,KPta_eq,Pi,default0,vPta_max)
global time
result = vPta_max*(1/(KPta_AcCoA_i*KPta_Pi_m))*(AcCoA*Pi-AcP*CoA/KPta_eq)/(1+AcCoA/KPta_AcCoA_i+Pi/KPta_Pi_i+AcP/KPta_AcP_i+CoA/KPta_CoA_i+AcCoA*Pi/(KPta_AcCoA_i*KPta_Pi_m)+AcP*CoA/(KPta_AcP_m*KPta_CoA_i))/default0;
return

function [result] = Function_for_vG_maeB(SS_Mez,default0,kdegr,mu)
global time
result = (mu+kdegr)*SS_Mez/default0;
return

function [result] = Function_for_vG_gltA(CrpcAMP,KgltA_Crp,default0,kexpr,mu,ngltA,vgltA_Crp_bound,vgltA_Crp_unbound)
global time
result = mu*kexpr*((1-(CrpcAMP)^(ngltA)/((CrpcAMP)^(ngltA)+(KgltA_Crp)^(ngltA)))*vgltA_Crp_unbound+(CrpcAMP)^(ngltA)/((CrpcAMP)^(ngltA)+(KgltA_Crp)^(ngltA))*vgltA_Crp_bound)/default0;
return

function [result] = Function_for_vE_Fba(FBP,Fba,GAP,KFba_DHAP,KFba_FBP,KFba_GAP,KFba_GAP_inh,KFba_eq,VFba_blf,default0,kFba_cat)
global time
result = Fba*kFba_cat*(FBP-(GAP)^(2)/KFba_eq)/(KFba_FBP+FBP+KFba_GAP*GAP/(KFba_eq*VFba_blf)+KFba_DHAP*GAP/(KFba_eq*VFba_blf)+FBP*GAP/KFba_GAP_inh+(GAP)^(2)/(KFba_eq*VFba_blf))/default0;
return

function [result] = Function_for_vG_mdh(CrpcAMP,Kmdh_Crp,default0,kexpr,mu,vmdh_Crp_bound,vmdh_Crp_unbound)
global time
result = mu*kexpr*((1-CrpcAMP/(CrpcAMP+Kmdh_Crp))*vmdh_Crp_unbound+CrpcAMP/(CrpcAMP+Kmdh_Crp)*vmdh_Crp_bound)/default0;
return

function [result] = Function_for_vG_icdA(Cra,KicdA_Cra,default0,kexpr,mu,vicdA_Cra_bound,vicdA_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KicdA_Cra))*vicdA_Cra_unbound+Cra/(Cra+KicdA_Cra)*vicdA_Cra_bound)/default0;
return

function [result] = Function_for_vE_SDH(FUM,KSDH_SUC_m,KSDH_eq,SUC,default0,vSDH1_max,vSDH2_max)
global time
result = vSDH1_max*vSDH2_max*(SUC-FUM/KSDH_eq)/(KSDH_SUC_m*vSDH2_max+vSDH2_max*SUC+vSDH1_max*FUM/KSDH_eq)/default0;
return

function [result] = Function_for_vE_Pgl(H,KPgl_6PGL_m,KPgl_6PG_m,KPgl_eq,KPgl_h1,KPgl_h2,default0,sixPG,sixPGL,vPgl_max)
global time
result = vPgl_max*(sixPGL-sixPG/KPgl_eq)/((1+H/KPgl_h1+KPgl_h2/H)*(KPgl_6PGL_m+sixPGL+KPgl_6PGL_m/KPgl_6PG_m*sixPG))/default0;
return

function [result] = Function_for_vE_Tal(E4P,F6P,GAP,KTal_eq,S7P,default0,vTal_max)
global time
result = vTal_max*(GAP*S7P-E4P*F6P/KTal_eq)/default0;
return

function [result] = Function_for_vE_Pgi(F6P,G6P,KPgi_F6P,KPgi_F6P_6pginh,KPgi_G6P,KPgi_G6P_6pginh,KPgi_eq,default0,sixPG,vPgi_max)
global time
result = vPgi_max*(G6P-F6P/KPgi_eq)/(KPgi_G6P*(1+F6P/(KPgi_F6P*(1+sixPG/KPgi_F6P_6pginh))+sixPG/KPgi_G6P_6pginh)+G6P)/default0;
return

function [result] = Function_for_vD_FUM(FUM,default0,mu)
global time
result = mu*FUM/default0;
return

function [result] = Function_for_vD_G6P(G6P,default0,mu)
global time
result = mu*G6P/default0;
return

function [result] = Function_for_vD_Fbp(Fbp,default0,kdegr,mu)
global time
result = (mu+kdegr)*Fbp/default0;
return

function [result] = Function_for_vD_GLC(GLC,default0,mu)
global time
result = mu*GLC/default0;
return

function [result] = Function_for_vE_R5PI(KR5PI_eq,R5P,RU5P,default0,vR5PI_max)
global time
result = vR5PI_max*(RU5P-R5P/KR5PI_eq)/default0;
return

function [result] = Function_for_vG_aceK(Cra,CrpcAMP,Factor_aceK,GOX,IclRtotal,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR,KaceBAK_PYRprime,LaceBAK,PYR,aceBAK_DNA,default0,kaceBAK_cat_IclR,kexpr,mu,vaceBAK_Cra_bound,vaceBAK_Cra_unbound,vaceBAK_Crp_bound,vaceBAK_Crp_unbound)
global time
result = Factor_aceK*mu*kexpr*((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)/default0;
return

function [result] = Function_for_vG_glk(Cra,Kglk_Cra,default0,kexpr,mu,vglk_Cra_bound,vglk_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+Kglk_Cra))*vglk_Cra_unbound+Cra/(Cra+Kglk_Cra)*vglk_Cra_bound)/default0;
return

function [result] = Function_for_vD_GAPDH(GAPDH,default0,kdegr,mu)
global time
result = (mu+kdegr)*GAPDH/default0;
return

function [result] = Function_for_vD_GLCex(D,GLCex,default0)
global time
result = D*GLCex/default0;
return

function [result] = Function_for_vD_GAP(GAP,default0,mu)
global time
result = mu*GAP/default0;
return

function [result] = Function_for_vD_Fba(Fba,default0,kdegr,mu)
global time
result = (mu+kdegr)*Fba/default0;
return

function [result] = Function_for_vD_FBP(FBP,default0,mu)
global time
result = mu*FBP/default0;
return

function [result] = Function_for_vE_Mez(AcCoA,KMez_AcCoA,KMez_MAL,KMez_cAMP,LMez,MAL,Mez,cAMP,default0,kMez_cat,nMez)
global time
result = Mez*kMez_cat*MAL/KMez_MAL*(1+MAL/KMez_MAL)^(nMez-1)/((1+MAL/KMez_MAL)^(nMez)+LMez*(1+AcCoA/KMez_AcCoA+cAMP/KMez_cAMP)^(nMez))/default0;
return

function [result] = Function_for_vD_Ppc(Ppc,default0,kdegr,mu)
global time
result = (mu+kdegr)*Ppc/default0;
return

function [result] = Function_for_vD_aKGDH(aKGDH,default0,kdegr,mu)
global time
result = (mu+kdegr)*aKGDH/default0;
return

function [result] = Function_for_vE_Ppc(FBP,KPpc_FBP,KPpc_PEP,LPpc,PEP,Ppc,default0,kPpc_cat,nPpc)
global time
result = Ppc*kPpc_cat*PEP/KPpc_PEP*(1+PEP/KPpc_PEP)^(nPpc-1)/((1+PEP/KPpc_PEP)^(nPpc)+LPpc/(1+FBP/KPpc_FBP)^(nPpc))/default0;
return

function [result] = Function_for_vE_G6PDH(G6P,KG6PDH_G6P,KG6PDH_NADP,KG6PDH_NADPH_g6pinh,KG6PDH_NADPH_nadpinh,NADP,NADPH,default0,vG6PDH_max)
global time
result = vG6PDH_max*G6P*NADP/((G6P+KG6PDH_G6P)*(1+NADPH/KG6PDH_NADPH_g6pinh)*(KG6PDH_NADP*(1+NADPH/KG6PDH_NADPH_nadpinh)+NADP))/default0;
return

function [result] = Function_for_vE_Glk(ATP,G6P,GLC,Glk,KGlk_ATP_m,KGlk_G6P_i,KGlk_GLC_m,default0,kGlk_cat)
global time
result = Glk*kGlk_cat*(GLC/KGlk_GLC_m)*(ATP/(KGlk_ATP_m*(1+G6P/KGlk_G6P_i)))/(1+GLC/KGlk_GLC_m+ATP/(KGlk_ATP_m*(1+G6P/KGlk_G6P_i))+GLC*ATP/(KGlk_GLC_m*KGlk_ATP_m*(1+G6P/KGlk_G6P_i))+G6P/KGlk_G6P_i)/default0;
return

function [result] = Function_for_vE_MDH(KMDH_MAL_I,KMDH_MAL_m,KMDH_NADH_I,KMDH_NADH_m,KMDH_NAD_I,KMDH_NAD_II,KMDH_NAD_m,KMDH_OAA_I,KMDH_OAA_II,KMDH_OAA_m,KMDH_eq,MAL,NAD,NADH,OAA,default0,vMDH1_max,vMDH2_max)
global time
result = vMDH1_max*vMDH2_max*(NAD*MAL-NADH*OAA/KMDH_eq)/(KMDH_NAD_I*KMDH_MAL_m*vMDH2_max+KMDH_MAL_m*vMDH2_max*NAD+KMDH_NAD_m*vMDH2_max*MAL+vMDH2_max*NAD*MAL+KMDH_OAA_m*vMDH1_max*NADH/KMDH_eq+KMDH_NADH_m*vMDH1_max*OAA/KMDH_eq+vMDH1_max*NADH*OAA/KMDH_eq+vMDH1_max*KMDH_OAA_m*NAD*NADH/(KMDH_eq*KMDH_NAD_I)+vMDH2_max*KMDH_NAD_m*MAL*OAA/KMDH_OAA_I+vMDH2_max*NAD*MAL*NADH/KMDH_NADH_I+vMDH1_max*MAL*NADH*OAA/(KMDH_eq*KMDH_MAL_I)+vMDH2_max*NAD*MAL*OAA/KMDH_OAA_II+vMDH1_max*NAD*NADH*OAA/(KMDH_NAD_II*KMDH_eq)+KMDH_NAD_I*vMDH2_max*NAD*MAL*NADH*OAA/(KMDH_NAD_II*KMDH_OAA_m*KMDH_NADH_I))/default0;
return

function [result] = Function_for_vD_RU5P(RU5P,default0,mu)
global time
result = mu*RU5P/default0;
return

function [result] = Function_for_vE_GAPDH(GAP,GAPDH,KGAPDH_GAP,KGAPDH_NAD,KGAPDH_NADH,KGAPDH_PGP,KGAPDH_eq,NAD,NADH,PEP,default0,kGAPDH_cat)
global time
result = GAPDH*kGAPDH_cat*(GAP*NAD-PEP*NADH/KGAPDH_eq)/((KGAPDH_GAP*(1+PEP/KGAPDH_PGP)+GAP)*(KGAPDH_NAD*(1+NADH/KGAPDH_NADH)+NAD))/default0;
return

function [result] = Function_for_vD_OAA(OAA,default0,mu)
global time
result = mu*OAA/default0;
return

function [result] = Function_for_vE_Fbp(FBP,Fbp,KFbp_FBP,KFbp_PEP,LFbp,PEP,default0,kFbp_cat,nFbp)
global time
result = Fbp*kFbp_cat*FBP/KFbp_FBP*(1+FBP/KFbp_FBP)^(nFbp-1)/((1+FBP/KFbp_FBP)^(nFbp)+LFbp/(1+PEP/KFbp_PEP)^(nFbp))/default0;
return

function [result] = Function_for_vE_MS(AcCoA,GOX,KMS_AcCoA,KMS_GOX,KMS_GOX_AcCoA,MS,default0,kMS_cat)
global time
result = MS*kMS_cat*GOX*AcCoA/(KMS_GOX_AcCoA*KMS_AcCoA+KMS_AcCoA*GOX+KMS_GOX*AcCoA+GOX*AcCoA)/default0;
return

function [result] = Function_for_vD_R5P(R5P,default0,mu)
global time
result = mu*R5P/default0;
return

function [result] = Function_for_vD_Pfk(Pfk,default0,kdegr,mu)
global time
result = (mu+kdegr)*Pfk/default0;
return

function [result] = Function_for_vD_PYR(PYR,default0,mu)
global time
result = mu*PYR/default0;
return

function [result] = Function_for_vD_X(D,X,default0)
global time
result = D*X/default0;
return

function [result] = Function_for_vE_Icl(GAP,ICIT,Icl,KIcl_3PG,KIcl_ICIT,KIcl_PEP,KIcl_aKG,LIcl,PEP,aKG,default0,kIcl_cat,nIcl)
global time
result = Icl*kIcl_cat*ICIT/KIcl_ICIT*(1+ICIT/KIcl_ICIT)^(nIcl-1)/((1+ICIT/KIcl_ICIT)^(nIcl)+LIcl*(1+PEP/KIcl_PEP+GAP/KIcl_3PG+aKG/KIcl_aKG)^(nIcl))/default0;
return

function [result] = Function_for_vE_Pps(KPps_PEP,KPps_PYR,LPps,PEP,PYR,Pps,default0,kPps_cat,nPps)
global time
result = Pps*kPps_cat*PYR/KPps_PYR*(1+PYR/KPps_PYR)^(nPps-1)/((1+PYR/KPps_PYR)^(nPps)+LPps*(1+PEP/KPps_PEP)^(nPps))/default0;
return

function [result] = Function_for_vD_SDH(SDH,default0,kdegr,mu)
global time
result = (mu+kdegr)*SDH/default0;
return

function [result] = Function_for_vD_Pck(Pck,default0,kdegr,mu)
global time
result = (mu+kdegr)*Pck/default0;
return

function [result] = Function_for_vD_PEP(PEP,default0,mu)
global time
result = mu*PEP/default0;
return

function [result] = Function_for_vD_PDH(PDH,default0,kdegr,mu)
global time
result = (mu+kdegr)*PDH/default0;
return

function [result] = Function_for_vD_Pyk(Pyk,default0,kdegr,mu)
global time
result = (mu+kdegr)*Pyk/default0;
return

function [result] = Function_for_vE_Ack_medium(ACEex,ADP,ATP,AcP,KAck_ACE_m,KAck_ADP_m,KAck_ATP_m,KAck_AcP_m,KAck_eq,X,default0,rho,vAck_max)
global time
result = vAck_max*(1/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1+ADP/KAck_ADP_m+ATP/KAck_ATP_m))*X/rho/default0;
return

function [result] = Function_for_vE_Acs(ACEex,Acs,KAcs_ACE,default0,kAcs_cat)
global time
result = Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE)/default0;
return

function [result] = Function_for_vE_AceKki(AceK,GAP,GOX,ICDH,ICIT,KAceK_3PG,KAceK_GOX,KAceK_ICDH,KAceK_ICIT,KAceK_OAA,KAceK_PEP,KAceK_PYR,KAceK_aKG,LAceK,OAA,PEP,PYR,aKG,default0,kAceKki_cat,nAceK)
global time
result = AceK*kAceKki_cat*ICDH/KAceK_ICDH*(1+ICDH/KAceK_ICDH)^(nAceK-1)/((1+ICDH/KAceK_ICDH)^(nAceK)+LAceK*(1+ICIT/KAceK_ICIT+GOX/KAceK_GOX+OAA/KAceK_OAA+aKG/KAceK_aKG+PEP/KAceK_PEP+GAP/KAceK_3PG+PYR/KAceK_PYR)^(nAceK))/default0;
return

function [result] = Function_for_vD_SUC(SUC,default0,mu)
global time
result = mu*SUC/default0;
return

function [result] = Function_for_vD_Pps(Pps,default0,kdegr,mu)
global time
result = (mu+kdegr)*Pps/default0;
return

function [result] = Function_for_vD_S7P(S7P,default0,mu)
global time
result = mu*S7P/default0;
return

function [result] = Function_for_vD_aKG(aKG,default0,mu)
global time
result = mu*aKG/default0;
return

function [result] = Function_for_vE_Acs_medium(ACEex,Acs,KAcs_ACE,X,default0,kAcs_cat,rho)
global time
result = Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE)*X/rho/default0;
return

function [result] = Function_for_vE_6PGDH(ATP,K6PGDH_6PG,K6PGDH_ATP_inh,K6PGDH_NADP,K6PGDH_NADPH_inh,NADP,NADPH,default0,sixPG,v6PGDH_max)
global time
result = v6PGDH_max*sixPG*NADP/((sixPG+K6PGDH_6PG)*(NADP+K6PGDH_NADP*(1+NADPH/K6PGDH_NADPH_inh)*(1+ATP/K6PGDH_ATP_inh)))/default0;
return

function [result] = Function_for_vE_Ack(ACEex,ADP,ATP,AcP,KAck_ACE_m,KAck_ADP_m,KAck_ATP_m,KAck_AcP_m,KAck_eq,default0,vAck_max)
global time
result = vAck_max*(1/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1+ADP/KAck_ADP_m+ATP/KAck_ATP_m))/default0;
return

function [result] = Function_for_vE_Cya(EIIAP,KCya_EIIAP,default0,vCya_max)
global time
result = vCya_max*EIIAP/(EIIAP+KCya_EIIAP)/default0;
return

function [result] = Function_for_vG_aceA(Cra,CrpcAMP,GOX,IclRtotal,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR,KaceBAK_PYRprime,LaceBAK,PYR,aceBAK_DNA,default0,kaceBAK_cat_IclR,kexpr,mu,vaceBAK_Cra_bound,vaceBAK_Cra_unbound,vaceBAK_Crp_bound,vaceBAK_Crp_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)/default0;
return

function [result] = Function_for_vG_fbaA(Cra,CrpcAMP,KfbaA_Cra,KfbaA_Crp,default0,kexpr,mu,vfbaA_Cra_bound,vfbaA_Cra_unbound,vfbaA_Crp_bound,vfbaA_Crp_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KfbaA_Cra))*vfbaA_Cra_unbound+Cra/(Cra+KfbaA_Cra)*vfbaA_Cra_bound+(1-CrpcAMP/(CrpcAMP+KfbaA_Crp))*vfbaA_Crp_unbound+CrpcAMP/(CrpcAMP+KfbaA_Crp)*vfbaA_Crp_bound)/default0;
return

function [result] = Function_for_vG_acs(CrpcAMP,Kacs_Crp,default0,kexpr,mu,nacs,vacs_Crp_bound,vacs_Crp_unbound)
global time
result = mu*kexpr*((1-(CrpcAMP)^(nacs)/((CrpcAMP)^(nacs)+(Kacs_Crp)^(nacs)))*vacs_Crp_unbound+(CrpcAMP)^(nacs)/((CrpcAMP)^(nacs)+(Kacs_Crp)^(nacs))*vacs_Crp_bound)/default0;
return

function [result] = Function_for_vD_X5P(X5P,default0,mu)
global time
result = mu*X5P/default0;
return

function [result] = Function_for_vD_GLCfeed(D,GLCfeed,default0)
global time
result = D*GLCfeed/default0;
return

function [result] = Function_for_vE_Edd(KDPG,KEdd_6PG_m,KEdd_KDPG_m,KEdd_eq,QEdd_pH,default0,sixPG,vEdd_max)
global time
result = vEdd_max*QEdd_pH*(sixPG-KDPG/KEdd_eq)/(KEdd_6PG_m+sixPG+KEdd_6PG_m*KDPG/KEdd_KDPG_m)/default0;
return

function [result] = Function_for_vG_fbp(Cra,Kfbp_Cra,default0,kexpr,mu,vfbp_Cra_bound,vfbp_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+Kfbp_Cra))*vfbp_Cra_unbound+Cra/(Cra+Kfbp_Cra)*vfbp_Cra_bound)/default0;
return

function [result] = Function_for_vG_fumABC(CrpcAMP,KfumABC_Crp,default0,kexpr,mu,nfumABC,vfumABC_Crp_bound,vfumABC_Crp_unbound)
global time
result = mu*kexpr*((1-(CrpcAMP)^(nfumABC)/((CrpcAMP)^(nfumABC)+(KfumABC_Crp)^(nfumABC)))*vfumABC_Crp_unbound+(CrpcAMP)^(nfumABC)/((CrpcAMP)^(nfumABC)+(KfumABC_Crp)^(nfumABC))*vfumABC_Crp_bound)/default0;
return

function [result] = Function_for_vE_Eda(GAP,KDPG,KEda_GAP_m,KEda_KDPG_m,KEda_PYR_m,KEda_eq,PYR,QEda_pH,default0,vEda_max)
global time
result = vEda_max*QEda_pH*(KDPG-GAP*PYR/KEda_eq)/(KEda_KDPG_m+KDPG+KEda_KDPG_m*(PYR/KEda_PYR_m+GAP/KEda_GAP_m+PYR*GAP/(KEda_PYR_m*KEda_GAP_m)))/default0;
return

function [result] = Function_for_vG_gapA(Cra,CrpcAMP,KgapA_Cra,KgapA_Crp,default0,kexpr,mu,vgapA_Cra_bound,vgapA_Cra_unbound,vgapA_Crp_bound,vgapA_Crp_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KgapA_Cra))*vgapA_Cra_unbound+Cra/(Cra+KgapA_Cra)*vgapA_Cra_bound+(1-CrpcAMP/(CrpcAMP+KgapA_Crp))*vgapA_Crp_unbound+CrpcAMP/(CrpcAMP+KgapA_Crp)*vgapA_Crp_bound)/default0;
return

function [result] = Function_for_vE_AceKph(AceK,GAP,ICDHP,KAceK_3PG,KAceK_ICDHP,KAceK_OAA,KAceK_PEP,KAceK_PYR,KAceK_aKG,LAceK,OAA,PEP,PYR,aKG,default0,kAceKph_cat,nAceK)
global time
result = AceK*kAceKph_cat*ICDHP/KAceK_ICDHP*(1+ICDHP/KAceK_ICDHP)^(nAceK-1)/((1+ICDHP/KAceK_ICDHP)^(nAceK)+LAceK/(1+OAA/KAceK_OAA+aKG/KAceK_aKG+PEP/KAceK_PEP+GAP/KAceK_3PG+PYR/KAceK_PYR)^(nAceK))/default0;
return

function [result] = Function_for_vD_cAMP(cAMP,default0,mu)
global time
result = mu*cAMP/default0;
return

function [result] = Function_for_vG_aceB(Cra,CrpcAMP,Factor_aceB,GOX,IclRtotal,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR,KaceBAK_PYRprime,LaceBAK,PYR,aceBAK_DNA,default0,kaceBAK_cat_IclR,kexpr,mu,vaceBAK_Cra_bound,vaceBAK_Cra_unbound,vaceBAK_Crp_bound,vaceBAK_Crp_unbound)
global time
result = Factor_aceB*mu*kexpr*((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)/default0;
return

function [result] = Function_for_vPTS1(EIIA,EIIAP,PEP,PYR,default0,kPTS1,kmPTS1)
global time
result = (kPTS1*PEP*EIIA-kmPTS1*PYR*EIIAP)/default0;
return

function [result] = Function_for_vE_Pyk(ADP,AMP,ATP,FBP,KPyk_ADP,KPyk_AMP,KPyk_ATP,KPyk_FBP,KPyk_PEP,LPyk,PEP,Pyk,default0,kPyk_cat,nPyk)
global time
result = Pyk*kPyk_cat*PEP*(PEP/KPyk_PEP+1)^(nPyk-1)*ADP/(KPyk_PEP*(LPyk*power((1+ATP/KPyk_ATP)/(FBP/KPyk_FBP+AMP/KPyk_AMP+1),nPyk)+(PEP/KPyk_PEP+1)^(nPyk))*(ADP+KPyk_ADP))/default0;
return

function [result] = Function_for_vE_Ru5P(KRu5P_eq,RU5P,X5P,default0,vRu5P_max)
global time
result = vRu5P_max*(RU5P-X5P/KRu5P_eq)/default0;
return

function [result] = Function_for_vE_Pfk(A,ADP,ATP,B,F6P,KPfk_ADP_c,KPfk_ATP_s,KPfk_F6P_s,LPfk,Pfk,default0,kPfk_cat,nPfk)
global time
result = Pfk*kPfk_cat*ATP*F6P/((ATP+KPfk_ATP_s*(1+ADP/KPfk_ADP_c))*(F6P+KPfk_F6P_s*A/B)*(1+LPfk/(1+F6P*B/(KPfk_F6P_s*A))^(nPfk)))/default0;
return

function [result] = Function_for_vG_sdhCDAB(CrpcAMP,KsdhCDAB_Crp,default0,kexpr,mu,nsdhCDAB,vsdhCDAB_Crp_bound,vsdhCDAB_Crp_unbound)
global time
result = mu*kexpr*((1-(CrpcAMP)^(nsdhCDAB)/((CrpcAMP)^(nsdhCDAB)+(KsdhCDAB_Crp)^(nsdhCDAB)))*vsdhCDAB_Crp_unbound+(CrpcAMP)^(nsdhCDAB)/((CrpcAMP)^(nsdhCDAB)+(KsdhCDAB_Crp)^(nsdhCDAB))*vsdhCDAB_Crp_bound)/default0;
return

function [result] = Function_for_vBM_GAP(GAP,alphaACE,alphaGLC,default0,kBM_ACE_GAP,kBM_GLC_GAP)
global time
result = (alphaGLC*kBM_GLC_GAP+alphaACE*kBM_ACE_GAP)*GAP/default0;
return

function [result] = Function_for_vD_ICDHP(ICDHP,default0,kdegr,mu)
global time
result = (mu+kdegr)*ICDHP/default0;
return

function [result] = Function_for_vD_MDH(MDH,default0,kdegr,mu)
global time
result = (mu+kdegr)*MDH/default0;
return

function [result] = Function_for_vPTS4(EIIAP,GLCex,KPTS_EIIA,KPTS_GLC,default0,vPTS4_max)
global time
result = vPTS4_max*EIIAP*GLCex/((KPTS_EIIA+EIIAP)*(KPTS_GLC+GLCex))/default0;
return

function [result] = Function_for_vgrowth(X,default0,mu)
global time
result = mu*X/default0;
return

function [result] = Function_for_vPTS4_medium(EIIAP,GLCex,KPTS_EIIA,KPTS_GLC,X,default0,rho,vPTS4_max)
global time
result = vPTS4_max*EIIAP*GLCex/((KPTS_EIIA+EIIAP)*(KPTS_GLC+GLCex))*X/rho/default0;
return

function [result] = Function_for_vNonPTS_medium(EIIA,GLCex,KNonPTS_I,KNonPTS_S,X,default0,rho,vNonPTS_max)
global time
result = vNonPTS_max*GLCex/(KNonPTS_S+(1+EIIA/KNonPTS_I)*GLCex)*X/rho/default0;
return

function [result] = Function_for_vNonPTS(EIIA,GLCex,KNonPTS_I,KNonPTS_S,default0,vNonPTS_max)
global time
result = vNonPTS_max*GLCex/(KNonPTS_S+(1+EIIA/KNonPTS_I)*GLCex)/default0;
return

function [result] = Function_for_vD_GOX(GOX,default0,mu)
global time
result = mu*GOX/default0;
return

function [result] = Function_for_vE_cAMPdegr(KcAMPdegr_cAMP,cAMP,default0,vcAMPdegr_max)
global time
result = vcAMPdegr_max*cAMP/(cAMP+KcAMPdegr_cAMP)/default0;
return

function [result] = Function_for_vBM_E4P(E4P,alphaACE,alphaGLC,default0,kBM_ACE_E4P,kBM_GLC_E4P)
global time
result = (alphaGLC*kBM_GLC_E4P+alphaACE*kBM_ACE_E4P)*E4P/default0;
return

function [result] = Function_for_vE_TktB(E4P,F6P,GAP,KTktB_eq,X5P,default0,vTktB_max)
global time
result = vTktB_max*(X5P*E4P-F6P*GAP/KTktB_eq)/default0;
return

function [result] = Function_for_vG_sucAB(CrpcAMP,KsucAB_Crp,default0,kexpr,mu,nsucAB,vsucAB_Crp_bound,vsucAB_Crp_unbound)
global time
result = mu*kexpr*((1-(CrpcAMP)^(nsucAB)/((CrpcAMP)^(nsucAB)+(KsucAB_Crp)^(nsucAB)))*vsucAB_Crp_unbound+(CrpcAMP)^(nsucAB)/((CrpcAMP)^(nsucAB)+(KsucAB_Crp)^(nsucAB))*vsucAB_Crp_bound)/default0;
return

function [result] = Function_for_vE_Pck(ADP,ATP,KPck_ADP_i,KPck_ATP_I,KPck_ATP_i,KPck_OAA,KPck_OAA_I,KPck_PEP,KPck_PEP_i,OAA,PEP,Pck,default0,kPck_cat)
global time
result = Pck*kPck_cat*OAA*ATP/ADP/(KPck_OAA*ATP/ADP+OAA*ATP/ADP+KPck_ATP_i*KPck_OAA/KPck_ADP_i+KPck_ATP_i*KPck_OAA/(KPck_PEP*KPck_ADP_i)*PEP+KPck_ATP_i*KPck_OAA/(KPck_PEP_i*KPck_ATP_I)*ATP/ADP*PEP+KPck_ATP_i*KPck_OAA/(KPck_ADP_i*KPck_OAA_I)*OAA)/default0;
return

function [result] = Function_for_vG_pfkA(Cra,KpfkA_Cra,default0,kexpr,mu,vpfkA_Cra_bound,vpfkA_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KpfkA_Cra))*vpfkA_Cra_unbound+Cra/(Cra+KpfkA_Cra)*vpfkA_Cra_bound)/default0;
return

function [result] = Function_for_vG_pckA(Cra,KpckA_Cra,default0,kexpr,mu,vpckA_Cra_bound,vpckA_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KpckA_Cra))*vpckA_Cra_unbound+Cra/(Cra+KpckA_Cra)*vpckA_Cra_bound)/default0;
return

function [result] = Function_for_vE_aKGDH(CoA,KaKGDH_CoA_m,KaKGDH_NADH_I,KaKGDH_NAD_m,KaKGDH_SUC_I,KaKGDH_Z,KaKGDH_aKG_I,KaKGDH_aKG_m,NAD,NADH,SUC,aKG,aKGDH,default0,kaKGDH_cat)
global time
result = aKGDH*kaKGDH_cat*aKG*CoA*NAD/(KaKGDH_NAD_m*aKG*CoA+KaKGDH_CoA_m*aKG*NAD+KaKGDH_aKG_m*CoA*NAD+aKG*CoA*NAD+KaKGDH_aKG_m*KaKGDH_Z*SUC*NADH/KaKGDH_SUC_I+KaKGDH_NAD_m*aKG*CoA*NADH/KaKGDH_NADH_I+KaKGDH_CoA_m*aKG*NAD*SUC/KaKGDH_SUC_I+KaKGDH_aKG_m*KaKGDH_Z*aKG*SUC*NADH/(KaKGDH_aKG_I*KaKGDH_SUC_I))/default0;
return

function [result] = Function_for_vD_ICDH(ICDH,default0,kdegr,mu)
global time
result = (mu+kdegr)*ICDH/default0;
return

function [result] = Function_for_vD_MS(MS,default0,kdegr,mu)
global time
result = (mu+kdegr)*MS/default0;
return

function [result] = Function_for_vD_Icl(Icl,default0,kdegr,mu)
global time
result = (mu+kdegr)*Icl/default0;
return

function [result] = Function_for_vG_pdh(Kpdh_PdhR,PdhR,default0,kexpr,mu,vpdh_PdhR_bound,vpdh_PdhR_unbound)
global time
result = mu*kexpr*((1-PdhR/(PdhR+Kpdh_PdhR))*vpdh_PdhR_unbound+PdhR/(PdhR+Kpdh_PdhR)*vpdh_PdhR_bound)/default0;
return

function [result] = Function_for_vD_Mez(Mez,default0,kdegr,mu)
global time
result = (mu+kdegr)*Mez/default0;
return

function [result] = Function_for_vG_ppc(SS_Ppc,default0,kdegr,mu)
global time
result = (mu+kdegr)*SS_Ppc/default0;
return

function [result] = Function_for_vD_Glk(Glk,default0,kdegr,mu)
global time
result = (mu+kdegr)*Glk/default0;
return

function [result] = Function_for_vG_ppsA(Cra,KppsA_Cra,default0,kexpr,mu,vppsA_Cra_bound,vppsA_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KppsA_Cra))*vppsA_Cra_unbound+Cra/(Cra+KppsA_Cra)*vppsA_Cra_bound)/default0;
return

function [result] = Function_for_vBM_AcCoA(AcCoA,alphaACE,alphaGLC,default0,kBM_ACE_AcCoA,kBM_GLC_AcCoA)
global time
result = (alphaGLC*kBM_GLC_AcCoA+alphaACE*kBM_ACE_AcCoA)*AcCoA/default0;
return

function [result] = Function_for_vD_KDPG(KDPG,default0,mu)
global time
result = mu*KDPG/default0;
return

function [result] = Function_for_vG_pykF(Cra,KpykF_Cra,default0,kexpr,mu,vpykF_Cra_bound,vpykF_Cra_unbound)
global time
result = mu*kexpr*((1-Cra/(Cra+KpykF_Cra))*vpykF_Cra_unbound+Cra/(Cra+KpykF_Cra)*vpykF_Cra_bound)/default0;
return

function [result] = Function_for_vBM_PEP(PEP,alphaACE,alphaGLC,default0,kBM_ACE_PEP,kBM_GLC_PEP)
global time
result = (alphaGLC*kBM_GLC_PEP+alphaACE*kBM_ACE_PEP)*PEP/default0;
return

function [result] = Function_for_vD_MAL(MAL,default0,mu)
global time
result = mu*MAL/default0;
return

function [result] = Function_for_vD_6PG(default0,mu,sixPG)
global time
result = mu*sixPG/default0;
return

function [result] = Function_for_vD_AcCoA(AcCoA,default0,mu)
global time
result = mu*AcCoA/default0;
return

function [result] = Function_for_vD_AcP(AcP,default0,mu)
global time
result = mu*AcP/default0;
return

function [result] = Function_for_vD_Acs(Acs,default0,kdegr,mu)
global time
result = (mu+kdegr)*Acs/default0;
return

function [result] = Function_for_vD_6PGL(default0,mu,sixPGL)
global time
result = mu*sixPGL/default0;
return

function [result] = Function_for_vBM_aKG(aKG,alphaACE,alphaGLC,default0,kBM_ACE_aKG,kBM_GLC_aKG)
global time
result = (alphaGLC*kBM_GLC_aKG+alphaACE*kBM_ACE_aKG)*aKG/default0;
return

function [result] = Function_for_vBM_OAA(OAA,alphaACE,alphaGLC,default0,kBM_ACE_OAA,kBM_GLC_OAA)
global time
result = (alphaGLC*kBM_GLC_OAA+alphaACE*kBM_ACE_OAA)*OAA/default0;
return

function [result] = Function_for_vD_ACEex(ACEex,D,default0)
global time
result = D*ACEex/default0;
return

function [result] = Function_for_vD_ICIT(ICIT,default0,mu)
global time
result = mu*ICIT/default0;
return

function [result] = Function_for_vBM_SUC(SUC,alphaACE,alphaGLC,default0,kBM_ACE_SUC,kBM_GLC_SUC)
global time
result = (alphaGLC*kBM_GLC_SUC+alphaACE*kBM_ACE_SUC)*SUC/default0;
return

function [result] = Function_for_vE_PDH(AcCoA,CoA,KPDH_AcCoA_m,KPDH_CoA_m,KPDH_NADH_m,KPDH_NAD_m,KPDH_PYR_m,KPDH_i,NAD,NADH,PDH,PYR,default0,kPDH_cat)
global time
result = PDH*kPDH_cat*(1/(1+KPDH_i*NADH/NAD))*(PYR/KPDH_PYR_m)*(NAD/KPDH_NAD_m)*(CoA/KPDH_CoA_m)/((1+PYR/KPDH_PYR_m)*(1+NAD/KPDH_NAD_m+NADH/KPDH_NADH_m)*(1+CoA/KPDH_CoA_m+AcCoA/KPDH_AcCoA_m))/default0;
return

function [result] = Function_for_vE_ICDH(ICDH,ICIT,KICDH_ICIT,KICDH_PEP,LICDH,PEP,default0,kICDH_cat,nICDH)
global time
result = ICDH*kICDH_cat*ICIT/KICDH_ICIT*(1+ICIT/KICDH_ICIT)^(nICDH-1)/((1+ICIT/KICDH_ICIT)^(nICDH)+LICDH*(1+PEP/KICDH_PEP)^(nICDH))/default0;
return

function [result] = Function_for_vBM_PYR(PYR,alphaACE,alphaGLC,default0,kBM_ACE_PYR,kBM_GLC_PYR)
global time
result = (alphaGLC*kBM_GLC_PYR+alphaACE*kBM_ACE_PYR)*PYR/default0;
return

function [result] = Function_for_vBM_R5P(R5P,alphaACE,alphaGLC,default0,kBM_ACE_R5P,kBM_GLC_R5P)
global time
result = (alphaGLC*kBM_GLC_R5P+alphaACE*kBM_ACE_R5P)*R5P/default0;
return

function [result] = Function_for_vE_Fum(FUM,KFum_FUM_m,KFum_eq,MAL,default0,vFum1_max,vFum2_max)
global time
result = vFum1_max*vFum2_max*(FUM-MAL/KFum_eq)/(KFum_FUM_m*vFum2_max+vFum2_max*FUM+vFum1_max*MAL/KFum_eq)/default0;
return

function [result] = Function_for_vD_AceK(AceK,default0,kdegr,mu)
global time
result = (mu+kdegr)*AceK/default0;
return

function [result] = Function_for_vBM_F6P(F6P,alphaACE,alphaGLC,default0,kBM_ACE_F6P,kBM_GLC_F6P)
global time
result = (alphaGLC*kBM_GLC_F6P+alphaACE*kBM_ACE_F6P)*F6P/default0;
return

function [result] = Function_for_vBM_FUM(FUM,alphaACE,alphaGLC,default0,kBM_ACE_FUM,kBM_GLC_FUM)
global time
result = (alphaGLC*kBM_GLC_FUM+alphaACE*kBM_ACE_FUM)*FUM/default0;
return


