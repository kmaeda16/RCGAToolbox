#include "mex.h"

const int NRSTATES = 52;
const int NRPARAMETERS = 340;
const int NRVARIABLES = 26;
const int NRREACTIONS = 129;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[52] = {
	1.11524e-10,0.422054,0.0237175,0.00107739,0.000406765,0.0200874,0.0261385,0.00367648,0.650133,0.0992524,0.0372591,0.0651811,0.000460524,0.0152583,0.913391,0.173332,0.0652304,0.012,22.2222,0,
	3.77198e-10,0.00548198,0.0154742,0.0128716,0.00736962,0.0367581,1.38744,0.417978,0.00976235,0.0187144,0.0133654,0.00254143,0.00183994,1.16312,0.144336,0.00360495,0.00967106,0.00323499,0.0173145,0.0107065,
	0.156479,0.334105,0.141549,0.193564,0.116135,0.011,0.486039,0.0010271,0.00819793,0.246959,0.391757,0.0109159};
char *defaultICs_nonnum[1];

double defaultParam[340] = {
	0.56,0.28,9.6,1.4,0.0003,0.0115,0,0.0769,0,0.509128,0.0293109,8.3e-05,0.0999956,3.03764,0.00774325,0.0180077,1.0927,0.415957,0.731925,8.19069,
	0.0658771,0.0345864,0.232436,0.0311645,0.723859,6.33002,0.296006,0.0881512,0.0742685,222.087,0.0233726,0.172622,0.0127217,0.0177911,0.271462,0.0322704,0.368586,0.00242536,1.37717,0.350113,
	11.6552,0.500041,0.347444,0.648821,1000.01,0.0879915,0.13302,0.0879908,0.771367,0.33082,0.00255028,0.225362,0.0554577,18.9482,0.0699847,0.00607658,0.192067,0.0251647,0.104745,0.449961,
	0.0200034,0.0995525,0.217878,0.799718,14.9942,0.219961,9.55731e-05,0.209919,0.492465,0.0257752,0.025628,0.202044,1.32753,1.33028,0.0242784,0.0168477,0.309937,0.330748,0.099998,0.25207,
	0.396105,0.269873,0.557468,0.381655,0.361087,0.311119,2.99104,0.00293864,2.40366,0.009166,2.34579,0.0080012,0.00698447,0.248943,0.399784,1.20561,21.7028,0.00660253,0.0127905,0.0388343,
	0.0396311,0.0415227,0.669902,0.532385,0.0701171,0.059769,0.0856307,238.935,0.250106,0.359964,6.05077,0.0207231,0.160031,0.0231651,1.93645,0.199983,0.19997,2.45964,0.200129,0.717321,
	0.0229252,9.96181,42.7993,0.00225241,9.71558e-06,0.126971,0.0385597,0.000604094,0.0006405,0.088354,0.0678373,0.523402,0.0244299,1.61647,0.699177,0.0326788,0.213774,0.209251,22.4085,0.190009,
	0.370704,0.484213,1.62946,0.220018,8.67636,1.17897,1.2001,10,0.00391074,0.0180044,0.0700101,0.533941,1.41777,1.78176,0.99931,0.00018881,4.56257,1.13147e-06,0.00566376,0.249242,
	0.00858996,0.00846341,0.0662446,0.00373129,0.03915,0.000122915,0.122193,0.00494581,0.0296522,2.91919e-08,0.0312059,5.3399e-05,0.255211,0.000275018,1.24734e-05,2.41644e-08,0.00062288,0.000117808,0.045294,0.117486,
	7.59785e+07,1.11817e+07,126.46,33725.8,197779,1.26952e+06,3.8947e+06,4.03223e-80,1000.15,1243.52,2.6,0.083,0.0021,0.12,3.30385,1.50583,6.66e-05,10,0.0334723,0.0133738,
	0.000852337,0.00323703,2.00107,5.15e-07,1.30699e-05,9.19676e+15,8.85104e+12,116474,164.317,284.114,321.291,143.731,118.486,734.492,7011.8,202.218,12876.3,240.315,150.808,324.664,
	2656.7,603.434,966.423,1091.82,52.0836,375.616,19606.2,708.44,707.651,247.339,2467.94,4673.85,792390,1.60675e+06,2.41475e+06,681983,681983,8.79789e+07,1.98397e+06,211215,
	1.09342e+06,328277,328277,13791.5,1.65343e+06,6.59379e+07,42555.9,2.60105e+06,1.49637e+10,1.0471e+07,1070.78,17530.3,22311.9,22311.9,2.04051e+08,3.36933,0.265324,4.85939,12994.3,2,
	2,1,4,2,4,1.92774,1,4,3,2,4,2.31,0.74,1.07,0.74,0.74,7.5,7.49978,7.15931,13.1433,
	3.42317,564,193585,533045,13.2427,582.202,944.737,53780.5,4679.47,10119.9,3.62077e+06,45257.5,5915.4,42445.8,23336.1,42057.1,7656.02,285535,0.0596534,7.53474e-05,
	1.59563e-05,0.00262305,0.000592461,0,9.91179,0,0.0173187,0.0150909,0,0.000701946,0,0.0540414,0,0,0.0156125,0.0242728,0,0.00139121,0.184677,0.0505062,
	0,0.0159016,0.00401129,0.158286,0,0.0110092,0,9.45249e-06,0.00156527,0.00272437,0.0754715,0.113449,0,8.98315e-05,0.0038128,0.357936,0,0.0282708,0,1};
char *stateNames[52] = {
	"ACEex","AcCoA","AcP","AceK","Acs","CS","E4P","EIIAP","F6P","FBP","FUM","Fba","Fbp","Fum","G6P","GAP","GAPDH","GLC","GLCex","GLCfeed",
	"GOX","Glk","ICDH","ICDHP","ICIT","Icl","KDPG","MAL","MDH","MS","Mez","OAA","PDH","PEP","PYR","Pck","Pfk","Ppc","Pps","Pyk",
	"R5P","RU5P","S7P","SDH","SUC","X","X5P","aKG","aKGDH","cAMP","sixPG","sixPGL"};
char *parameterNames[340] = {
	"ADP","AMP","ATP","CoA","Cratotal","Crptotal","D","EIIAtotal","EXTERNAL","Factor_aceB","Factor_aceK","IclRtotal","K6PGDH_6PG","K6PGDH_ATP_inh","K6PGDH_NADP","K6PGDH_NADPH_inh","KAceK_3PG","KAceK_GOX","KAceK_ICDH","KAceK_ICDHP",
	"KAceK_ICIT","KAceK_OAA","KAceK_PEP","KAceK_PYR","KAceK_aKG","KAck_ACE_m","KAck_ADP_m","KAck_ATP_m","KAck_AcP_m","KAck_eq","KAcs_ACE","KCS_AcCoA","KCS_OAA","KCS_OAA_AcCoA","KCS_aKG","KCraFBP","KCrpcAMP","KCya_EIIAP","KEda_GAP_m","KEda_KDPG_m",
	"KEda_PYR_m","KEda_eq","KEdd_6PG_m","KEdd_KDPG_m","KEdd_eq","KFba_DHAP","KFba_FBP","KFba_GAP","KFba_GAP_inh","KFba_eq","KFbp_FBP","KFbp_PEP","KFum_FUM_m","KFum_eq","KG6PDH_G6P","KG6PDH_NADP","KG6PDH_NADPH_g6pinh","KG6PDH_NADPH_nadpinh","KGAPDH_GAP","KGAPDH_NAD",
	"KGAPDH_NADH","KGAPDH_PGP","KGAPDH_eq","KGlk_ATP_m","KGlk_G6P_i","KGlk_GLC_m","KICDH_ICIT","KICDH_PEP","KIcl_3PG","KIcl_ICIT","KIcl_PEP","KIcl_aKG","KMDH_MAL_I","KMDH_MAL_m","KMDH_NADH_I","KMDH_NADH_m","KMDH_NAD_I","KMDH_NAD_II","KMDH_NAD_m","KMDH_OAA_I",
	"KMDH_OAA_II","KMDH_OAA_m","KMDH_eq","KMS_AcCoA","KMS_GOX","KMS_GOX_AcCoA","KMez_AcCoA","KMez_MAL","KMez_cAMP","KNonPTS_I","KNonPTS_S","KPDH_AcCoA_m","KPDH_CoA_m","KPDH_NADH_m","KPDH_NAD_m","KPDH_PYR_m","KPDH_i","KPTS_EIIA","KPTS_GLC","KPck_ADP_i",
	"KPck_ATP_I","KPck_ATP_i","KPck_OAA","KPck_OAA_I","KPck_PEP","KPck_PEP_i","KPdhRPYR","KPfk_ADP_a","KPfk_ADP_b","KPfk_ADP_c","KPfk_AMP_a","KPfk_AMP_b","KPfk_ATP_s","KPfk_F6P_s","KPfk_PEP","KPgi_F6P","KPgi_F6P_6pginh","KPgi_G6P","KPgi_G6P_6pginh","KPgi_eq",
	"KPgl_6PGL_m","KPgl_6PG_m","KPgl_eq","KPgl_h1","KPgl_h2","KPpc_FBP","KPpc_PEP","KPps_PEP","KPps_PYR","KPta_AcCoA_i","KPta_AcP_i","KPta_AcP_m","KPta_CoA_i","KPta_Pi_i","KPta_Pi_m","KPta_eq","KPyk_ADP","KPyk_AMP","KPyk_ATP","KPyk_FBP",
	"KPyk_PEP","KR5PI_eq","KRu5P_eq","KSDH_SUC_m","KSDH_eq","KTal_eq","KTktA_eq","KTktB_eq","KaKGDH_CoA_m","KaKGDH_NADH_I","KaKGDH_NAD_m","KaKGDH_SUC_I","KaKGDH_Z","KaKGDH_aKG_I","KaKGDH_aKG_m","KaceBAK_Cra","KaceBAK_Crp","KaceBAK_DNA","KaceBAK_GOX","KaceBAK_PYR",
	"KaceBAK_PYRprime","Kacs_Crp","KcAMPdegr_cAMP","KfbaA_Cra","KfbaA_Crp","Kfbp_Cra","KfumABC_Crp","KgapA_Cra","KgapA_Crp","Kglk_Cra","KgltA_Crp","KicdA_Cra","Kmdh_Crp","KpckA_Cra","Kpdh_PdhR","KpfkA_Cra","KppsA_Cra","KpykF_Cra","KsdhCDAB_Crp","KsucAB_Crp",
	"LAceK","LFbp","LICDH","LIcl","LMez","LPfk","LPpc","LPps","LPyk","LaceBAK","NAD","NADH","NADP","NADPH","POratio","POratio_prime","PdhRtotal","Pi","SS_Mez_ACE","SS_Mez_GLC",
	"SS_Ppc_ACE","SS_Ppc_GLC","VFba_blf","aceBAK_DNA","kATP","kAceKki_cat","kAceKph_cat","kAcs_cat","kBM_ACE_AcCoA","kBM_ACE_E4P","kBM_ACE_F6P","kBM_ACE_FUM","kBM_ACE_G6P","kBM_ACE_GAP","kBM_ACE_OAA","kBM_ACE_PEP","kBM_ACE_PYR","kBM_ACE_R5P","kBM_ACE_SUC","kBM_ACE_aKG",
	"kBM_GLC_AcCoA","kBM_GLC_E4P","kBM_GLC_F6P","kBM_GLC_FUM","kBM_GLC_G6P","kBM_GLC_GAP","kBM_GLC_OAA","kBM_GLC_PEP","kBM_GLC_PYR","kBM_GLC_R5P","kBM_GLC_SUC","kBM_GLC_aKG","kCS_cat","kFba_cat","kFbp_cat","kFum1_cat","kFum2_cat","kGAPDH_cat","kGlk_cat","kICDH_cat",
	"kIcl_cat","kMDH1_cat","kMDH2_cat","kMS_cat","kMez_cat","kPDH_cat","kPTS1","kPck_cat","kPfk_cat","kPpc_cat","kPps_cat","kPyk_cat","kSDH1_cat","kSDH2_cat","kaKGDH_cat","kaceBAK_cat_IclR","kdegr","kexpr","kmPTS1","nAceK",
	"nCraFBP","nCrpcAMP","nFbp","nICDH","nIcl","nMez","nPdhRPYR","nPfk","nPpc","nPps","nPyk","nacs","nfumABC","ngltA","nsdhCDAB","nsucAB","pH","pH_Eda_m","pH_Edd_m","pK_Eda",
	"pK_Edd","rho","v6PGDH_max","vAck_max","vCya_max","vEda_max","vEdd_max","vG6PDH_max","vNonPTS_max","vPTS4_max","vPgi_max","vPgl_max","vPta_max","vR5PI_max","vRu5P_max","vTal_max","vTktA_max","vTktB_max","vaceBAK_Cra_bound","vaceBAK_Cra_unbound",
	"vaceBAK_Crp_bound","vaceBAK_Crp_unbound","vacs_Crp_bound","vacs_Crp_unbound","vcAMPdegr_max","vfbaA_Cra_bound","vfbaA_Cra_unbound","vfbaA_Crp_bound","vfbaA_Crp_unbound","vfbp_Cra_bound","vfbp_Cra_unbound","vfumABC_Crp_bound","vfumABC_Crp_unbound","vgapA_Cra_bound","vgapA_Cra_unbound","vgapA_Crp_bound","vgapA_Crp_unbound","vglk_Cra_bound","vglk_Cra_unbound","vgltA_Crp_bound",
	"vgltA_Crp_unbound","vicdA_Cra_bound","vicdA_Cra_unbound","vmdh_Crp_bound","vmdh_Crp_unbound","vpckA_Cra_bound","vpckA_Cra_unbound","vpdh_PdhR_bound","vpdh_PdhR_unbound","vpfkA_Cra_bound","vpfkA_Cra_unbound","vppsA_Cra_bound","vppsA_Cra_unbound","vpykF_Cra_bound","vpykF_Cra_unbound","vsdhCDAB_Crp_bound","vsdhCDAB_Crp_unbound","vsucAB_Crp_bound","vsucAB_Crp_unbound","default0"};
char *variableNames[26] = {
	"CraFBP","B","A","Cra","CrpcAMP","Crp","H","EIIA","alphaGLC","alphaACE","SS_Ppc","SS_Mez","OP_NADH","OP_FADH2","PdhRPYR","PdhR","QEda_pH","QEdd_pH","vMDH1_max","vFum2_max",
	"vSDH2_max","vSDH1_max","vATP","mu","vMDH2_max","vFum1_max"};
char *variableFormulas[26] = {
	"Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))","1.0+ADP/KPfk_ADP_a+AMP/KPfk_AMP_a","1.0+PEP/KPfk_PEP+ADP/KPfk_ADP_b+AMP/KPfk_AMP_b","Cratotal-CraFBP","Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))","Crptotal-CrpcAMP","pow(10.0,-pH)*1000.0","EIIAtotal-EIIAP","GLCex/(GLCex+KPTS_GLC)","ACEex/(ACEex+KAcs_ACE)*(1.0-alphaGLC)","alphaGLC*SS_Ppc_GLC+alphaACE*SS_Ppc_ACE","alphaGLC*SS_Mez_GLC+alphaACE*SS_Mez_ACE","(vE_GAPDH+vE_PDH+vE_aKGDH+vE_MDH)*POratio","vE_SDH*POratio_prime","PdhRtotal*pow(PYR,nPdhRPYR)/(pow(PYR,nPdhRPYR)+pow(KPdhRPYR,nPdhRPYR))","PdhRtotal-PdhRPYR","1.0+2.0*pow(10.0,pH_Eda_m-pK_Eda)/(1.0+pow(10.0,pH-pK_Eda)+pow(10.0,2.0*pH_Eda_m-pH-pK_Eda))","1.0+2.0*pow(10.0,pH_Edd_m-pK_Edd)/(1.0+pow(10.0,pH-pK_Edd)+pow(10.0,2.0*pH_Edd_m-pH-pK_Edd))","MDH*kMDH1_cat","Fum*kFum2_cat",
	"SDH*kSDH2_cat","SDH*kSDH1_cat","OP_NADH+OP_FADH2-vE_Glk-vE_Pfk+vE_GAPDH+vE_Pyk-vE_Pps+vE_Ack-vE_Acs+vE_aKGDH-vE_Pck-vE_AceKki-vE_Cya","kATP*vATP","MDH*kMDH2_cat","Fum*kFum1_cat"};
char *reactionNames[129] = {
	"vBM_AcCoA","vBM_E4P","vBM_F6P","vBM_FUM","vBM_G6P","vBM_GAP","vBM_OAA","vBM_PEP","vBM_PYR","vBM_R5P","vBM_SUC","vBM_aKG","vD_6PG","vD_6PGL","vD_ACEex","vD_AcCoA","vD_AcP","vD_AceK","vD_Acs","vD_CS",
	"vD_E4P","vD_F6P","vD_FBP","vD_FUM","vD_Fba","vD_Fbp","vD_Fum","vD_G6P","vD_GAP","vD_GAPDH","vD_GLC","vD_GLCex","vD_GLCfeed","vD_GOX","vD_Glk","vD_ICDH","vD_ICDHP","vD_ICIT","vD_Icl","vD_KDPG",
	"vD_MAL","vD_MDH","vD_MS","vD_Mez","vD_OAA","vD_PDH","vD_PEP","vD_PYR","vD_Pck","vD_Pfk","vD_Ppc","vD_Pps","vD_Pyk","vD_R5P","vD_RU5P","vD_S7P","vD_SDH","vD_SUC","vD_X","vD_X5P",
	"vD_aKG","vD_aKGDH","vD_cAMP","vE_6PGDH","vE_AceKki","vE_AceKph","vE_Ack","vE_Ack_medium","vE_Acs","vE_Acs_medium","vE_CS","vE_Cya","vE_Eda","vE_Edd","vE_Fba","vE_Fbp","vE_Fum","vE_G6PDH","vE_GAPDH","vE_Glk",
	"vE_ICDH","vE_Icl","vE_MDH","vE_MS","vE_Mez","vE_PDH","vE_Pck","vE_Pfk","vE_Pgi","vE_Pgl","vE_Ppc","vE_Pps","vE_Pta","vE_Pyk","vE_R5PI","vE_Ru5P","vE_SDH","vE_Tal","vE_TktA","vE_TktB",
	"vE_aKGDH","vE_cAMPdegr","vG_aceA","vG_aceB","vG_aceK","vG_acs","vG_fbaA","vG_fbp","vG_fumABC","vG_gapA","vG_glk","vG_gltA","vG_icdA","vG_maeB","vG_mdh","vG_pckA","vG_pdh","vG_pfkA","vG_ppc","vG_ppsA",
	"vG_pykF","vG_sdhCDAB","vG_sucAB","vNonPTS","vNonPTS_medium","vPTS1","vPTS4","vPTS4_medium","vgrowth"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
