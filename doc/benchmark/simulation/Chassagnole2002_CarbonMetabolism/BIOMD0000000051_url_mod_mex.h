#include "mex.h"

const int NRSTATES = 18;
const int NRPARAMETERS = 137;
const int NRVARIABLES = 7;
const int NRREACTIONS = 48;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[18] = {
	2.67,2,3.48,2.67,0.6,0.653,0.808,0.272,0.276,0.218,0.098,0.138,0.398,0.167,0.008,2.13,0.399,0.111};
char *defaultICs_nonnum[1];

double defaultParam[137] = {
	7829.78,3082.3,0.01,245.3,3.66,2.15,650.988,0.1725,2.9,0.266,0.2,0.2,0.839824,0.196,1.038,0.0136,1.3802,14.4,6.43,0.0246,
	0.01,1840.58,0.123,4.14,0.325,3.26,3.89,3.2,128,19.1,5.62907e+06,11.1,10.8716,1.05,9.47338,1.2,86.5586,10,0.00043711,17.4146,
	0.144,1.75,0.088,2,0.088,0.6,921.594,0.63,0.683,1.04e-05,0.252,1.09,68.6747,1.39,2.8,0.3,0.001037,0.0116204,1,3021.77,
	1934.4,0.185,0.653,0.0468,0.473,0.0257121,1,89.0497,0.188,0.2,0.369,330.448,6.73,0.1,0.135,0.0611315,0.31,4,1000,22.5,
	0.19,0.2,0.26,0.107021,0.7,4.21,4.07,0.019539,1,0.0736186,1,0.107953,2.6,2.2,0.035,0.0053,6.05953,3.68,1159,0.0022627,
	16.2324,37.5,0.0506,0.0138,208,4.83841,4,6.73903,1.4,0.0129005,0.1,0.00752546,0.119,1.2,4.42,3.2,2.78e-05,2.78e-05,2.78e-05,2.78e-05,
	2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,2.78e-05,110.96,1,1};
char *stateNames[18] = {
	"cpep","cglcex","cg6p","cpyr","cf6p","cg1p","cpg","cfdp","csed7p","cgap","ce4p","cxyl5p","crib5p","cdhap","cpgp","cpg3","cpg2","cribu5p"};
char *parameterNames[137] = {
	"vPTS_rmaxPTS","vPTS_KPTSa1","vPTS_KPTSa2","vPTS_KPTSa3","vPTS_nPTSg6p","vPTS_KPTSg6p","vPGI_rmaxPGI","vPGI_KPGIeq","vPGI_KPGIg6p","vPGI_KPGIf6p","vPGI_KPGIf6ppginh","vPGI_KPGIg6ppginh","vPGM_rmaxPGM","vPGM_KPGMeq","vPGM_KPGMg6p","vPGM_KPGMg1p","vG6PDH_rmaxG6PDH","vG6PDH_KG6PDHg6p","vG6PDH_KG6PDHnadphg6pinh","vG6PDH_KG6PDHnadp",
	"vG6PDH_KG6PDHnadphnadpinh","vPFK_rmaxPFK","vPFK_KPFKatps","vPFK_KPFKadpc","vPFK_KPFKf6ps","vPFK_KPFKpep","vPFK_KPFKadpb","vPFK_KPFKampb","vPFK_KPFKadpa","vPFK_KPFKampa","vPFK_LPFK","vPFK_nPFK","vTA_rmaxTA","vTA_KTAeq","vTKA_rmaxTKa","vTKA_KTKaeq","vTKB_rmaxTKb","vTKB_KTKbeq","vMURSyNTH_rmaxMurSynth","vALDO_rmaxALDO",
	"vALDO_kALDOeq","vALDO_kALDOfdp","vALDO_kALDOgap","vALDO_VALDOblf","vALDO_kALDOdhap","vALDO_kALDOgapinh","vGAPDH_rmaxGAPDH","vGAPDH_KGAPDHeq","vGAPDH_KGAPDHgap","vGAPDH_KGAPDHpgp","vGAPDH_KGAPDHnad","vGAPDH_KGAPDHnadh","vTIS_rmaxTIS","vTIS_kTISeq","vTIS_kTISdhap","vTIS_kTISgap","vTRPSYNTH_rmaxTrpSynth","vG3PDH_rmaxG3PDH","vG3PDH_KG3PDHdhap","vPGK_rmaxPGK",
	"vPGK_KPGKeq","vPGK_KPGKadp","vPGK_KPGKatp","vPGK_KPGKpgp","vPGK_KPGKpg3","vsersynth_rmaxSerSynth","vsersynth_KSerSynthpg3","vrpGluMu_rmaxPGluMu","vrpGluMu_KPGluMueq","vrpGluMu_KPGluMupg3","vrpGluMu_KPGluMupg2","vENO_rmaxENO","vENO_KENOeq","vENO_KENOpg2","vENO_KENOpep","vPK_rmaxPK","vPK_KPKpep","vPK_nPK","vPK_LPK","vPK_KPKatp",
	"vPK_KPKfdp","vPK_KPKamp","vPK_KPKadp","vpepCxylase_rmaxpepCxylase","vpepCxylase_KpepCxylasefdp","vpepCxylase_npepCxylasefdp","vpepCxylase_KpepCxylasepep","vSynth1_rmaxSynth1","vSynth1_KSynth1pep","vSynth2_rmaxSynth2","vSynth2_KSynth2pyr","vDAHPS_rmaxDAHPS","vDAHPS_nDAHPSe4p","vDAHPS_nDAHPSpep","vDAHPS_KDAHPSe4p","vDAHPS_KDAHPSpep","vPDH_rmaxPDH","vPDH_nPDH","vPDH_KPDHpyr","vMethSynth_rmaxMetSynth",
	"vPGDH_rmaxPGDH","vPGDH_KPGDHpg","vPGDH_KPGDHnadp","vPGDH_KPGDHnadphinh","vPGDH_KPGDHatpinh","vR5PI_rmaxR5PI","vR5PI_KR5PIeq","vRu5P_rmaxRu5P","vRu5P_KRu5Peq","vPPK_rmaxRPPK","vPPK_KRPPKrib5p","vG1PAT_rmaxG1PAT","vG1PAT_KG1PATfdp","vG1PAT_nG1PATfdp","vG1PAT_KG1PATatp","vG1PAT_KG1PATg1p","vG6P_mu","vf6P_mu","vfdP_mu","vGAP_mu",
	"vDHAP_mu","vPGP_mu","vPG3_mu","vpg2_mu","vPEP_mu","vRibu5p_mu","vRIB5P_mu","vXYL5P_mu","vSED7P_mu","vpyr_mu","vPG_mu","vE4P_mu","vGLP_mu","vEXTER_Dil","vEXTER_cfeed","extracellular","cytosol"};
char *variableNames[7] = {
	"catp","cadp","camp","cnadph","cnadp","cnadh","cnad"};
char *variableFormulas[7] = {
	"4.27-4.163*(time/(0.657+1.43*time+0.0364*pow(time,2.0)))","0.582+1.73*pow(2.731,-0.15*time)*(0.12*time+0.000214*pow(time,3.0))","0.123+7.25*(time/(7.25+1.47*time+0.17*pow(time,2.0)))+1.073/(1.29+8.05*time)","0.062+0.332*pow(2.718,-0.464*time)*(0.0166*pow(time,1.58)+0.000166*pow(time,4.73)+0.1312*pow(10.0,-9.0)*pow(time,7.89)+0.1362*pow(10.0,-12.0)*pow(time,11.0)+0.1233*pow(10.0,-15.0)*pow(time,14.2))","0.159-0.00554*(time/(2.8-0.271*time+0.01*pow(time,2.0)))+0.182/(4.82+0.526*time)","0.0934+0.00111*pow(2.371,-0.123*time)*(0.844*time+0.104*pow(time,3.0))","1.314+1.314*pow(2.73,-0.0435*time-0.342)-(time+7.871)*(pow(2.73,-0.0218*time-0.171)/(8.481+time))"};
char *reactionNames[48] = {
	"vPTS","vPGI","vPGM","vG6PDH","vPFK","vTA","vTKA","vTKB","vMURSyNTH","vALDO","vGAPDH","vTIS","vTRPSYNTH","vG3PDH","vPGK","vsersynth","vrpGluMu","vENO","vPK","vpepCxylase",
	"vSynth1","vSynth2","vDAHPS","vPDH","vMethSynth","vPGDH","vR5PI","vRu5P","vPPK","vG1PAT","vG6P","vf6P","vfdP","vGAP","vDHAP","vPGP","vPG3","vpg2","vPEP","vRibu5p",
	"vRIB5P","vXYL5P","vSED7P","vpyr","vPG","vE4P","vGLP","vEXTER"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
