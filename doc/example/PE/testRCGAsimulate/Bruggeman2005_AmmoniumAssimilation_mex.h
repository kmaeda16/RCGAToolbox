#include "mex.h"

const int NRSTATES = 21;
const int NRPARAMETERS = 96;
const int NRVARIABLES = 4;
const int NRREACTIONS = 16;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[21] = {
	0.003,0.5,0,0.05,1,0,0,0.01,0.014,0,0.05,0.2,0.15,1,0.05,1,2.685,2.685,0.1,1,
	0.1};
char *defaultICs_nonnum[1];

double defaultParam[96] = {
	10,0.0006,137,0.07,0.0018,0.003,0.0035,0.04,0.1135,0.0006,5.5,0.0023,8.4,0.07,1e-22,0.5166,0.5974,0.0387,0.5,1.052e-05,
	0.9714,0.001703,1e-22,2.766,3.323,0.2148,1e-22,1e-22,1e-22,0.02316,0.8821,8.491,0.8791,0.5,2.274e-06,0.04444,1.805e-05,0.0002015,360,0.32,
	1.1,10,0.04,0.042,1290,2.5,85,0.175,0.007,0.0015,11,0.0037,0.65,600,10,2.3667,0.1012,10.8688,1.1456,19.2166,
	460,0.35,4.1,0.1,0.0585,3.7,5.65,460,120,8,1e+10,0.5,70,2,1e+10,0.25,20,1,1e+10,0.5,
	30,0.3,1e+10,0.5,100,0.5,5.37,0.014,0.003,0.005,0.15,0.15,0.025,0.15,0.15,1};
char *stateNames[21] = {
	"PII","UTP","PIIUMP","PPi","GLN","PIIUMP2","PIIUMP3","UMP","GS","AMP","NH4","KG","NADPH","GLU","NADP","AZGLU","ATP","ADP","AZglu","AZGLN",
	"AZgln"};
char *parameterNames[96] = {
	"P_i","UT","kcatut","Kglnut","Kutipii","Kutpii","Kutpiiump","Kututp","Kutippi","UR","kcatur","Kurpiiump","Kurump","Kglnur","a1","b1","c1","d1","Vad","Kadpiikg",
	"Kadgln","Kadgs","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","Vdead","Kdeadpiikg","Kdeadgln","Kdeadpiiu","Kdeadgsa","Vgdh","Kgdhkg",
	"Kgdhnh","Kgdhglu","Kgdhnadph","Kgdhnadp","Keqgdh","Kgdhazglu","Vgog","Kgoggln","Kgogkg","Kgognadph","Kgogglu","Kgognadp","Kgogaz","Vgs","aamp","bamp","camp","damp","n1amp","n2amp",
	"Kgseq","Kgsatp","Kgsglu","Kgsnh","Kgsadp","Kgspi","Kgsgln","Keq","Vgludem","Kgludemglu","Kgludemeq","Kgludemazglu","Vglndem","Kglndemgln","Kglndemeq","Kglndemazgln","Vazglndem","Kazglndemazgln","Kazglndemeq","Kazglndemazinter",
	"Vazgludem","Kazgludemazglu","Kazgludemeq","Kazgludemazinter","Vadp","Kadp","ATPtot","GStot","PIItot","Kd1","Kd2","Kd3","Kd1piiump","Kd2piiump","Kd3piiump","compartment"};
char *variableNames[4] = {
	"vAPP_GS","nAMP","PIIKG1","PIIUMP3KG3"};
char *variableFormulas[4] = {
	"aamp*camp/((1.0+pow(12.0,n1amp)*pow(AMP/(bamp*GStot),n1amp))*(1.0+pow(12.0,n2amp)*pow(AMP/(damp*GStot),n2amp)))*Vgs","12.0*AMP/GStot","3.0*PII*KG/Kd1/(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))","PIIUMP3*pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)/(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump))"};
char *reactionNames[16] = {
	"vut1","vur1","vut2","vur2","vut3","vur3","vad","vdead","vgdh","vgog","vgs","vgludem","vazgludem","vglndem","vazglndem","vatpase"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
