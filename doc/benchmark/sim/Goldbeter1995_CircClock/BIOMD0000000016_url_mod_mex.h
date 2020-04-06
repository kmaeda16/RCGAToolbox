#include "mex.h"

const int NRSTATES = 6;
const int NRPARAMETERS = 21;
const int NRVARIABLES = 1;
const int NRREACTIONS = 10;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[6] = {
	0,0.1,0.25,0.25,0.25,0.25};
char *defaultICs_nonnum[1];

double defaultParam[21] = {
	0.76,1,4,0.38,3.2,2,1.58,2,5,2,2.5,2,1.9,1.3,0.5,0.65,0.95,0.2,1e-15,1e-15,
	1e-15};
char *stateNames[6] = {
	"EmptySet","M","P0","P1","P2","Pn"};
char *parameterNames[21] = {
	"rM_Vs","rM_KI","rM_n","rTL_ks","rP01_V1","rP01_K1","rP10_V2","rP10_K2","rP12_V3","rP12_K3","rP21_V4","rP21_K4","rP2n_k1","rPn2_k2","rmRNAd_Km","rmRNAd_Vm","rVd_Vd","rVd_Kd","default0","CYTOPLASM",
	"compartment_0000004"};
char *variableNames[1] = {
	"Pt"};
char *variableFormulas[1] = {
	"P0+P1+P2+Pn"};
char *reactionNames[10] = {
	"rM","rTL","rP01","rP10","rP12","rP21","rP2n","rPn2","rmRNAd","rVd"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
