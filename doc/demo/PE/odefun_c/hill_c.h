#include "mex.h"

const int NRSTATES = 3;
const int NRPARAMETERS = 9;
const int NRVARIABLES = 0;
const int NRREACTIONS = 4;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[3] = {
	0,0,0};
char *defaultICs_nonnum[1];

double defaultParam[9] = {
	0,5,5.5,4,0.5,0.1,0.1,0.01,1};
char *stateNames[3] = {
	"S1","S2","S3"};
char *parameterNames[9] = {
	"S4","S0","J1_Vmax","J1_n","J1_K","J2_J2_k","J3_J3_k","J0_J0_k","compart"};
char *variableNames[1];
char *variableFormulas[1];
char *reactionNames[4] = {
	"J1","J2","J3","J0"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
