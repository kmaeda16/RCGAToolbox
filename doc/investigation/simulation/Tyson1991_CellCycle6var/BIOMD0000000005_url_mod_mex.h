#include "mex.h"

const int NRSTATES = 6;
const int NRPARAMETERS = 12;
const int NRVARIABLES = 2;
const int NRREACTIONS = 9;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[6] = {
	0,0.75,0,0.25,0,0};
char *defaultICs_nonnum[1];

double defaultParam[12] = {
	0,1,1e+06,1000,200,0,0.015,0,0.6,180,0.018,1};
char *stateNames[6] = {
	"C2","CP","M","pM","Y","YP"};
char *parameterNames[12] = {
	"EmptySet","Reaction1_k6","Reaction2_k8notP","Reaction3_k9","Reaction4_k3","Reaction5_k5notP","Reaction6_k1aa","Reaction7_k2","Reaction8_k7","Reaction9_k4","Reaction9_k4prime","cell"};
char *variableNames[2] = {
	"YT","CT"};
char *variableFormulas[2] = {
	"(Y/cell)+(YP/cell)+(M/cell)+(pM/cell)","(C2/cell)+(CP/cell)+(M/cell)+(pM/cell)"};
char *reactionNames[9] = {
	"Reaction1","Reaction2","Reaction3","Reaction4","Reaction5","Reaction6","Reaction7","Reaction8","Reaction9"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
