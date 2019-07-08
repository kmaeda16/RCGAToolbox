#include "mex.h"

const int NRSTATES = 11;
const int NRPARAMETERS = 37;
const int NRVARIABLES = 0;
const int NRREACTIONS = 15;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[11] = {
	0.4,0.36409,0.29457,0.66667,0.57254,0.41758,1.419,0.93464,0.05,0.1,0};
char *defaultICs_nonnum[1];

double defaultParam[37] = {
	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.1,0.1,
	0.1,1,1,1,0.1,0.1,0.1,1,1,1,2,2,2,2,2,2,1};
char *stateNames[11] = {
	"E1","E2","E3","G1","G2","G3","M1","M2","P","S","VOID"};
char *parameterNames[37] = {
	"K4","K5","K6","Ka1","Ka2","Ka3","Ki1","Ki2","Ki3","Km1","Km2","Km3","Km4","Km5","Km6","V1","V2","V3","V4","V5",
	"V6","k1","k2","k3","k4","k5","k6","kcat1","kcat2","kcat3","na1","na2","na3","ni1","ni2","ni3","default0"};
char *variableNames[1];
char *variableFormulas[1];
char *reactionNames[15] = {
	"vdeg1","vdeg2","vdeg3","vdeg4","vdeg5","vdeg6","vmet1","vmet2","vmet3","vtl1","vtl2","vtl3","vts1","vts2","vts3"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
