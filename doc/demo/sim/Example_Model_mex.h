#include "mex.h"

const int NRSTATES = 2;
const int NRPARAMETERS = 7;
const int NRVARIABLES = 1;
const int NRREACTIONS = 3;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[2] = {
	0,0};
char *defaultICs_nonnum[1];

double defaultParam[7] = {
	0.1,1,1,1,1,1,1};
char *stateNames[2] = {
	"X1","X2"};
char *parameterNames[7] = {
	"X0","k1","k2","k3","K2","K3","rootCompartment"};
char *variableNames[1] = {
	"X12"};
char *variableFormulas[1] = {
	"(X1/rootCompartment)+(X2/rootCompartment)"};
char *reactionNames[3] = {
	"v1","v2","v3"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
