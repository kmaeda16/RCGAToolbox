#include "mex.h"

const int NRSTATES = 9;
const int NRPARAMETERS = 9;
const int NRVARIABLES = 0;
const int NRREACTIONS = 10;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[9] = {
	0.005387,0,0,0,0,0,0,0,24.6378};
char *defaultICs_nonnum[1];

double defaultParam[9] = {
	5.49137,0.000582,0,0.000177,0.1,100,269.804,46.3493,1};
char *stateNames[9] = {
	"E","EI","EJ","EP","ES","I","M","P","S"};
char *parameterNames[9] = {
	"kcat","kde","kdm","ki","kmd","kon","kp","ks","default0"};
char *variableNames[1];
char *variableFormulas[1];
char *reactionNames[10] = {
	"ved","vef","veid","veif","vejf","vepd","vepf","vesd","vesf","vpf"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
