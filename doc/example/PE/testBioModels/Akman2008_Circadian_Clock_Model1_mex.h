#include "mex.h"

const int NRSTATES = 5;
const int NRPARAMETERS = 18;
const int NRVARIABLES = 1;
const int NRREACTIONS = 8;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[5] = {
	2.46246,2.71231,1.844,2.74225,0.725579};
char *defaultICs_nonnum[1];

double defaultParam[18] = {
	1.22363,5.04543,6.3958,0.885376,0.0846004,0.313846,0.161111,0.222637,0.331485,0.29484,0.13975,0.272306,0.295421,0,6,18,1,1};
char *stateNames[5] = {
	"FC","FCp","FN","FNp","MF"};
char *parameterNames[18] = {
	"vs","ki","n","vm","km","ks","vd","k1n","k2n","ksp","vdp","k1np","k2np","amp","dawn","dusk","nucleus","cytoplasm"};
char *variableNames[1] = {
	"Tot_FRQ"};
char *variableFormulas[1] = {
	"FC+FCp+FN+FNp"};
char *reactionNames[8] = {
	"MFtrn","MFdeg","FCtrl","FCdeg","FCtrs","FCptrl","FCpdeg","FCptrs"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
