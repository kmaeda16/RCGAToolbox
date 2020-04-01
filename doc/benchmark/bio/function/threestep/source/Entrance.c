#include "mex.h"
#include "matrix.h"

#define N_GENE 36 /* 13 */
#define N_CONSTRAINT 24 /* 9 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double x[N_GENE+1], f, g[N_CONSTRAINT+1];
	double *pointer;
	int i;

	if ( nrhs != 1 ) {
        mexErrMsgIdAndTxt("MATLAB:arrayFillGetPr:rhs","This function takes one input arguments.");
    }

	pointer =  mxGetPr(prhs[0]);
	for(i=1; i<=N_GENE; i++) x[i] = pointer[i-1];

	benchmark(x, N_GENE, N_CONSTRAINT, &f, g);

	plhs[0] = mxCreateDoubleScalar(mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, N_CONSTRAINT, mxREAL);

	pointer = mxGetPr(plhs[0]);
	(*pointer) = f;

	pointer = mxGetPr(plhs[1]);
	for(i=1; i<=N_CONSTRAINT; i++) pointer[i-1] = g[i];

	return;
}
