#include <stdio.h>
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

void PrintFinalStats(void *cvode_mem);
int check_flag(void *flagvalue, char *funcname, int opt);

int Simulation(N_Vector y, N_Vector param, int neq, realtype t0, realtype tintvl, int nout, realtype reltol, N_Vector abstol,
	long int mxsteps, int maxnef, int maxcor, int maxncf, int flag_stats, N_Vector T_out, DlsMat Y_out, long int *kount,
	int (*f)(realtype, N_Vector, N_Vector, void *))
{
	realtype t, tout;
	int flag, i, iout;
	void *cvode_mem;

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	flag = CVodeInit(cvode_mem, f, t0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);

	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	flag = CVDense(cvode_mem, neq);
	if (check_flag(&flag, "CVDense", 1)) return(1);

	flag = CVodeSetUserData(cvode_mem, param);
	if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

	flag = CVodeSetMaxNumSteps(cvode_mem, mxsteps);
	if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

	flag = CVodeSetMaxErrTestFails(cvode_mem, maxnef);
	if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

	flag = CVodeSetMaxNonlinIters(cvode_mem, maxcor);
	if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

	flag = CVodeSetMaxConvFails(cvode_mem, maxncf);
	if (check_flag(&flag, "CVodeSetMaxConvFails", 1)) return(1);

	flag = CVodeSetMaxHnilWarns(cvode_mem, -1);
	if (check_flag(&flag, "CVodeSetMaxHnilWarns", 1)) return(1);

	flag = CVodeSetNoInactiveRootWarn(cvode_mem);
	if (check_flag(&flag, "CVodeSetNoInactiveRootWarn", 1)) return(1);

	Ith(T_out,0+1) = t = t0;
	for(i=1;i<=neq;i++) IJth(Y_out,0+1,i) = Ith(y,i);
	iout = 1; tout = t0 + tintvl;
	while(1) {
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		Ith(T_out,iout+1) = t;
		for(i=1;i<=neq;i++) IJth(Y_out,iout+1,i) = Ith(y,i);
		if (iout == nout) break;
		if (check_flag(&flag, "CVode", 1)) {
			*kount = iout + 1;
			if(flag_stats) PrintFinalStats(cvode_mem);
			CVodeFree(&cvode_mem);
			return(1);
		}
		if (flag == CV_SUCCESS) {
			iout++;
			tout += tintvl;
		}
	}
	*kount = iout + 1;
	if(flag_stats) PrintFinalStats(cvode_mem);
	CVodeFree(&cvode_mem);
	return(0);
}


void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}


/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}


void writeTimeCourse(N_Vector T_out, DlsMat Y_out, int neq, long int npoint, char *filename)
{
	int i,j;
	FILE *out;

	if((out= fopen(filename,"w"))==NULL){
		printf("File error in function writeTimeCourse\n");
	}else{
		fprintf(out,"Time\t");
		for(j=1;j<=neq;j++){
			fprintf(out,"%d\t",j);
		}
		fprintf(out,"\n");
		for(i=1;i<=npoint;i++){
			fprintf(out,"%e\t",Ith(T_out,i));
			for(j=1;j<=neq;j++){
				fprintf(out,"%e\t",IJth(Y_out,i,j));
			}
			fprintf(out,"\n");
		}
		fclose(out);
	}

}
