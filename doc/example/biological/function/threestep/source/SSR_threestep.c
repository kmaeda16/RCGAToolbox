#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definitions of mathematical functions */
#include "mex.h"

#define N_VAR   8  /* number of variables */
#define N_PARAM 38 /* number of parameters */
#define T0      RCONST(0) /* initial time */
#define TINTVL  RCONST(2) /* output time factor */
#define NOUT    60        /* number of output times */

#define RTOL    RCONST(1.0e-6) /* scalar relative tolerance */
#define ATOL    RCONST(1.0e-6) /* vector absolute tolerance components */
#define MXSTEPS 2000           /* maximum number of steps to be taken by the solver in its attempt to reach the next output time */
#define MAXNEF  20             /* maximum number of error test failures permitted in attempting one step */
#define MAXCOR  3              /* maximum number of nonlinear solver iterations permitted per step */ 
#define MAXNCF  10             /* maximum number of nonlinear solver convergence failures permitted during one step */

#define FLAG_STATS 0                    /* 1: print statistics of simulation, 0: do not */
#define OUT_TIMECOURSE "TimeCourse.dat" /* name of output file */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


int Simulation(N_Vector y, N_Vector param, int neq, realtype t0, realtype tintvl, int nout, realtype reltol, N_Vector abstol,
	long int mxsteps, int maxnef, int maxcor, int maxncf, int flag_stats, N_Vector T_out, DlsMat Y_out, long int *kount,
	int (*f)(realtype, N_Vector, N_Vector, void *));
int check_flag(void *flagvalue, char *funcname, int opt);
void writeTimeCourse(N_Vector T_out, DlsMat Y_out, int neq, long int npoint, char *filename);
int threestep(realtype t, N_Vector y, N_Vector ydot, void *user_data);
void setInitConc(N_Vector y);
void setParamSet(N_Vector param);


/*-----------------------------------------------------------------------------*/

#define N_EXPDATA 16
#define N_ROW 61
#define N_COL 9

void readExpdata();
void readMatrix(char filename[], int n_row, int n_col, double Matrix[N_ROW+1][N_COL+1]);
double ExpData[N_EXPDATA+1][N_ROW+1][N_COL+1];

/*-----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES

int N_GENE = 36;
int N_CONSTRAINT = 24;
double ALLOWABLE_ERROR = -1e+10;

void setBounds(double *lb, double *ub, int n_gene);


void benchmark(double *x, int n_gene, int n_constraint, double *f, double *g)
{
	N_Vector y, param, abstol, T_out;
	DlsMat Y_out;
	long int kount;
	int i, j , k;
	double ModelInput[N_EXPDATA+1][2+1];

	  readExpdata();

	(*f) = 0;

	y = N_VNew_Serial(N_VAR);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) exit(1);
	param = N_VNew_Serial(N_PARAM);
	if (check_flag((void *)param, "N_VNew_Serial", 0)) exit(1);
	abstol = N_VNew_Serial(N_VAR); 
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) exit(1);
	T_out = N_VNew_Serial(NOUT+1);
	if (check_flag((void *)T_out, "N_VNew_Serial", 0)) exit(1);
	Y_out = NewDenseMat(NOUT+1,N_VAR);
	if (check_flag((void *)Y_out, "NewDenseMat", 0)) exit(1);

	for(i=1;i<=N_VAR;i++) Ith(abstol,i) = ATOL;

	ModelInput[ 1][1] =  0.1;     ModelInput[ 1][2] =  0.05;    /* Experiment 1 */
	ModelInput[ 2][1] =  0.46416; ModelInput[ 2][2] =  0.05;    /* Experiment 2 */
	ModelInput[ 3][1] =  2.1544;  ModelInput[ 3][2] =  0.05;    /* Experiment 3 */
	ModelInput[ 4][1] = 10;       ModelInput[ 4][2] =  0.05;    /* Experiment 4 */
	ModelInput[ 5][1] =  0.1;     ModelInput[ 5][2] =  0.13572; /* Experiment 5 */
	ModelInput[ 6][1] =  0.46416; ModelInput[ 6][2] =  0.13572; /* Experiment 6 */
	ModelInput[ 7][1] =  2.1544;  ModelInput[ 7][2] =  0.13572; /* Experiment 7 */
	ModelInput[ 8][1] = 10;       ModelInput[ 8][2] =  0.13572; /* Experiment 8 */
	ModelInput[ 9][1] =  0.1;     ModelInput[ 9][2] =  0.36840; /* Experiment 9 */
	ModelInput[10][1] =  0.46416; ModelInput[10][2] =  0.36840; /* Experiment 10 */
	ModelInput[11][1] =  2.1544;  ModelInput[11][2] =  0.36840; /* Experiment 11 */
	ModelInput[12][1] = 10;       ModelInput[12][2] =  0.36840; /* Experiment 12 */
	ModelInput[13][1] =  0.1;     ModelInput[13][2] =  1.0;     /* Experiment 13 */
	ModelInput[14][1] =  0.46416; ModelInput[14][2] =  1.0;     /* Experiment 14 */
	ModelInput[15][1] =  2.1544;  ModelInput[15][2] =  1.0;     /* Experiment 15 */
	ModelInput[16][1] = 10;       ModelInput[16][2] =  1.0;     /* Experiment 16 */
				   
	setParamSet(param);
	for(i=1; i<=36; i++) Ith(param,2+i) = x[i];
	for(i=1; i<=16; i++){
		Ith(param,1) = ModelInput[i][1]; /* S */
		Ith(param,2) = ModelInput[i][2]; /* P */
		setInitConc(y);
		if(Simulation(y,param,N_VAR,T0,TINTVL,NOUT,RTOL,abstol,MXSTEPS,MAXNEF,MAXCOR,MAXNCF,FLAG_STATS,T_out,Y_out,&kount,threestep)){
			(*f) += 1e+20;
		}else{
			/* writeTimeCourse(T_out,Y_out,N_VAR,kount,OUT_TIMECOURSE); exit(1); */
			for(j=1; j<=N_ROW; j++){
				for(k=1; k<=N_VAR; k++){
					(*f) += pow( IJth(Y_out,j,k) - ExpData[i][j][k+1], 2.0 );
					/* printf("%e %e\n",Ith(T_out,j),ExpData[i][j][1]); */
					/* printf("%e %e\n",IJth(Y_out,j,1),ExpData[i][j][2]); */
				}
			}
			/* printf("%d %e\n",i,(*f)); */
		}
		/*
		// Check if the ExpData was correctly read
		printf("------------ %d ------------\n",i);
		for(j=1; j<=N_ROW; j++){
			for(k=1; k<=N_COL; k++) printf("%e\t",ExpData[i][j][k]);
		printf("\n");
		}
		*/
	}

	g[ 1] = fabs( log2( x[ 7] / x[ 1] ) ) - 1.0; /* | log2( V2  / V1  ) | - 1  <= 0 */
	g[ 2] = fabs( log2( x[13] / x[ 1] ) ) - 1.0; /* | log2( V3  / V1  ) | - 1  <= 0 */
	g[ 3] = fabs( log2( x[ 8] / x[ 2] ) ) - 1.0; /* | log2( Ki2 / Ki1 ) | - 1  <= 0 */
	g[ 4] = fabs( log2( x[14] / x[ 2] ) ) - 1.0; /* | log2( Ki3 / Ki1 ) | - 1  <= 0 */
	g[ 5] = fabs( log2( x[ 9] / x[ 3] ) ) - 1.0; /* | log2( ni2 / ni1 ) | - 1  <= 0 */
	g[ 6] = fabs( log2( x[15] / x[ 3] ) ) - 1.0; /* | log2( ni3 / ni1 ) | - 1  <= 0 */
	g[ 7] = fabs( log2( x[10] / x[ 4] ) ) - 1.0; /* | log2( Ka2 / Ka1 ) | - 1  <= 0 */
	g[ 8] = fabs( log2( x[16] / x[ 4] ) ) - 1.0; /* | log2( Ka3 / Ka1 ) | - 1  <= 0 */
	g[ 9] = fabs( log2( x[11] / x[ 5] ) ) - 1.0; /* | log2( na2 / na1 ) | - 1  <= 0 */
	g[10] = fabs( log2( x[17] / x[ 5] ) ) - 1.0; /* | log2( na3 / na1 ) | - 1  <= 0 */
	g[11] = fabs( log2( x[12] / x[ 6] ) ) - 1.0; /* | log2( k2  / k1  ) | - 1  <= 0 */
	g[12] = fabs( log2( x[18] / x[ 6] ) ) - 1.0; /* | log2( k3  / k1  ) | - 1  <= 0 */
	g[13] = fabs( log2( x[22] / x[19] ) ) - 1.0; /* | log2( V5  / V4  ) | - 1  <= 0 */
	g[14] = fabs( log2( x[25] / x[19] ) ) - 1.0; /* | log2( V6  / V4  ) | - 1  <= 0 */
	g[15] = fabs( log2( x[23] / x[20] ) ) - 1.0; /* | log2( K5  / K4  ) | - 1  <= 0 */
	g[16] = fabs( log2( x[26] / x[20] ) ) - 1.0; /* | log2( K6  / K4  ) | - 1  <= 0 */
	g[17] = fabs( log2( x[24] / x[21] ) ) - 1.0; /* | log2( k5  / k4  ) | - 1  <= 0 */
	g[18] = fabs( log2( x[27] / x[21] ) ) - 1.0; /* | log2( k6  / k4  ) | - 1  <= 0 */
	g[19] = fabs( log2( x[31] / x[28] ) ) - 1.0; /* | log2( kcat2 / kcat1 ) | - 1  <= 0 */
	g[20] = fabs( log2( x[34] / x[28] ) ) - 1.0; /* | log2( kcat3 / kcat1 ) | - 1  <= 0 */
	g[21] = fabs( log2( x[32] / x[29] ) ) - 1.0; /* | log2( Km3 / Km1 ) | - 1  <= 0 */
	g[22] = fabs( log2( x[35] / x[29] ) ) - 1.0; /* | log2( Km5 / Km1 ) | - 1  <= 0 */
	g[23] = fabs( log2( x[33] / x[30] ) ) - 1.0; /* | log2( Km4 / Km2 ) | - 1  <= 0 */
	g[24] = fabs( log2( x[36] / x[30] ) ) - 1.0; /* | log2( Km6 / Km2 ) | - 1  <= 0 */
	
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(param);
	N_VDestroy_Serial(abstol);
	N_VDestroy_Serial(T_out);
	DestroyMat(Y_out);

	return;
}


void benchmark_decode(double *gene, int n_gene, double *x)
{
	double lb[n_gene], ub[n_gene];
	int i;

	setBounds(lb,ub,n_gene);

	for(i=1; i<=n_gene; i++)
		x[i] = pow( 10.0, gene[i] * ( ub[i-1] - lb[i-1] ) + lb[i-1] );

	return;
}


void setBounds(double *lb, double *ub, int n_gene)
{
	/* SSR_threestep */

	int i;

	for(i=1; i<=n_gene ;i++){
		lb[i-1] = -12; ub[i-1] = 6;
	}
	lb[ 3-1] = -1; ub[ 3-1] = 1;
	lb[ 5-1] = -1; ub[ 5-1] = 1;
	lb[ 9-1] = -1; ub[ 9-1] = 1;
	lb[11-1] = -1; ub[11-1] = 1;
	lb[15-1] = -1; ub[15-1] = 1;
	lb[17-1] = -1; ub[17-1] = 1;

	return;
}


void readExpdata()
{
	char str[100];
	int i;

	for(i=1; i<=N_EXPDATA; i++){
		sprintf(str,"%s%d%s","common/pseudoexpdata",i,".dat");
		readMatrix(str, N_ROW, N_COL, ExpData[i]);
	}
}


void readMatrix(char filename[], int n_row, int n_col, double Matrix[N_ROW+1][N_COL+1])
{
	FILE *in;
	int i, j;


again:
	if( ( in = fopen(filename,"r") ) == NULL ){
		fprintf(stderr,"File error in function readMatrix\n");
		goto again; /* exit(1); */
	}else{
		for(i=1; i<=n_row; i++){
			for(j=1; j<=n_col; j++){
				fscanf(in,"%lf",&Matrix[i][j]);
			}
		}
		fclose(in);
	}
}
