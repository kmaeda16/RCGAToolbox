#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definitions of mathematical functions */

#define N_VAR   9  /* number of variables */
#define N_PARAM 8  /* number of parameters */
#define T0      RCONST(0)   /* initial time */
#define TINTVL  RCONST(0.4) /* output time factor */
#define NOUT    9000        /* number of output times */

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
int hiv(realtype t, N_Vector y, N_Vector ydot, void *user_data);
void setInitConc(N_Vector y);
void setParamSet(N_Vector param);


/*-----------------------------------------------------------------------------*/

#define N_EXPDATA 5
#define N_ROW 300
#define N_COL 3

void readExpdata();
void readMatrix(char filename[], int n_row, int n_col, double Matrix[N_ROW+1][N_COL+1]);
double ExpData[N_EXPDATA+1][N_ROW+1][N_COL+1];
double interp(double T[], double Data[], int n, double time);

/*-----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES

int N_GENE = 20;
int N_CONSTRAINT = 12;
double ALLOWABLE_ERROR = -1e+10;

void setBounds(double *lb, double *ub, int n_gene);


void benchmark(double *x, int n_gene, int n_constraint, double *f, double *g)
{
	N_Vector y, param, abstol, T_out;
	DlsMat Y_out;
	long int kount;
	int i, j;
	double ModelInput[N_EXPDATA+1][4+1];
	double TT[NOUT+2], YY[NOUT+2], eP[N_ROW+1];

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

	/* ModelInput[][1] = I, ModelInput[][2] = S, ModelInput[][3] = E, ModelInput[][4] = offset */
	ModelInput[1][1] = 0;      ModelInput[1][2] = x[ 6]; ModelInput[1][3] = x[11]; ModelInput[1][4] = x[16]; /* Experiment 1 */
	ModelInput[2][1] = 0.0015; ModelInput[2][2] = x[ 7]; ModelInput[2][3] = x[12]; ModelInput[2][4] = x[17]; /* Experiment 2 */
	ModelInput[3][1] = 0.003;  ModelInput[3][2] = x[ 8]; ModelInput[3][3] = x[13]; ModelInput[3][4] = x[18]; /* Experiment 3 */
	ModelInput[4][1] = 0.004;  ModelInput[4][2] = x[ 9]; ModelInput[4][3] = x[14]; ModelInput[4][4] = x[19]; /* Experiment 4 */
	ModelInput[5][1] = 0.004;  ModelInput[5][2] = x[10]; ModelInput[5][3] = x[15]; ModelInput[5][4] = x[20]; /* Experiment 5 */

	
	/* Known optimum */
	/*
	ModelInput[1][1] = 0;      ModelInput[1][2] = 24.637840; ModelInput[1][3] = 0.005387; ModelInput[1][4] = -0.004763; // Experiment 1
	ModelInput[2][1] = 0.0015; ModelInput[2][2] = 23.456802; ModelInput[2][3] = 0.005183; ModelInput[2][4] = -0.004950; // Experiment 2
	ModelInput[3][1] = 0.003;  ModelInput[3][2] = 27.159763; ModelInput[3][3] = 0.006000; ModelInput[3][4] = -0.017078; // Experiment 3
	ModelInput[4][1] = 0.004;  ModelInput[4][2] = 16.190568; ModelInput[4][3] = 0.004119; ModelInput[4][4] = -0.007473; // Experiment 4
	ModelInput[5][1] = 0.004;  ModelInput[5][2] = 24.672660; ModelInput[5][3] = 0.003051; ModelInput[5][4] =  0.002483; // Experiment 5
	*/
	
	setParamSet(param);
	for(i=1; i<=5; i++) Ith(param,3+i) = x[i];
	for(i=1; i<=5; i++){
		setInitConc(y);
		Ith(y,3) = ModelInput[i][2]; /* S */
		Ith(y,4) = ModelInput[i][1]; /* I */
		Ith(y,7) = ModelInput[i][3]; /* E */
		if(Simulation(y,param,N_VAR,T0,TINTVL,NOUT,RTOL,abstol,MXSTEPS,MAXNEF,MAXCOR,MAXNCF,FLAG_STATS,T_out,Y_out,&kount,hiv)){
			(*f) += 1e+20;
		}else{
			/* writeTimeCourse(T_out,Y_out,N_VAR,kount,OUT_TIMECOURSE); exit(1); */
			for(j=1; j<=kount; j++){
				TT[j] = Ith(T_out,j);
				YY[j] = IJth(Y_out,j,2);
			}
			for(j=1; j<=N_ROW; j++){
				eP[j] = 0.024 * interp(TT,YY,kount,ExpData[i][j][2]) + ModelInput[i][4];
				/* printf("%e\t%e\t%e\n",ExpData[i][j][2],ExpData[i][j][3],eP[j]); */
				(*f) += pow( eP[j] - ExpData[i][j][3], 2.0 );
				/* (*f) += pow( (eP[j] - ExpData[i][j][3])/ExpData[i][j][3], 2.0 ); */
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
		/* if(i==1) exit(1); */
		
	}

	g[ 1] = x[1] - x[3]; /* ks - kp <= 0 */
	g[ 2] = x[4] - x[1]; /*  ki - ks <= 0 */
	g[ 3] = x[2] - x[1]; /*  kcat - ks <= 0 */
	g[ 4] = x[5] - x[2]; /*  kde - kcat <= 0 */
	g[ 5] = fabs( log2( x[ 7] / x[ 6] ) ) - 1.0; /* | log2( S0B / S0A ) | - 1  <= 0 */
	g[ 6] = fabs( log2( x[ 8] / x[ 6] ) ) - 1.0; /* | log2( S0C / S0A ) | - 1  <= 0 */
	g[ 7] = fabs( log2( x[ 9] / x[ 6] ) ) - 1.0; /* | log2( S0D / S0A ) | - 1  <= 0 */
	g[ 8] = fabs( log2( x[10] / x[ 6] ) ) - 1.0; /* | log2( S0E / S0A ) | - 1  <= 0 */
	g[ 9] = fabs( log2( x[12] / x[11] ) ) - 1.0; /* | log2( E0B / E0A ) | - 1  <= 0 */
	g[10] = fabs( log2( x[13] / x[11] ) ) - 1.0; /* | log2( E0C / E0A ) | - 1  <= 0 */
	g[11] = fabs( log2( x[14] / x[11] ) ) - 1.0; /* | log2( E0D / E0A ) | - 1  <= 0 */
	g[12] = fabs( log2( x[15] / x[11] ) ) - 1.0; /* | log2( E0E / E0A ) | - 1  <= 0 */
	
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
		x[i] = gene[i] * ( ub[i-1] - lb[i-1] ) + lb[i-1];

	return;
}


void setBounds(double *lb, double *ub, int n_gene)
{
	/* SSR_hiv */

	int i;
	
	for(i=1; i<=5; i++){
		lb[i-1] = 1e-6; ub[i-1] = 1e+6;
	}
	for(i=6; i<=10; i++){
		lb[i-1] = 10; ub[i-1] = 40;
	}
	for(i=11; i<=15; i++){
		lb[i-1] = 0.002; ub[i-1] = 0.006;
	}
	for(i=16; i<=20; i++){
		lb[i-1] = -0.1; ub[i-1] = 0.1;
	}

	return;
}


void readExpdata()
{
	char str[100];
	int i;

	for(i=1; i<=N_EXPDATA; i++){
		sprintf(str,"%s%d%s","common/expdata",i,".dat");
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


double interp(double T[], double Data[], int n, double time)
{
	int i;

	if(time < T[1] || T[n] < time) return(-1e+10);

	for(i=1; i<=n-1; i++){
		if(T[i] <= time && time <= T[i+1]){
			return( (time-T[i])/(T[i+1]-T[i])*(Data[i+1]-Data[i])+Data[i] );
		}
	}
	return(-1e+10);
}
