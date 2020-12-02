#include <stdio.h>
#include <math.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


int hiv(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_Vector param;
	double M, P, S, I, ES, EP, E, EI, EJ;
	double M_dot, P_dot, S_dot, I_dot, ES_dot, EP_dot, E_dot, EI_dot, EJ_dot;
	double kmd, kdm, kon, ks, kcat, kp, ki, kde;

	param = (N_Vector)user_data;

	M  = Ith(y,1);
	P  = Ith(y,2);
	S  = Ith(y,3);
	I  = Ith(y,4);
	ES = Ith(y,5);
	EP = Ith(y,6);
	E  = Ith(y,7);
	EI = Ith(y,8);
	EJ = Ith(y,9);
	
	kmd  = Ith(param,1); /* 0.1   // Fixed */
	kdm  = Ith(param,2); /* 0.001 // Fixed */
	kon  = Ith(param,3); /* 100   // Fixed */
	ks   = Ith(param,4); /* 46.349292 (known optimum)  // 10^-6 -- 10^+6 */
	kcat = Ith(param,5); /* 5.491365 (known optimum)   // 10^-6 -- 10^+6 */
	kp   = Ith(param,6); /* 269.804443 (known optimum) // 10^-6 -- 10^+6 */
	ki   = Ith(param,7); /* 0.000177 (known optimum)   // 10^-6 -- 10^+6 */
	kde  = Ith(param,8); /* 0.000582 (known optimum)   // 10^-6 -- 10^+6 */
	
	M_dot  = - 2 * kmd * M * M + 2 * kdm * E;
	P_dot  =   kcat * ES - kon * P * E + kp * EP;
	S_dot  = - kon * S * E + ks * ES;
	I_dot  = - kon * I * E + ki * EI;
	ES_dot = kon * S * E - ks * ES - kcat * ES;
	EP_dot = kon * P * E - kp * EP;
	E_dot  = kmd * M * M - kdm * E - kon * S * E + ks * ES
		+ kcat * ES - kon * P * E + kp * EP - kon * I * E + ki * EI;
	EI_dot = kon * I * E - ki * EI - kde * EI;
	EJ_dot = kde * EI;
	
	Ith(ydot,1) = M_dot;
	Ith(ydot,2) = P_dot;
	Ith(ydot,3) = S_dot;
	Ith(ydot,4) = I_dot;
	Ith(ydot,5) = ES_dot;
	Ith(ydot,6) = EP_dot;
	Ith(ydot,7) = E_dot;
	Ith(ydot,8) = EI_dot;
	Ith(ydot,9) = EJ_dot;

	return(0);
}


void setInitConc(N_Vector y)
{
	double I, S, E, offset;

	I = 0;      S = 24.637840; E = 0.005387; offset = -0.004763; /* Experiment 1 */
	/* I = 0.0015; S = 23.456802; E = 0.005183; offset = -0.004950; // Experiment 2 */
	/* I = 0.003;  S = 27.159763; E = 0.006000; offset = -0.017078; // Experiment 3 */
	/* I = 0.004;  S = 16.190568; E = 0.004119; offset = -0.007473; // Experiment 4 */
	/* I = 0.004;  S = 24.672660; E = 0.003051; offset =  0.002483; // Experiment 5 */
	
	Ith(y,1) = 0; /* M  */
	Ith(y,2) = 0; /* P  */
	Ith(y,3) = S; /* S  */
	Ith(y,4) = I; /* I  */
	Ith(y,5) = 0; /* ES */
	Ith(y,6) = 0; /* EP */
	Ith(y,7) = E; /* E  */
	Ith(y,8) = 0; /* EI */
	Ith(y,9) = 0; /* EJ */
}


void setParamSet(N_Vector param)
{
	Ith(param,1) = 0.1;   /* kmd // Fixed */
	Ith(param,2) = 0.001; /* kdm // Fixed */
	Ith(param,3) = 100;   /* kon // Fixed */
	Ith(param,4) = 46.349292;  /* (known optimum) // ks   // 10^-6 -- 10^+6 */
	Ith(param,5) = 5.491365;   /* (known optimum) // kcat // 10^-6 -- 10^+6 */
	Ith(param,6) = 269.804443; /* (known optimum) // kp   // 10^-6 -- 10^+6 */
	Ith(param,7) = 0.000177;   /* (known optimum) // ki   // 10^-6 -- 10^+6 */
	Ith(param,8) = 0.000582;   /* (known optimum) // kde  // 10^-6 -- 10^+6 */
}
