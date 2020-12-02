#include <stdio.h>
#include <math.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


int threestep(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_Vector param;
	double G1, G2, G3, E1, E2, E3, M1, M2;
	double G1_dot, G2_dot, G3_dot, E1_dot, E2_dot, E3_dot, M1_dot, M2_dot;
	double S, P, /* 1 - 2 */
		V1, Ki1, ni1, Ka1, na1, k1, /*  3 -  8 */
		V2, Ki2, ni2, Ka2, na2, k2, /*  9 - 14 */
		V3, Ki3, ni3, Ka3, na3, k3, /* 15 - 20 */
		V4, K4, k4, /* 21 - 23 */
		V5, K5, k5, /* 24 - 26 */
		V6, K6, k6, /* 27 - 29 */
		kcat1, Km1, Km2, /* 30 - 32 */
		kcat2, Km3, Km4, /* 33 - 35 */
		kcat3, Km5, Km6; /* 36 - 38 */

	param = (N_Vector)user_data;

	G1 = Ith(y,1);
	G2 = Ith(y,2);
	G3 = Ith(y,3);
	E1 = Ith(y,4);
	E2 = Ith(y,5);
	E3 = Ith(y,6);
	M1 = Ith(y,7);
	M2 = Ith(y,8);

	S     = Ith(param, 1); /* 0.1  (model input) // Fixed */
	P     = Ith(param, 2); /* 0.05 (model input) // Fixed */
	V1    = Ith(param, 3); /* 1   // 10^-12 -- 10^+6 */
	Ki1   = Ith(param, 4); /* 1   // 10^-12 -- 10^+6 */
	ni1   = Ith(param, 5); /* 2   // 10^-1  -- 10^+1 */
	Ka1   = Ith(param, 6); /* 1   // 10^-12 -- 10^+6 */
	na1   = Ith(param, 7); /* 2   // 10^-1  -- 10^+1 */
	k1    = Ith(param, 8); /* 1   // 10^-12 -- 10^+6 */
	V2    = Ith(param, 9); /* 1   // 10^-12 -- 10^+6 */
	Ki2   = Ith(param,10); /* 1   // 10^-12 -- 10^+6 */
	ni2   = Ith(param,11); /* 2   // 10^-1  -- 10^+1 */
	Ka2   = Ith(param,12); /* 1   // 10^-12 -- 10^+6 */
	na2   = Ith(param,13); /* 2   // 10^-1  -- 10^+1 */
	k2    = Ith(param,14); /* 1   // 10^-12 -- 10^+6 */
	V3    = Ith(param,15); /* 1   // 10^-12 -- 10^+6 */
	Ki3   = Ith(param,16); /* 1   // 10^-12 -- 10^+6 */
	ni3   = Ith(param,17); /* 2   // 10^-1  -- 10^+1 */
	Ka3   = Ith(param,18); /* 1   // 10^-12 -- 10^+6 */
	na3   = Ith(param,19); /* 2   // 10^-1  -- 10^+1 */
	k3    = Ith(param,20); /* 1   // 10^-12 -- 10^+6 */
	V4    = Ith(param,21); /* 0.1 // 10^-12 -- 10^+6 */
	K4    = Ith(param,22); /* 1   // 10^-12 -- 10^+6 */
	k4    = Ith(param,23); /* 0.1 // 10^-12 -- 10^+6 */
	V5    = Ith(param,24); /* 0.1 // 10^-12 -- 10^+6 */
	K5    = Ith(param,25); /* 1   // 10^-12 -- 10^+6 */
	k5    = Ith(param,26); /* 0.1 // 10^-12 -- 10^+6 */
	V6    = Ith(param,27); /* 0.1 // 10^-12 -- 10^+6 */
	K6    = Ith(param,28); /* 1   // 10^-12 -- 10^+6 */
	k6    = Ith(param,29); /* 0.1 // 10^-12 -- 10^+6 */
	kcat1 = Ith(param,30); /* 1   // 10^-12 -- 10^+6 */
	Km1   = Ith(param,31); /* 1   // 10^-12 -- 10^+6 */
	Km2   = Ith(param,32); /* 1   // 10^-12 -- 10^+6 */
	kcat2 = Ith(param,33); /* 1   // 10^-12 -- 10^+6 */
	Km3   = Ith(param,34); /* 1   // 10^-12 -- 10^+6 */
	Km4   = Ith(param,35); /* 1   // 10^-12 -- 10^+6 */
	kcat3 = Ith(param,36); /* 1   // 10^-12 -- 10^+6 */
	Km5   = Ith(param,37); /* 1   // 10^-12 -- 10^+6 */
	Km6   = Ith(param,38); /* 1   // 10^-12 -- 10^+6 */
	
	G1_dot = V1 / ( 1 + pow( P / Ki1, ni1 ) + pow( Ka1 / S,  na1) ) - k1 * G1;
	G2_dot = V2 / ( 1 + pow( P / Ki2, ni2 ) + pow( Ka2 / M1, na2) ) - k2 * G2;
	G3_dot = V3 / ( 1 + pow( P / Ki3, ni3 ) + pow( Ka3 / M2, na3) ) - k3 * G3;
	E1_dot = V4 * G1 / ( K4 + G1 ) - k4 * E1;
	E2_dot = V5 * G2 / ( K5 + G2 ) - k5 * E2;
	E3_dot = V6 * G3 / ( K6 + G3 ) - k6 * E3;
	M1_dot = kcat1 * E1 * ( 1 / Km1 ) * ( S  - M1 ) / ( 1 +  S / Km1 + M1 / Km2 )
		- kcat2 * E2 * ( 1 / Km3 ) * ( M1 - M2 ) / ( 1 + M1 / Km3 + M2 / Km4 );
	M2_dot = kcat2 * E2 * ( 1 / Km3 ) * ( M1 - M2 ) / ( 1 + M1 / Km3 + M2 / Km4 )
		- kcat3 * E3 * ( 1 / Km5 ) * ( M2 - P  ) / ( 1 + M2 / Km5 + P  / Km6 );

	Ith(ydot,1) = G1_dot;
	Ith(ydot,2) = G2_dot;
	Ith(ydot,3) = G3_dot;
	Ith(ydot,4) = E1_dot;
	Ith(ydot,5) = E2_dot;
	Ith(ydot,6) = E3_dot;
	Ith(ydot,7) = M1_dot;
	Ith(ydot,8) = M2_dot;

	return(0);
}


void setInitConc(N_Vector y)
{
	Ith(y,1) = 0.66667; /* G1 */
	Ith(y,2) = 0.57254; /* G2 */
	Ith(y,3) = 0.41758; /* G3 */
	Ith(y,4) = 0.4;     /* E1 */
	Ith(y,5) = 0.36409; /* E2 */
	Ith(y,6) = 0.29457; /* E3 */
	Ith(y,7) = 1.419;   /* M1 */
	Ith(y,8) = 0.93464; /* M2 */
}


void setParamSet(N_Vector param)
{
	Ith(param, 1) = 0.1;  /* S (model input) // Fixed */
	Ith(param, 2) = 0.05; /* P (model input) // Fixed */
	Ith(param, 3) = 1;    /* V1    // 10^-12 -- 10^+6 */
	Ith(param, 4) = 1;    /* Ki1   // 10^-12 -- 10^+6 */
	Ith(param, 5) = 2;    /* ni1   // 10^-1  -- 10^+1 */
	Ith(param, 6) = 1;    /* Ka1   // 10^-12 -- 10^+6 */
	Ith(param, 7) = 2;    /* na1   // 10^-1  -- 10^+1 */
	Ith(param, 8) = 1;    /* k1    // 10^-12 -- 10^+6 */
	Ith(param, 9) = 1;    /* V2    // 10^-12 -- 10^+6 */
	Ith(param,10) = 1;    /* Ki2   // 10^-12 -- 10^+6 */
	Ith(param,11) = 2;    /* ni2   // 10^-1  -- 10^+1 */
	Ith(param,12) = 1;    /* Ka2   // 10^-12 -- 10^+6 */
	Ith(param,13) = 2;    /* na2   // 10^-1  -- 10^+1 */
	Ith(param,14) = 1;    /* k2    // 10^-12 -- 10^+6 */
	Ith(param,15) = 1;    /* V3    // 10^-12 -- 10^+6 */
	Ith(param,16) = 1;    /* Ki3   // 10^-12 -- 10^+6 */
	Ith(param,17) = 2;    /* ni3   // 10^-1  -- 10^+1 */
	Ith(param,18) = 1;    /* Ka3   // 10^-12 -- 10^+6 */
	Ith(param,19) = 2;    /* na3   // 10^-1  -- 10^+1 */
	Ith(param,20) = 1;    /* k3    // 10^-12 -- 10^+6 */
	Ith(param,21) = 0.1;  /* V4    // 10^-12 -- 10^+6 */
	Ith(param,22) = 1;    /* K4    // 10^-12 -- 10^+6 */
	Ith(param,23) = 0.1;  /* k4    // 10^-12 -- 10^+6 */
	Ith(param,24) = 0.1;  /* V5    // 10^-12 -- 10^+6 */
	Ith(param,25) = 1;    /* K5    // 10^-12 -- 10^+6 */
	Ith(param,26) = 0.1;  /* k5    // 10^-12 -- 10^+6 */
	Ith(param,27) = 0.1;  /* V6    // 10^-12 -- 10^+6 */
	Ith(param,28) = 1;    /* K6    // 10^-12 -- 10^+6 */
	Ith(param,29) = 0.1;  /* k6    // 10^-12 -- 10^+6 */
	Ith(param,30) = 1;    /* kcat1 // 10^-12 -- 10^+6 */
	Ith(param,31) = 1;    /* Km1   // 10^-12 -- 10^+6 */
	Ith(param,32) = 1;    /* Km2   // 10^-12 -- 10^+6 */
	Ith(param,33) = 1;    /* kcat2 // 10^-12 -- 10^+6 */
	Ith(param,34) = 1;    /* Km3   // 10^-12 -- 10^+6 */
	Ith(param,35) = 1;    /* Km4   // 10^-12 -- 10^+6 */
	Ith(param,36) = 1;    /* kcat3 // 10^-12 -- 10^+6 */
	Ith(param,37) = 1;    /* Km5   // 10^-12 -- 10^+6 */
	Ith(param,38) = 1;    /* Km6   // 10^-12 -- 10^+6 */
}
