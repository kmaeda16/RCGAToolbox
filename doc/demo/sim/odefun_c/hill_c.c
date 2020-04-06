#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "hill_c.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;


void model(double t, double *S, double *dSdt, ParamData *paramdataPtr,
	int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
	double *p;
	double S1, S2, S3;
	double S4, S0, J1_Vmax, J1_n, J1_K, J2_J2_k, J3_J3_k, J0_J0_k, compart;
	double J1, J2, J3, J0;

	time = t;

	p = paramdataPtr->parametervector;

	S1 = S[0];
	S2 = S[1];
	S3 = S[2];

	S4 = 0; // {isSpecie:compart:concentration}
	S0 = 5; // {isSpecie:compart:concentration}
	J1_Vmax = 5.5; // {isParameter}
	J1_n = 4; // {isParameter}
	J1_K = 0.5; // {isParameter}
	J2_J2_k = 0.10000000000000001; // {isParameter}
	J3_J3_k = 0.10000000000000001; // {isParameter}
	J0_J0_k = 0.01; // {isParameter}
	compart = 1; // {isCompartment:}

	S4 = p[0]; // {isSpecie:compart:concentration}
	S0 = p[1]; // {isSpecie:compart:concentration}
	J1_Vmax = p[2]; // {isParameter}
	J1_n = p[3]; // {isParameter}
	J1_K = p[4]; // {isParameter}
	J2_J2_k = p[5]; // {isParameter}
	J3_J3_k = p[6]; // {isParameter}
	J0_J0_k = p[7]; // {isParameter}
	compart = p[8]; // {isCompartment:}

	J1 = J1_Vmax * pow(S1,J1_n) / (pow(J1_K,J1_n) + pow(S1,J1_n));
	J2 = J2_J2_k * S2;
	J3 = J3_J3_k * S3;
	J0 = J0_J0_k * S0;

	if(DOflag == DOFLAG_DDT){
    	dSdt[0] = (-J1 + J0) / compart; // {isSpecie:compart:concentration}
		dSdt[1] = (+J1 - J2) / compart; // {isSpecie:compart:concentration}
		dSdt[2] = (+J2 - J3) / compart; // {isSpecie:compart:concentration}
    }else if (DOflag == DOFLAG_VARREAC){
        reactionVector[0] = J1;
        reactionVector[1] = J2;
        reactionVector[2] = J3;
        reactionVector[3] = J0;
    }else if(DOflag == DOFLAG_EVENTS){
	}else if(DOflag == DOFLAG_EVENTASSIGN){
    }

}

void calc_ic_model(double *S_Init, ParamData *paramdataPtr)
{
	double S1, S2, S3;

	S1 = 0.0;
    S2 = 0.0;
    S3 = 0.0;

	S_Init[0] = S1;
	S_Init[1] = S2;
	S_Init[2] = S3;
}