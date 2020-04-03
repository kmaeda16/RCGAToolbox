#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "SBMLexampleLevel2_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double S1,S2,S3;
    double S4,S0,J1_Vmax,J1_n,J1_K,J2_J2_k,J3_J3_k,J0_J0_k,compart;
    double J1,J2,J3,J0;

    time = time_local;

    S1 = stateVector[0];
    S2 = stateVector[1];
    S3 = stateVector[2];
    S4 = paramdataPtr->parametervector[0]; /* 0 */
    S0 = paramdataPtr->parametervector[1]; /* 5 */
    J1_Vmax = paramdataPtr->parametervector[2]; /* 5.5 */
    J1_n = paramdataPtr->parametervector[3]; /* 4 */
    J1_K = paramdataPtr->parametervector[4]; /* 0.5 */
    J2_J2_k = paramdataPtr->parametervector[5]; /* 0.1 */
    J3_J3_k = paramdataPtr->parametervector[6]; /* 0.1 */
    J0_J0_k = paramdataPtr->parametervector[7]; /* 0.01 */
    compart = paramdataPtr->parametervector[8]; /* 1 */
    J1 = J1_Vmax*pow(S1,J1_n)/(pow(J1_K,J1_n)+pow(S1,J1_n));
    J2 = J2_J2_k*S2;
    J3 = J3_J3_k*S3;
    J0 = J0_J0_k*S0;
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = (-J1+J0)/compart;
    	DDTvector[1] = (+J1-J2)/compart;
    	DDTvector[2] = (+J2-J3)/compart;
    } else if (DOflag == DOFLAG_VARREAC) {
        reactionVector[0] = J1;
        reactionVector[1] = J2;
        reactionVector[2] = J3;
        reactionVector[3] = J0;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double S1,S2,S3;
    double S4,S0,J1_Vmax,J1_n,J1_K,J2_J2_k,J3_J3_k,J0_J0_k,compart;
    S4 = paramdataPtr->parametervector[0]; /* 0 */
    S0 = paramdataPtr->parametervector[1]; /* 5 */
    J1_Vmax = paramdataPtr->parametervector[2]; /* 5.5 */
    J1_n = paramdataPtr->parametervector[3]; /* 4 */
    J1_K = paramdataPtr->parametervector[4]; /* 0.5 */
    J2_J2_k = paramdataPtr->parametervector[5]; /* 0.1 */
    J3_J3_k = paramdataPtr->parametervector[6]; /* 0.1 */
    J0_J0_k = paramdataPtr->parametervector[7]; /* 0.01 */
    compart = paramdataPtr->parametervector[8]; /* 1 */
    S1 = 0.0;
    S2 = 0.0;
    S3 = 0.0;
    icVector[0] = S1;
    icVector[1] = S2;
    icVector[2] = S3;
}

