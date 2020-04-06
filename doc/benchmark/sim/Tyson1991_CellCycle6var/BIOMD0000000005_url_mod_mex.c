#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "BIOMD0000000005_url_mod_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double C2,CP,M,pM,Y,YP;
    double EmptySet,Reaction1_k6,Reaction2_k8notP,Reaction3_k9,Reaction4_k3,Reaction5_k5notP,Reaction6_k1aa,Reaction7_k2,Reaction8_k7,Reaction9_k4,Reaction9_k4prime,cell;
    double YT,CT;
    double Reaction1,Reaction2,Reaction3,Reaction4,Reaction5,Reaction6,Reaction7,Reaction8,Reaction9;

    time = time_local;

    C2 = stateVector[0];
    CP = stateVector[1];
    M = stateVector[2];
    pM = stateVector[3];
    Y = stateVector[4];
    YP = stateVector[5];
    EmptySet = paramdataPtr->parametervector[0]; /* 0 */
    Reaction1_k6 = paramdataPtr->parametervector[1]; /* 1 */
    Reaction2_k8notP = paramdataPtr->parametervector[2]; /* 1e+06 */
    Reaction3_k9 = paramdataPtr->parametervector[3]; /* 1000 */
    Reaction4_k3 = paramdataPtr->parametervector[4]; /* 200 */
    Reaction5_k5notP = paramdataPtr->parametervector[5]; /* 0 */
    Reaction6_k1aa = paramdataPtr->parametervector[6]; /* 0.015 */
    Reaction7_k2 = paramdataPtr->parametervector[7]; /* 0 */
    Reaction8_k7 = paramdataPtr->parametervector[8]; /* 0.6 */
    Reaction9_k4 = paramdataPtr->parametervector[9]; /* 180 */
    Reaction9_k4prime = paramdataPtr->parametervector[10]; /* 0.018 */
    cell = paramdataPtr->parametervector[11]; /* 1 */
    YT = (Y/cell)+(YP/cell)+(M/cell)+(pM/cell);
    CT = (C2/cell)+(CP/cell)+(M/cell)+(pM/cell);
    Reaction1 = cell*Reaction1_k6*(M/cell);
    Reaction2 = cell*(C2/cell)*Reaction2_k8notP;
    Reaction3 = cell*(CP/cell)*Reaction3_k9;
    Reaction4 = cell*(CP/cell)*Reaction4_k3*(Y/cell);
    Reaction5 = cell*Reaction5_k5notP*(M/cell);
    Reaction6 = cell*Reaction6_k1aa;
    Reaction7 = cell*Reaction7_k2*(Y/cell);
    Reaction8 = cell*Reaction8_k7*(YP/cell);
    Reaction9 = cell*(pM/cell)*(Reaction9_k4prime+Reaction9_k4*pow((M/cell)/CT,2.0));
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = +Reaction1-Reaction2+Reaction3;
    	DDTvector[1] = +Reaction2-Reaction3-Reaction4;
    	DDTvector[2] = -Reaction1-Reaction5+Reaction9;
    	DDTvector[3] = +Reaction4+Reaction5-Reaction9;
    	DDTvector[4] = -Reaction4+Reaction6-Reaction7;
    	DDTvector[5] = +Reaction1-Reaction8;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = YT;
        variableVector[1] = CT;
        reactionVector[0] = Reaction1;
        reactionVector[1] = Reaction2;
        reactionVector[2] = Reaction3;
        reactionVector[3] = Reaction4;
        reactionVector[4] = Reaction5;
        reactionVector[5] = Reaction6;
        reactionVector[6] = Reaction7;
        reactionVector[7] = Reaction8;
        reactionVector[8] = Reaction9;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double C2,CP,M,pM,Y,YP;
    double EmptySet,Reaction1_k6,Reaction2_k8notP,Reaction3_k9,Reaction4_k3,Reaction5_k5notP,Reaction6_k1aa,Reaction7_k2,Reaction8_k7,Reaction9_k4,Reaction9_k4prime,cell;
    double YT,CT;
    EmptySet = paramdataPtr->parametervector[0]; /* 0 */
    Reaction1_k6 = paramdataPtr->parametervector[1]; /* 1 */
    Reaction2_k8notP = paramdataPtr->parametervector[2]; /* 1e+06 */
    Reaction3_k9 = paramdataPtr->parametervector[3]; /* 1000 */
    Reaction4_k3 = paramdataPtr->parametervector[4]; /* 200 */
    Reaction5_k5notP = paramdataPtr->parametervector[5]; /* 0 */
    Reaction6_k1aa = paramdataPtr->parametervector[6]; /* 0.015 */
    Reaction7_k2 = paramdataPtr->parametervector[7]; /* 0 */
    Reaction8_k7 = paramdataPtr->parametervector[8]; /* 0.6 */
    Reaction9_k4 = paramdataPtr->parametervector[9]; /* 180 */
    Reaction9_k4prime = paramdataPtr->parametervector[10]; /* 0.018 */
    cell = paramdataPtr->parametervector[11]; /* 1 */
    YT = (Y/cell)+(YP/cell)+(M/cell)+(pM/cell);
    CT = (C2/cell)+(CP/cell)+(M/cell)+(pM/cell);
    C2 = 0.0;
    CP = 0.75;
    M = 0.0;
    pM = 0.25;
    Y = 0.0;
    YP = 0.0;
    icVector[0] = C2;
    icVector[1] = CP;
    icVector[2] = M;
    icVector[3] = pM;
    icVector[4] = Y;
    icVector[5] = YP;
}

