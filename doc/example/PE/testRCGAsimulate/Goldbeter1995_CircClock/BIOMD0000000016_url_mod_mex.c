#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "BIOMD0000000016_url_mod_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double EmptySet,M,P0,P1,P2,Pn;
    double rM_Vs,rM_KI,rM_n,rTL_ks,rP01_V1,rP01_K1,rP10_V2,rP10_K2,rP12_V3,rP12_K3,rP21_V4,rP21_K4,rP2n_k1,rPn2_k2,rmRNAd_Km,rmRNAd_Vm,rVd_Vd,rVd_Kd,default0,CYTOPLASM;
    double compartment_0000004;
    double Pt;
    double rM,rTL,rP01,rP10,rP12,rP21,rP2n,rPn2,rmRNAd,rVd;

    time = time_local;

    EmptySet = stateVector[0];
    M = stateVector[1];
    P0 = stateVector[2];
    P1 = stateVector[3];
    P2 = stateVector[4];
    Pn = stateVector[5];
    rM_Vs = paramdataPtr->parametervector[0]; /* 0.76 */
    rM_KI = paramdataPtr->parametervector[1]; /* 1 */
    rM_n = paramdataPtr->parametervector[2]; /* 4 */
    rTL_ks = paramdataPtr->parametervector[3]; /* 0.38 */
    rP01_V1 = paramdataPtr->parametervector[4]; /* 3.2 */
    rP01_K1 = paramdataPtr->parametervector[5]; /* 2 */
    rP10_V2 = paramdataPtr->parametervector[6]; /* 1.58 */
    rP10_K2 = paramdataPtr->parametervector[7]; /* 2 */
    rP12_V3 = paramdataPtr->parametervector[8]; /* 5 */
    rP12_K3 = paramdataPtr->parametervector[9]; /* 2 */
    rP21_V4 = paramdataPtr->parametervector[10]; /* 2.5 */
    rP21_K4 = paramdataPtr->parametervector[11]; /* 2 */
    rP2n_k1 = paramdataPtr->parametervector[12]; /* 1.9 */
    rPn2_k2 = paramdataPtr->parametervector[13]; /* 1.3 */
    rmRNAd_Km = paramdataPtr->parametervector[14]; /* 0.5 */
    rmRNAd_Vm = paramdataPtr->parametervector[15]; /* 0.65 */
    rVd_Vd = paramdataPtr->parametervector[16]; /* 0.95 */
    rVd_Kd = paramdataPtr->parametervector[17]; /* 0.2 */
    default0 = paramdataPtr->parametervector[18]; /* 1e-15 */
    CYTOPLASM = paramdataPtr->parametervector[19]; /* 1e-15 */
    compartment_0000004 = paramdataPtr->parametervector[20]; /* 1e-15 */
    Pt = P0+P1+P2+Pn;
    rM = default0*rM_Vs*pow(rM_KI,rM_n)/(pow(rM_KI,rM_n)+pow(Pn,rM_n));
    rTL = rTL_ks*M*default0;
    rP01 = CYTOPLASM*rP01_V1*P0/(rP01_K1+P0);
    rP10 = CYTOPLASM*rP10_V2*P1/(rP10_K2+P1);
    rP12 = CYTOPLASM*rP12_V3*P1/(rP12_K3+P1);
    rP21 = CYTOPLASM*rP21_V4*P2/(rP21_K4+P2);
    rP2n = rP2n_k1*P2*CYTOPLASM;
    rPn2 = rPn2_k2*Pn*compartment_0000004;
    rmRNAd = rmRNAd_Vm*M*CYTOPLASM/(rmRNAd_Km+M);
    rVd = CYTOPLASM*rVd_Vd*P2/(rVd_Kd+P2);
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = 0.0;
    	DDTvector[1] = (+rM-rmRNAd)/CYTOPLASM;
    	DDTvector[2] = (+rTL-rP01+rP10)/CYTOPLASM;
    	DDTvector[3] = (+rP01-rP10-rP12+rP21)/CYTOPLASM;
    	DDTvector[4] = (+rP12-rP21-rP2n+rPn2-rVd)/CYTOPLASM;
    	DDTvector[5] = (+rP2n-rPn2)/compartment_0000004;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = Pt;
        reactionVector[0] = rM;
        reactionVector[1] = rTL;
        reactionVector[2] = rP01;
        reactionVector[3] = rP10;
        reactionVector[4] = rP12;
        reactionVector[5] = rP21;
        reactionVector[6] = rP2n;
        reactionVector[7] = rPn2;
        reactionVector[8] = rmRNAd;
        reactionVector[9] = rVd;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double EmptySet,M,P0,P1,P2,Pn;
    double rM_Vs,rM_KI,rM_n,rTL_ks,rP01_V1,rP01_K1,rP10_V2,rP10_K2,rP12_V3,rP12_K3,rP21_V4,rP21_K4,rP2n_k1,rPn2_k2,rmRNAd_Km,rmRNAd_Vm,rVd_Vd,rVd_Kd,default0,CYTOPLASM;
    double compartment_0000004;
    double Pt;
    rM_Vs = paramdataPtr->parametervector[0]; /* 0.76 */
    rM_KI = paramdataPtr->parametervector[1]; /* 1 */
    rM_n = paramdataPtr->parametervector[2]; /* 4 */
    rTL_ks = paramdataPtr->parametervector[3]; /* 0.38 */
    rP01_V1 = paramdataPtr->parametervector[4]; /* 3.2 */
    rP01_K1 = paramdataPtr->parametervector[5]; /* 2 */
    rP10_V2 = paramdataPtr->parametervector[6]; /* 1.58 */
    rP10_K2 = paramdataPtr->parametervector[7]; /* 2 */
    rP12_V3 = paramdataPtr->parametervector[8]; /* 5 */
    rP12_K3 = paramdataPtr->parametervector[9]; /* 2 */
    rP21_V4 = paramdataPtr->parametervector[10]; /* 2.5 */
    rP21_K4 = paramdataPtr->parametervector[11]; /* 2 */
    rP2n_k1 = paramdataPtr->parametervector[12]; /* 1.9 */
    rPn2_k2 = paramdataPtr->parametervector[13]; /* 1.3 */
    rmRNAd_Km = paramdataPtr->parametervector[14]; /* 0.5 */
    rmRNAd_Vm = paramdataPtr->parametervector[15]; /* 0.65 */
    rVd_Vd = paramdataPtr->parametervector[16]; /* 0.95 */
    rVd_Kd = paramdataPtr->parametervector[17]; /* 0.2 */
    default0 = paramdataPtr->parametervector[18]; /* 1e-15 */
    CYTOPLASM = paramdataPtr->parametervector[19]; /* 1e-15 */
    compartment_0000004 = paramdataPtr->parametervector[20]; /* 1e-15 */
    Pt = P0+P1+P2+Pn;
    EmptySet = 0.0;
    M = 0.1;
    P0 = 0.25;
    P1 = 0.25;
    P2 = 0.25;
    Pn = 0.25;
    icVector[0] = EmptySet;
    icVector[1] = M;
    icVector[2] = P0;
    icVector[3] = P1;
    icVector[4] = P2;
    icVector[5] = Pn;
}

