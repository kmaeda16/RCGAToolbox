#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "HIV_mod_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

static double Function_for_vesd(double ES,double default0,double ks)
{
    return ks*ES/default0;
}

static double Function_for_veif(double E,double I,double default0,double kon)
{
    return kon*I*E/default0;
}

static double Function_for_ved(double E,double default0,double kdm)
{
    return kdm*E/default0;
}

static double Function_for_vepf(double E,double P,double default0,double kon)
{
    return kon*P*E/default0;
}

static double Function_for_vepd(double EP,double default0,double kp)
{
    return kp*EP/default0;
}

static double Function_for_vejf(double EI,double default0,double kde)
{
    return kde*EI/default0;
}

static double Function_for_vesf(double E,double S,double default0,double kon)
{
    return kon*S*E/default0;
}

static double Function_for_veid(double EI,double default0,double ki)
{
    return ki*EI/default0;
}

static double Function_for_vpf(double ES,double default0,double kcat)
{
    return kcat*ES/default0;
}

static double Function_for_vef(double M,double default0,double kmd)
{
    return kmd*M*M/default0;
}

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double E,EI,EJ,EP,ES,I,M,P,S;
    double kcat,kde,kdm,ki,kmd,kon,kp,ks,default0;
    double ved,vef,veid,veif,vejf,vepd,vepf,vesd,vesf,vpf;

    time = time_local;

    E = stateVector[0];
    EI = stateVector[1];
    EJ = stateVector[2];
    EP = stateVector[3];
    ES = stateVector[4];
    I = stateVector[5];
    M = stateVector[6];
    P = stateVector[7];
    S = stateVector[8];
    kcat = paramdataPtr->parametervector[0]; /* 5.49137 */
    kde = paramdataPtr->parametervector[1]; /* 0.000582 */
    kdm = paramdataPtr->parametervector[2]; /* 0 */
    ki = paramdataPtr->parametervector[3]; /* 0.000177 */
    kmd = paramdataPtr->parametervector[4]; /* 0.1 */
    kon = paramdataPtr->parametervector[5]; /* 100 */
    kp = paramdataPtr->parametervector[6]; /* 269.804 */
    ks = paramdataPtr->parametervector[7]; /* 46.3493 */
    default0 = paramdataPtr->parametervector[8]; /* 1 */
    ved = default0*Function_for_ved(E,default0,kdm);
    vef = default0*Function_for_vef(M,default0,kmd);
    veid = default0*Function_for_veid(EI,default0,ki);
    veif = default0*Function_for_veif(E,I,default0,kon);
    vejf = default0*Function_for_vejf(EI,default0,kde);
    vepd = default0*Function_for_vepd(EP,default0,kp);
    vepf = default0*Function_for_vepf(E,P,default0,kon);
    vesd = default0*Function_for_vesd(ES,default0,ks);
    vesf = default0*Function_for_vesf(E,S,default0,kon);
    vpf = default0*Function_for_vpf(ES,default0,kcat);
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = (-ved+vef+veid-veif+vepd-vepf+vesd-vesf+vpf)/default0;
    	DDTvector[1] = (-veid+veif-vejf)/default0;
    	DDTvector[2] = (+vejf)/default0;
    	DDTvector[3] = (-vepd+vepf)/default0;
    	DDTvector[4] = (-vesd+vesf-vpf)/default0;
    	DDTvector[5] = (+veid-veif)/default0;
    	DDTvector[6] = (+2.0*ved-2.0*vef)/default0;
    	DDTvector[7] = (+vepd-vepf+vpf)/default0;
    	DDTvector[8] = (+vesd-vesf)/default0;
    } else if (DOflag == DOFLAG_VARREAC) {
        reactionVector[0] = ved;
        reactionVector[1] = vef;
        reactionVector[2] = veid;
        reactionVector[3] = veif;
        reactionVector[4] = vejf;
        reactionVector[5] = vepd;
        reactionVector[6] = vepf;
        reactionVector[7] = vesd;
        reactionVector[8] = vesf;
        reactionVector[9] = vpf;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double E,EI,EJ,EP,ES,I,M,P,S;
    double kcat,kde,kdm,ki,kmd,kon,kp,ks,default0;
    kcat = paramdataPtr->parametervector[0]; /* 5.49137 */
    kde = paramdataPtr->parametervector[1]; /* 0.000582 */
    kdm = paramdataPtr->parametervector[2]; /* 0 */
    ki = paramdataPtr->parametervector[3]; /* 0.000177 */
    kmd = paramdataPtr->parametervector[4]; /* 0.1 */
    kon = paramdataPtr->parametervector[5]; /* 100 */
    kp = paramdataPtr->parametervector[6]; /* 269.804 */
    ks = paramdataPtr->parametervector[7]; /* 46.3493 */
    default0 = paramdataPtr->parametervector[8]; /* 1 */
    E = 0.005387;
    EI = 0.0;
    EJ = 0.0;
    EP = 0.0;
    ES = 0.0;
    I = 0.0;
    M = 0.0;
    P = 0.0;
    S = 24.6378;
    icVector[0] = E;
    icVector[1] = EI;
    icVector[2] = EJ;
    icVector[3] = EP;
    icVector[4] = ES;
    icVector[5] = I;
    icVector[6] = M;
    icVector[7] = P;
    icVector[8] = S;
}

