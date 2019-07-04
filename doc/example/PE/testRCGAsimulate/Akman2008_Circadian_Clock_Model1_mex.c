#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "Akman2008_Circadian_Clock_Model1_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double FC,FCp,FN,FNp,MF;
    double vs,ki,n,vm,km,ks,vd,k1n,k2n,ksp,vdp,k1np,k2np,amp,dawn,dusk,nucleus,cytoplasm;
    double Tot_FRQ;
    double MFtrn,MFdeg,FCtrl,FCdeg,FCtrs,FCptrl,FCpdeg,FCptrs;

    time = time_local;

    FC = stateVector[0];
    FCp = stateVector[1];
    FN = stateVector[2];
    FNp = stateVector[3];
    MF = stateVector[4];
    vs = paramdataPtr->parametervector[0]; /* 1.22363 */
    ki = paramdataPtr->parametervector[1]; /* 5.04543 */
    n = paramdataPtr->parametervector[2]; /* 6.3958 */
    vm = paramdataPtr->parametervector[3]; /* 0.885376 */
    km = paramdataPtr->parametervector[4]; /* 0.0846004 */
    ks = paramdataPtr->parametervector[5]; /* 0.313846 */
    vd = paramdataPtr->parametervector[6]; /* 0.161111 */
    k1n = paramdataPtr->parametervector[7]; /* 0.222637 */
    k2n = paramdataPtr->parametervector[8]; /* 0.331485 */
    ksp = paramdataPtr->parametervector[9]; /* 0.29484 */
    vdp = paramdataPtr->parametervector[10]; /* 0.13975 */
    k1np = paramdataPtr->parametervector[11]; /* 0.272306 */
    k2np = paramdataPtr->parametervector[12]; /* 0.295421 */
    amp = paramdataPtr->parametervector[13]; /* 0 */
    dawn = paramdataPtr->parametervector[14]; /* 6 */
    dusk = paramdataPtr->parametervector[15]; /* 18 */
    nucleus = paramdataPtr->parametervector[16]; /* 1 */
    cytoplasm = paramdataPtr->parametervector[17]; /* 1 */
    Tot_FRQ = FC+FCp+FN+FNp;
    MFtrn = (vs+amp*((1.0+tanh(2.0*(time-24*floor(time/24.0)-dawn)))*(1.0-tanh(2.0*(time-24*floor(time/24.0)-dusk)))/4.0))*pow(ki,n)/(pow(ki,n)+pow(FN+FNp,n));
    MFdeg = vm*MF/(km+MF);
    FCtrl = ks*MF;
    FCdeg = vd*FC;
    FCtrs = k1n*FC-k2n*FN;
    FCptrl = ksp*MF;
    FCpdeg = vdp*FCp;
    FCptrs = k1np*FCp-k2np*FNp;
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = (+FCtrl-FCdeg-FCtrs)/cytoplasm;
    	DDTvector[1] = (+FCptrl-FCpdeg-FCptrs)/cytoplasm;
    	DDTvector[2] = (+FCtrs)/nucleus;
    	DDTvector[3] = (+FCptrs)/nucleus;
    	DDTvector[4] = (+MFtrn-MFdeg)/nucleus;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = Tot_FRQ;
        reactionVector[0] = MFtrn;
        reactionVector[1] = MFdeg;
        reactionVector[2] = FCtrl;
        reactionVector[3] = FCdeg;
        reactionVector[4] = FCtrs;
        reactionVector[5] = FCptrl;
        reactionVector[6] = FCpdeg;
        reactionVector[7] = FCptrs;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double FC,FCp,FN,FNp,MF;
    double vs,ki,n,vm,km,ks,vd,k1n,k2n,ksp,vdp,k1np,k2np,amp,dawn,dusk,nucleus,cytoplasm;
    double Tot_FRQ;
    vs = paramdataPtr->parametervector[0]; /* 1.22363 */
    ki = paramdataPtr->parametervector[1]; /* 5.04543 */
    n = paramdataPtr->parametervector[2]; /* 6.3958 */
    vm = paramdataPtr->parametervector[3]; /* 0.885376 */
    km = paramdataPtr->parametervector[4]; /* 0.0846004 */
    ks = paramdataPtr->parametervector[5]; /* 0.313846 */
    vd = paramdataPtr->parametervector[6]; /* 0.161111 */
    k1n = paramdataPtr->parametervector[7]; /* 0.222637 */
    k2n = paramdataPtr->parametervector[8]; /* 0.331485 */
    ksp = paramdataPtr->parametervector[9]; /* 0.29484 */
    vdp = paramdataPtr->parametervector[10]; /* 0.13975 */
    k1np = paramdataPtr->parametervector[11]; /* 0.272306 */
    k2np = paramdataPtr->parametervector[12]; /* 0.295421 */
    amp = paramdataPtr->parametervector[13]; /* 0 */
    dawn = paramdataPtr->parametervector[14]; /* 6 */
    dusk = paramdataPtr->parametervector[15]; /* 18 */
    nucleus = paramdataPtr->parametervector[16]; /* 1 */
    cytoplasm = paramdataPtr->parametervector[17]; /* 1 */
    Tot_FRQ = FC+FCp+FN+FNp;
    FC = 2.46246;
    FCp = 2.71231;
    FN = 1.844;
    FNp = 2.74225;
    MF = 0.725579;
    icVector[0] = FC;
    icVector[1] = FCp;
    icVector[2] = FN;
    icVector[3] = FNp;
    icVector[4] = MF;
}

