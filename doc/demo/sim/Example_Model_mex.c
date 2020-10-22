#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "ExampleModel_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double X1,X2;
    double X0,k1,k2,k3,K2,K3,rootCompartment;
    double X12;
    double v1,v2,v3;

    time = time_local;

    X1 = stateVector[0];
    X2 = stateVector[1];
    X0 = paramdataPtr->parametervector[0]; /* 0.1 */
    k1 = paramdataPtr->parametervector[1]; /* 1 */
    k2 = paramdataPtr->parametervector[2]; /* 1 */
    k3 = paramdataPtr->parametervector[3]; /* 1 */
    K2 = paramdataPtr->parametervector[4]; /* 1 */
    K3 = paramdataPtr->parametervector[5]; /* 1 */
    rootCompartment = paramdataPtr->parametervector[6]; /* 1 */
    X12 = (X1/rootCompartment)+(X2/rootCompartment);
    v1 = k1*X0;
    v2 = k2*(X1/rootCompartment)/(K2+(X1/rootCompartment));
    v3 = k3*(X2/rootCompartment)/(K3+(X2/rootCompartment));
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = +v1-v2;
    	DDTvector[1] = +v2-v3;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = X12;
        reactionVector[0] = v1;
        reactionVector[1] = v2;
        reactionVector[2] = v3;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double X1,X2;
    double X0,k1,k2,k3,K2,K3,rootCompartment;
    double X12;
    X0 = paramdataPtr->parametervector[0]; /* 0.1 */
    k1 = paramdataPtr->parametervector[1]; /* 1 */
    k2 = paramdataPtr->parametervector[2]; /* 1 */
    k3 = paramdataPtr->parametervector[3]; /* 1 */
    K2 = paramdataPtr->parametervector[4]; /* 1 */
    K3 = paramdataPtr->parametervector[5]; /* 1 */
    rootCompartment = paramdataPtr->parametervector[6]; /* 1 */
    X12 = (X1/rootCompartment)+(X2/rootCompartment);
    X1 = 0.0;
    X2 = 0.0;
    icVector[0] = X1;
    icVector[1] = X2;
}

