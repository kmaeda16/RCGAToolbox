#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "Threestep_mod_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

static double Function_for_vdeg1(double G1,double default0,double k1)
{
    return k1*G1/default0;
}

static double Function_for_vdeg5(double E2,double default0,double k5)
{
    return k5*E2/default0;
}

static double Function_for_vdeg3(double G3,double default0,double k3)
{
    return k3*G3/default0;
}

static double Function_for_vmet1(double E1,double Km1,double Km2,double M1,double S,double default0,double kcat1)
{
    return kcat1*E1*(1.0/Km1)*(S-M1)/(1.0+S/Km1+M1/Km2)/default0;
}

static double Function_for_vmet2(double E2,double Km3,double Km4,double M1,double M2,double default0,double kcat2)
{
    return kcat2*E2*(1.0/Km3)*(M1-M2)/(1.0+M1/Km3+M2/Km4)/default0;
}

static double Function_for_vmet3(double E3,double Km5,double Km6,double M2,double P,double default0,double kcat3)
{
    return kcat3*E3*(1.0/Km5)*(M2-P)/(1.0+M2/Km5+P/Km6)/default0;
}

static double Function_for_vtl1(double G1,double K4,double V4,double default0)
{
    return V4*G1/(K4+G1)/default0;
}

static double Function_for_vtl2(double G2,double K5,double V5,double default0)
{
    return V5*G2/(K5+G2)/default0;
}

static double Function_for_vdeg2(double G2,double default0,double k2)
{
    return k2*G2/default0;
}

static double Function_for_vtl3(double G3,double K6,double V6,double default0)
{
    return V6*G3/(K6+G3)/default0;
}

static double Function_for_vdeg4(double E1,double default0,double k4)
{
    return k4*E1/default0;
}

static double Function_for_vts1(double Ka1,double Ki1,double P,double S,double V1,double default0,double na1,double ni1)
{
    return V1/(1.0+pow(P/Ki1,ni1)+pow(Ka1/S,na1))/default0;
}

static double Function_for_vdeg6(double E3,double default0,double k6)
{
    return k6*E3/default0;
}

static double Function_for_vts3(double Ka3,double Ki3,double M2,double P,double V3,double default0,double na3,double ni3)
{
    return V3/(1.0+pow(P/Ki3,ni3)+pow(Ka3/M2,na3))/default0;
}

static double Function_for_vts2(double Ka2,double Ki2,double M1,double P,double V2,double default0,double na2,double ni2)
{
    return V2/(1.0+pow(P/Ki2,ni2)+pow(Ka2/M1,na2))/default0;
}

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double E1,E2,E3,G1,G2,G3,M1,M2,P,S,VOID;
    double K4,K5,K6,Ka1,Ka2,Ka3,Ki1,Ki2,Ki3,Km1,Km2,Km3,Km4,Km5,Km6,V1,V2,V3,V4,V5;
    double V6,k1,k2,k3,k4,k5,k6,kcat1,kcat2,kcat3,na1,na2,na3,ni1,ni2,ni3,default0;
    double vdeg1,vdeg2,vdeg3,vdeg4,vdeg5,vdeg6,vmet1,vmet2,vmet3,vtl1,vtl2,vtl3,vts1,vts2,vts3;

    time = time_local;

    E1 = stateVector[0];
    E2 = stateVector[1];
    E3 = stateVector[2];
    G1 = stateVector[3];
    G2 = stateVector[4];
    G3 = stateVector[5];
    M1 = stateVector[6];
    M2 = stateVector[7];
    P = stateVector[8];
    S = stateVector[9];
    VOID = stateVector[10];
    K4 = paramdataPtr->parametervector[0]; /* 1 */
    K5 = paramdataPtr->parametervector[1]; /* 1 */
    K6 = paramdataPtr->parametervector[2]; /* 1 */
    Ka1 = paramdataPtr->parametervector[3]; /* 1 */
    Ka2 = paramdataPtr->parametervector[4]; /* 1 */
    Ka3 = paramdataPtr->parametervector[5]; /* 1 */
    Ki1 = paramdataPtr->parametervector[6]; /* 1 */
    Ki2 = paramdataPtr->parametervector[7]; /* 1 */
    Ki3 = paramdataPtr->parametervector[8]; /* 1 */
    Km1 = paramdataPtr->parametervector[9]; /* 1 */
    Km2 = paramdataPtr->parametervector[10]; /* 1 */
    Km3 = paramdataPtr->parametervector[11]; /* 1 */
    Km4 = paramdataPtr->parametervector[12]; /* 1 */
    Km5 = paramdataPtr->parametervector[13]; /* 1 */
    Km6 = paramdataPtr->parametervector[14]; /* 1 */
    V1 = paramdataPtr->parametervector[15]; /* 1 */
    V2 = paramdataPtr->parametervector[16]; /* 1 */
    V3 = paramdataPtr->parametervector[17]; /* 1 */
    V4 = paramdataPtr->parametervector[18]; /* 0.1 */
    V5 = paramdataPtr->parametervector[19]; /* 0.1 */
    V6 = paramdataPtr->parametervector[20]; /* 0.1 */
    k1 = paramdataPtr->parametervector[21]; /* 1 */
    k2 = paramdataPtr->parametervector[22]; /* 1 */
    k3 = paramdataPtr->parametervector[23]; /* 1 */
    k4 = paramdataPtr->parametervector[24]; /* 0.1 */
    k5 = paramdataPtr->parametervector[25]; /* 0.1 */
    k6 = paramdataPtr->parametervector[26]; /* 0.1 */
    kcat1 = paramdataPtr->parametervector[27]; /* 1 */
    kcat2 = paramdataPtr->parametervector[28]; /* 1 */
    kcat3 = paramdataPtr->parametervector[29]; /* 1 */
    na1 = paramdataPtr->parametervector[30]; /* 2 */
    na2 = paramdataPtr->parametervector[31]; /* 2 */
    na3 = paramdataPtr->parametervector[32]; /* 2 */
    ni1 = paramdataPtr->parametervector[33]; /* 2 */
    ni2 = paramdataPtr->parametervector[34]; /* 2 */
    ni3 = paramdataPtr->parametervector[35]; /* 2 */
    default0 = paramdataPtr->parametervector[36]; /* 1 */
    vdeg1 = default0*Function_for_vdeg1(G1,default0,k1);
    vdeg2 = default0*Function_for_vdeg2(G2,default0,k2);
    vdeg3 = default0*Function_for_vdeg3(G3,default0,k3);
    vdeg4 = default0*Function_for_vdeg4(E1,default0,k4);
    vdeg5 = default0*Function_for_vdeg5(E2,default0,k5);
    vdeg6 = default0*Function_for_vdeg6(E3,default0,k6);
    vmet1 = default0*Function_for_vmet1(E1,Km1,Km2,M1,S,default0,kcat1);
    vmet2 = default0*Function_for_vmet2(E2,Km3,Km4,M1,M2,default0,kcat2);
    vmet3 = default0*Function_for_vmet3(E3,Km5,Km6,M2,P,default0,kcat3);
    vtl1 = default0*Function_for_vtl1(G1,K4,V4,default0);
    vtl2 = default0*Function_for_vtl2(G2,K5,V5,default0);
    vtl3 = default0*Function_for_vtl3(G3,K6,V6,default0);
    vts1 = default0*Function_for_vts1(Ka1,Ki1,P,S,V1,default0,na1,ni1);
    vts2 = default0*Function_for_vts2(Ka2,Ki2,M1,P,V2,default0,na2,ni2);
    vts3 = default0*Function_for_vts3(Ka3,Ki3,M2,P,V3,default0,na3,ni3);
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = (-vdeg4+vtl1)/default0;
    	DDTvector[1] = (-vdeg5+vtl2)/default0;
    	DDTvector[2] = (-vdeg6+vtl3)/default0;
    	DDTvector[3] = (-vdeg1+vts1)/default0;
    	DDTvector[4] = (-vdeg2+vts2)/default0;
    	DDTvector[5] = (-vdeg3+vts3)/default0;
    	DDTvector[6] = (+vmet1-vmet2)/default0;
    	DDTvector[7] = (+vmet2-vmet3)/default0;
    	DDTvector[8] = 0.0;
    	DDTvector[9] = 0.0;
    	DDTvector[10] = 0.0;
    } else if (DOflag == DOFLAG_VARREAC) {
        reactionVector[0] = vdeg1;
        reactionVector[1] = vdeg2;
        reactionVector[2] = vdeg3;
        reactionVector[3] = vdeg4;
        reactionVector[4] = vdeg5;
        reactionVector[5] = vdeg6;
        reactionVector[6] = vmet1;
        reactionVector[7] = vmet2;
        reactionVector[8] = vmet3;
        reactionVector[9] = vtl1;
        reactionVector[10] = vtl2;
        reactionVector[11] = vtl3;
        reactionVector[12] = vts1;
        reactionVector[13] = vts2;
        reactionVector[14] = vts3;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double E1,E2,E3,G1,G2,G3,M1,M2,P,S,VOID;
    double K4,K5,K6,Ka1,Ka2,Ka3,Ki1,Ki2,Ki3,Km1,Km2,Km3,Km4,Km5,Km6,V1,V2,V3,V4,V5;
    double V6,k1,k2,k3,k4,k5,k6,kcat1,kcat2,kcat3,na1,na2,na3,ni1,ni2,ni3,default0;
    K4 = paramdataPtr->parametervector[0]; /* 1 */
    K5 = paramdataPtr->parametervector[1]; /* 1 */
    K6 = paramdataPtr->parametervector[2]; /* 1 */
    Ka1 = paramdataPtr->parametervector[3]; /* 1 */
    Ka2 = paramdataPtr->parametervector[4]; /* 1 */
    Ka3 = paramdataPtr->parametervector[5]; /* 1 */
    Ki1 = paramdataPtr->parametervector[6]; /* 1 */
    Ki2 = paramdataPtr->parametervector[7]; /* 1 */
    Ki3 = paramdataPtr->parametervector[8]; /* 1 */
    Km1 = paramdataPtr->parametervector[9]; /* 1 */
    Km2 = paramdataPtr->parametervector[10]; /* 1 */
    Km3 = paramdataPtr->parametervector[11]; /* 1 */
    Km4 = paramdataPtr->parametervector[12]; /* 1 */
    Km5 = paramdataPtr->parametervector[13]; /* 1 */
    Km6 = paramdataPtr->parametervector[14]; /* 1 */
    V1 = paramdataPtr->parametervector[15]; /* 1 */
    V2 = paramdataPtr->parametervector[16]; /* 1 */
    V3 = paramdataPtr->parametervector[17]; /* 1 */
    V4 = paramdataPtr->parametervector[18]; /* 0.1 */
    V5 = paramdataPtr->parametervector[19]; /* 0.1 */
    V6 = paramdataPtr->parametervector[20]; /* 0.1 */
    k1 = paramdataPtr->parametervector[21]; /* 1 */
    k2 = paramdataPtr->parametervector[22]; /* 1 */
    k3 = paramdataPtr->parametervector[23]; /* 1 */
    k4 = paramdataPtr->parametervector[24]; /* 0.1 */
    k5 = paramdataPtr->parametervector[25]; /* 0.1 */
    k6 = paramdataPtr->parametervector[26]; /* 0.1 */
    kcat1 = paramdataPtr->parametervector[27]; /* 1 */
    kcat2 = paramdataPtr->parametervector[28]; /* 1 */
    kcat3 = paramdataPtr->parametervector[29]; /* 1 */
    na1 = paramdataPtr->parametervector[30]; /* 2 */
    na2 = paramdataPtr->parametervector[31]; /* 2 */
    na3 = paramdataPtr->parametervector[32]; /* 2 */
    ni1 = paramdataPtr->parametervector[33]; /* 2 */
    ni2 = paramdataPtr->parametervector[34]; /* 2 */
    ni3 = paramdataPtr->parametervector[35]; /* 2 */
    default0 = paramdataPtr->parametervector[36]; /* 1 */
    E1 = 0.4;
    E2 = 0.36409;
    E3 = 0.29457;
    G1 = 0.66667;
    G2 = 0.57254;
    G3 = 0.41758;
    M1 = 1.419;
    M2 = 0.93464;
    P = 0.05;
    S = 0.1;
    VOID = 0.0;
    icVector[0] = E1;
    icVector[1] = E2;
    icVector[2] = E3;
    icVector[3] = G1;
    icVector[4] = G2;
    icVector[5] = G3;
    icVector[6] = M1;
    icVector[7] = M2;
    icVector[8] = P;
    icVector[9] = S;
    icVector[10] = VOID;
}

