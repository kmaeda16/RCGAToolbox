#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "Bruggeman2005_AmmoniumAssimilation_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double PII,UTP,PIIUMP,PPi,GLN,PIIUMP2,PIIUMP3,UMP,GS,AMP,NH4,KG,NADPH,GLU,NADP,AZGLU,ATP,ADP,AZglu,AZGLN;
    double AZgln;
    double P_i,UT,kcatut,Kglnut,Kutipii,Kutpii,Kutpiiump,Kututp,Kutippi,UR,kcatur,Kurpiiump,Kurump,Kglnur,a1,b1,c1,d1,Vad,Kadpiikg;
    double Kadgln,Kadgs,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,Vdead,Kdeadpiikg,Kdeadgln,Kdeadpiiu,Kdeadgsa,Vgdh,Kgdhkg;
    double Kgdhnh,Kgdhglu,Kgdhnadph,Kgdhnadp,Keqgdh,Kgdhazglu,Vgog,Kgoggln,Kgogkg,Kgognadph,Kgogglu,Kgognadp,Kgogaz,Vgs,aamp,bamp,camp,damp,n1amp,n2amp;
    double Kgseq,Kgsatp,Kgsglu,Kgsnh,Kgsadp,Kgspi,Kgsgln,Keq,Vgludem,Kgludemglu,Kgludemeq,Kgludemazglu,Vglndem,Kglndemgln,Kglndemeq,Kglndemazgln,Vazglndem,Kazglndemazgln,Kazglndemeq,Kazglndemazinter;
    double Vazgludem,Kazgludemazglu,Kazgludemeq,Kazgludemazinter,Vadp,Kadp,ATPtot,GStot,PIItot,Kd1,Kd2,Kd3,Kd1piiump,Kd2piiump,Kd3piiump,compartment;
    double vAPP_GS,nAMP,PIIKG1,PIIUMP3KG3;
    double vut1,vur1,vut2,vur2,vut3,vur3,vad,vdead,vgdh,vgog,vgs,vgludem,vazgludem,vglndem,vazglndem,vatpase;

    time = time_local;

    PII = stateVector[0];
    UTP = stateVector[1];
    PIIUMP = stateVector[2];
    PPi = stateVector[3];
    GLN = stateVector[4];
    PIIUMP2 = stateVector[5];
    PIIUMP3 = stateVector[6];
    UMP = stateVector[7];
    GS = stateVector[8];
    AMP = stateVector[9];
    NH4 = stateVector[10];
    KG = stateVector[11];
    NADPH = stateVector[12];
    GLU = stateVector[13];
    NADP = stateVector[14];
    AZGLU = stateVector[15];
    ATP = stateVector[16];
    ADP = stateVector[17];
    AZglu = stateVector[18];
    AZGLN = stateVector[19];
    AZgln = stateVector[20];
    P_i = paramdataPtr->parametervector[0]; /* 10 */
    UT = paramdataPtr->parametervector[1]; /* 0.0006 */
    kcatut = paramdataPtr->parametervector[2]; /* 137 */
    Kglnut = paramdataPtr->parametervector[3]; /* 0.07 */
    Kutipii = paramdataPtr->parametervector[4]; /* 0.0018 */
    Kutpii = paramdataPtr->parametervector[5]; /* 0.003 */
    Kutpiiump = paramdataPtr->parametervector[6]; /* 0.0035 */
    Kututp = paramdataPtr->parametervector[7]; /* 0.04 */
    Kutippi = paramdataPtr->parametervector[8]; /* 0.1135 */
    UR = paramdataPtr->parametervector[9]; /* 0.0006 */
    kcatur = paramdataPtr->parametervector[10]; /* 5.5 */
    Kurpiiump = paramdataPtr->parametervector[11]; /* 0.0023 */
    Kurump = paramdataPtr->parametervector[12]; /* 8.4 */
    Kglnur = paramdataPtr->parametervector[13]; /* 0.07 */
    a1 = paramdataPtr->parametervector[14]; /* 1e-22 */
    b1 = paramdataPtr->parametervector[15]; /* 0.5166 */
    c1 = paramdataPtr->parametervector[16]; /* 0.5974 */
    d1 = paramdataPtr->parametervector[17]; /* 0.0387 */
    Vad = paramdataPtr->parametervector[18]; /* 0.5 */
    Kadpiikg = paramdataPtr->parametervector[19]; /* 1.052e-05 */
    Kadgln = paramdataPtr->parametervector[20]; /* 0.9714 */
    Kadgs = paramdataPtr->parametervector[21]; /* 0.001703 */
    e1 = paramdataPtr->parametervector[22]; /* 1e-22 */
    f1 = paramdataPtr->parametervector[23]; /* 2.766 */
    g1 = paramdataPtr->parametervector[24]; /* 3.323 */
    h1 = paramdataPtr->parametervector[25]; /* 0.2148 */
    i1 = paramdataPtr->parametervector[26]; /* 1e-22 */
    j1 = paramdataPtr->parametervector[27]; /* 1e-22 */
    k1 = paramdataPtr->parametervector[28]; /* 1e-22 */
    l1 = paramdataPtr->parametervector[29]; /* 0.02316 */
    m1 = paramdataPtr->parametervector[30]; /* 0.8821 */
    n1 = paramdataPtr->parametervector[31]; /* 8.491 */
    o1 = paramdataPtr->parametervector[32]; /* 0.8791 */
    Vdead = paramdataPtr->parametervector[33]; /* 0.5 */
    Kdeadpiikg = paramdataPtr->parametervector[34]; /* 2.274e-06 */
    Kdeadgln = paramdataPtr->parametervector[35]; /* 0.04444 */
    Kdeadpiiu = paramdataPtr->parametervector[36]; /* 1.805e-05 */
    Kdeadgsa = paramdataPtr->parametervector[37]; /* 0.0002015 */
    Vgdh = paramdataPtr->parametervector[38]; /* 360 */
    Kgdhkg = paramdataPtr->parametervector[39]; /* 0.32 */
    Kgdhnh = paramdataPtr->parametervector[40]; /* 1.1 */
    Kgdhglu = paramdataPtr->parametervector[41]; /* 10 */
    Kgdhnadph = paramdataPtr->parametervector[42]; /* 0.04 */
    Kgdhnadp = paramdataPtr->parametervector[43]; /* 0.042 */
    Keqgdh = paramdataPtr->parametervector[44]; /* 1290 */
    Kgdhazglu = paramdataPtr->parametervector[45]; /* 2.5 */
    Vgog = paramdataPtr->parametervector[46]; /* 85 */
    Kgoggln = paramdataPtr->parametervector[47]; /* 0.175 */
    Kgogkg = paramdataPtr->parametervector[48]; /* 0.007 */
    Kgognadph = paramdataPtr->parametervector[49]; /* 0.0015 */
    Kgogglu = paramdataPtr->parametervector[50]; /* 11 */
    Kgognadp = paramdataPtr->parametervector[51]; /* 0.0037 */
    Kgogaz = paramdataPtr->parametervector[52]; /* 0.65 */
    Vgs = paramdataPtr->parametervector[53]; /* 600 */
    aamp = paramdataPtr->parametervector[54]; /* 10 */
    bamp = paramdataPtr->parametervector[55]; /* 2.3667 */
    camp = paramdataPtr->parametervector[56]; /* 0.1012 */
    damp = paramdataPtr->parametervector[57]; /* 10.8688 */
    n1amp = paramdataPtr->parametervector[58]; /* 1.1456 */
    n2amp = paramdataPtr->parametervector[59]; /* 19.2166 */
    Kgseq = paramdataPtr->parametervector[60]; /* 460 */
    Kgsatp = paramdataPtr->parametervector[61]; /* 0.35 */
    Kgsglu = paramdataPtr->parametervector[62]; /* 4.1 */
    Kgsnh = paramdataPtr->parametervector[63]; /* 0.1 */
    Kgsadp = paramdataPtr->parametervector[64]; /* 0.0585 */
    Kgspi = paramdataPtr->parametervector[65]; /* 3.7 */
    Kgsgln = paramdataPtr->parametervector[66]; /* 5.65 */
    Keq = paramdataPtr->parametervector[67]; /* 460 */
    Vgludem = paramdataPtr->parametervector[68]; /* 120 */
    Kgludemglu = paramdataPtr->parametervector[69]; /* 8 */
    Kgludemeq = paramdataPtr->parametervector[70]; /* 1e+10 */
    Kgludemazglu = paramdataPtr->parametervector[71]; /* 0.5 */
    Vglndem = paramdataPtr->parametervector[72]; /* 70 */
    Kglndemgln = paramdataPtr->parametervector[73]; /* 2 */
    Kglndemeq = paramdataPtr->parametervector[74]; /* 1e+10 */
    Kglndemazgln = paramdataPtr->parametervector[75]; /* 0.25 */
    Vazglndem = paramdataPtr->parametervector[76]; /* 20 */
    Kazglndemazgln = paramdataPtr->parametervector[77]; /* 1 */
    Kazglndemeq = paramdataPtr->parametervector[78]; /* 1e+10 */
    Kazglndemazinter = paramdataPtr->parametervector[79]; /* 0.5 */
    Vazgludem = paramdataPtr->parametervector[80]; /* 30 */
    Kazgludemazglu = paramdataPtr->parametervector[81]; /* 0.3 */
    Kazgludemeq = paramdataPtr->parametervector[82]; /* 1e+10 */
    Kazgludemazinter = paramdataPtr->parametervector[83]; /* 0.5 */
    Vadp = paramdataPtr->parametervector[84]; /* 100 */
    Kadp = paramdataPtr->parametervector[85]; /* 0.5 */
    ATPtot = paramdataPtr->parametervector[86]; /* 5.37 */
    GStot = paramdataPtr->parametervector[87]; /* 0.014 */
    PIItot = paramdataPtr->parametervector[88]; /* 0.003 */
    Kd1 = paramdataPtr->parametervector[89]; /* 0.005 */
    Kd2 = paramdataPtr->parametervector[90]; /* 0.15 */
    Kd3 = paramdataPtr->parametervector[91]; /* 0.15 */
    Kd1piiump = paramdataPtr->parametervector[92]; /* 0.025 */
    Kd2piiump = paramdataPtr->parametervector[93]; /* 0.15 */
    Kd3piiump = paramdataPtr->parametervector[94]; /* 0.15 */
    compartment = paramdataPtr->parametervector[95]; /* 1 */
    vAPP_GS = aamp*camp/((1.0+pow(12.0,n1amp)*pow(AMP/(bamp*GStot),n1amp))*(1.0+pow(12.0,n2amp)*pow(AMP/(damp*GStot),n2amp)))*Vgs;
    nAMP = 12.0*AMP/GStot;
    PIIKG1 = 3.0*PII*KG/Kd1/(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3));
    PIIUMP3KG3 = PIIUMP3*pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)/(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump));
    vut1 = compartment*(kcatut*UT*UTP*PII/(Kutipii*Kututp*(1.0+GLN/Kglnut)*(1.0+UTP/Kututp+(PII+PIIUMP+PIIUMP2)/Kutipii+UTP*(PII+PIIUMP+PIIUMP2)/(Kutipii*Kututp)+PPi*UTP*(PII+PIIUMP+PIIUMP2)/(Kutipii*Kutippi*Kututp)+Kutpii*(PIIUMP+PIIUMP2+PIIUMP3)/(Kutipii*Kutpiiump))));
    vur1 = compartment*(kcatur*UR*PIIUMP/(Kurpiiump*(1.0+Kglnur/GLN)*(1.0+(1.0+UMP/Kurump)*(PIIUMP+PIIUMP2+PIIUMP3)/Kurpiiump)));
    vut2 = compartment*(kcatut*UT*UTP*PIIUMP/(Kutipii*Kututp*(1.0+GLN/Kglnut)*(1.0+UTP/Kututp+(PII+PIIUMP+PIIUMP2)/Kutipii+UTP*(PII+PIIUMP+PIIUMP2)/(Kutipii*Kututp)+PPi*UTP*(PII+PIIUMP+PIIUMP2)/(Kutipii*Kutippi*Kututp)+Kutpii*(PIIUMP+PIIUMP2+PIIUMP3)/(Kutipii*Kutpiiump))));
    vur2 = compartment*(kcatur*UR*PIIUMP2/(Kurpiiump*(1.0+Kglnur/GLN)*(1.0+(1.0+UMP/Kurump)*(PIIUMP+PIIUMP2+PIIUMP3)/Kurpiiump)));
    vut3 = compartment*(kcatut*UT*UTP*PIIUMP2/(Kutipii*Kututp*(1.0+GLN/Kglnut)*(1.0+UTP/Kututp+(PII+PIIUMP+PIIUMP2)/Kutipii+UTP*(PII+PIIUMP+PIIUMP2)/(Kutipii*Kututp)+PPi*UTP*(PII+PIIUMP+PIIUMP2)/(Kutipii*Kutippi*Kututp)+Kutpii*(PIIUMP+PIIUMP2+PIIUMP3)/(Kutipii*Kutpiiump))));
    vur3 = compartment*(kcatur*UR*PIIUMP3/(Kurpiiump*(1.0+Kglnur/GLN)*(1.0+(1.0+UMP/Kurump)*(PIIUMP+PIIUMP2+PIIUMP3)/Kurpiiump)));
    vad = compartment*(Vad*GS*(b1*GLN/Kadgln+3.0*a1*KG*PII/(Kadpiikg*Kd1*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3)))+3.0*c1*KG*GLN*PII/(Kadgln*Kadpiikg*Kd1*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))))/((Kadgs+GS)*(1.0+GLN/Kadgln+3.0*KG*PII/(Kadpiikg*Kd1*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3)))+3.0*KG*GLN*PII/(d1*Kadgln*Kadpiikg*Kd1*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))))));
    vdead = compartment*(Vdead*AMP*(f1*GLN/Kdeadgln+3.0*e1*KG*PII/(Kd1*Kdeadpiikg*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3)))+3.0*h1*KG*GLN*PII/(Kd1*Kdeadgln*Kdeadpiikg*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3)))+g1*pow(KG,3.0)*PIIUMP3/(Kd1piiump*Kd2piiump*Kd3piiump*Kdeadpiiu*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)))+j1*pow(KG,3.0)*GLN*PIIUMP3/(Kd1piiump*Kd2piiump*Kd3piiump*Kdeadgln*Kdeadpiiu*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)))+3.0*i1*pow(KG,4.0)*PII*PIIUMP3/(Kd1*Kd1piiump*Kd2piiump*Kd3piiump*Kdeadpiikg*Kdeadpiiu*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)))+3.0*k1*pow(KG,4.0)*GLN*PII*PIIUMP3/(Kd1*Kd1piiump*Kd2piiump*Kd3piiump*Kdeadgln*Kdeadpiikg*Kdeadpiiu*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump))))/((Kdeadgsa+AMP)*(1.0+GLN/Kdeadgln+3.0*KG*PII/(Kd1*Kdeadpiikg*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3)))+3.0*KG*GLN*PII/(Kd1*Kdeadgln*Kdeadpiikg*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))*l1)+pow(KG,3.0)*PIIUMP3/(Kd1piiump*Kd2piiump*Kd3piiump*Kdeadpiiu*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)))+pow(KG,3.0)*GLN*PIIUMP3/(Kd1piiump*Kd2piiump*Kd3piiump*Kdeadgln*Kdeadpiiu*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump))*n1)+3.0*pow(KG,4.0)*PII*PIIUMP3/(Kd1*Kd1piiump*Kd2piiump*Kd3piiump*Kdeadpiikg*Kdeadpiiu*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump))*m1)+3.0*pow(KG,4.0)*GLN*PII*PIIUMP3/(Kd1*Kd1piiump*Kd2piiump*Kd3piiump*Kdeadgln*Kdeadpiikg*Kdeadpiiu*(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3))*(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump))*o1))));
    vgdh = compartment*(Vgdh*(KG*NADPH*NH4-NADP*GLU/Keqgdh)/(Kgdhkg*Kgdhnadph*Kgdhnh*(1.0+NADP/Kgdhnadp+NADPH/Kgdhnadph)*(1.0+NH4/Kgdhnh)*(1.0+KG/Kgdhkg+GLU/Kgdhglu)));
    vgog = compartment*(KG*NADPH*Vgog*GLN/(Kgoggln*Kgogkg*Kgognadph*(1.0+NADP/Kgognadp+NADPH/Kgognadph)*(1.0+AZGLU/Kgogaz)*(1.0+KG/Kgogkg+GLU/Kgogglu)*(1.0+GLN/Kgoggln+GLU/Kgogglu)));
    vgs = compartment*(aamp*camp*Vgs*(-(P_i*ADP*GLN/Keq)+NH4*ATP*GLU)/(Kgsatp*Kgsglu*Kgsnh*(1.0+P_i/Kgspi+ADP/Kgsadp+P_i*ADP/(Kgsadp*Kgspi)+ATP/Kgsatp)*(1.0+NH4/Kgsnh+GLN/Kgsgln+NH4*GLN/(Kgsgln*Kgsnh)+GLU/Kgsglu+NH4*GLU/(Kgsglu*Kgsnh))*(1.0+pow(12.0,n1amp)*pow(AMP/(bamp*GStot),n1amp))*(1.0+pow(12.0,n2amp)*pow(AMP/(damp*GStot),n2amp))));
    vgludem = compartment*(Vgludem*(-(AZGLU/Kgludemeq)+GLU)/(Kgludemglu*(1.0+AZGLU/Kgludemazglu+GLU/Kgludemglu)));
    vazgludem = compartment*(Vazgludem*(-(AZglu/Kazgludemeq)+AZGLU)/(Kazgludemazglu*(1.0+AZglu/Kazgludemazinter+AZGLU/Kazgludemazglu)));
    vglndem = compartment*(Vglndem*(-(AZGLN/Kglndemeq)+GLN)/(Kglndemgln*(1.0+AZGLN/Kglndemazgln+GLN/Kglndemgln)));
    vazglndem = compartment*(Vazglndem*(-(AZgln/Kazglndemeq)+AZGLN)/(Kazglndemazgln*(1.0+AZgln/Kazglndemazinter+AZGLN/Kazglndemazgln)));
    vatpase = compartment*(Vadp*ADP/(Kadp+ADP));
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = (-vut1+vur1)/compartment;
    	DDTvector[1] = 0.0;
    	DDTvector[2] = (+vut1-vur1-vut2+vur2)/compartment;
    	DDTvector[3] = 0.0;
    	DDTvector[4] = (-vgog+vgs-vglndem)/compartment;
    	DDTvector[5] = (+vut2-vur2-vut3+vur3)/compartment;
    	DDTvector[6] = (+vut3-vur3)/compartment;
    	DDTvector[7] = 0.0;
    	DDTvector[8] = (-vad+vdead)/compartment;
    	DDTvector[9] = (+vad-vdead)/compartment;
    	DDTvector[10] = 0.0;
    	DDTvector[11] = 0.0;
    	DDTvector[12] = 0.0;
    	DDTvector[13] = (+vgdh+2.0*vgog-vgs-vgludem)/compartment;
    	DDTvector[14] = 0.0;
    	DDTvector[15] = (+vgludem-vazgludem)/compartment;
    	DDTvector[16] = (-vgs+vatpase)/compartment;
    	DDTvector[17] = (+vgs-vatpase)/compartment;
    	DDTvector[18] = 0.0;
    	DDTvector[19] = (+vglndem-vazglndem)/compartment;
    	DDTvector[20] = 0.0;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = vAPP_GS;
        variableVector[1] = nAMP;
        variableVector[2] = PIIKG1;
        variableVector[3] = PIIUMP3KG3;
        reactionVector[0] = vut1;
        reactionVector[1] = vur1;
        reactionVector[2] = vut2;
        reactionVector[3] = vur2;
        reactionVector[4] = vut3;
        reactionVector[5] = vur3;
        reactionVector[6] = vad;
        reactionVector[7] = vdead;
        reactionVector[8] = vgdh;
        reactionVector[9] = vgog;
        reactionVector[10] = vgs;
        reactionVector[11] = vgludem;
        reactionVector[12] = vazgludem;
        reactionVector[13] = vglndem;
        reactionVector[14] = vazglndem;
        reactionVector[15] = vatpase;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double PII,UTP,PIIUMP,PPi,GLN,PIIUMP2,PIIUMP3,UMP,GS,AMP,NH4,KG,NADPH,GLU,NADP,AZGLU,ATP,ADP,AZglu,AZGLN;
    double AZgln;
    double P_i,UT,kcatut,Kglnut,Kutipii,Kutpii,Kutpiiump,Kututp,Kutippi,UR,kcatur,Kurpiiump,Kurump,Kglnur,a1,b1,c1,d1,Vad,Kadpiikg;
    double Kadgln,Kadgs,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,Vdead,Kdeadpiikg,Kdeadgln,Kdeadpiiu,Kdeadgsa,Vgdh,Kgdhkg;
    double Kgdhnh,Kgdhglu,Kgdhnadph,Kgdhnadp,Keqgdh,Kgdhazglu,Vgog,Kgoggln,Kgogkg,Kgognadph,Kgogglu,Kgognadp,Kgogaz,Vgs,aamp,bamp,camp,damp,n1amp,n2amp;
    double Kgseq,Kgsatp,Kgsglu,Kgsnh,Kgsadp,Kgspi,Kgsgln,Keq,Vgludem,Kgludemglu,Kgludemeq,Kgludemazglu,Vglndem,Kglndemgln,Kglndemeq,Kglndemazgln,Vazglndem,Kazglndemazgln,Kazglndemeq,Kazglndemazinter;
    double Vazgludem,Kazgludemazglu,Kazgludemeq,Kazgludemazinter,Vadp,Kadp,ATPtot,GStot,PIItot,Kd1,Kd2,Kd3,Kd1piiump,Kd2piiump,Kd3piiump,compartment;
    double vAPP_GS,nAMP,PIIKG1,PIIUMP3KG3;
    P_i = paramdataPtr->parametervector[0]; /* 10 */
    UT = paramdataPtr->parametervector[1]; /* 0.0006 */
    kcatut = paramdataPtr->parametervector[2]; /* 137 */
    Kglnut = paramdataPtr->parametervector[3]; /* 0.07 */
    Kutipii = paramdataPtr->parametervector[4]; /* 0.0018 */
    Kutpii = paramdataPtr->parametervector[5]; /* 0.003 */
    Kutpiiump = paramdataPtr->parametervector[6]; /* 0.0035 */
    Kututp = paramdataPtr->parametervector[7]; /* 0.04 */
    Kutippi = paramdataPtr->parametervector[8]; /* 0.1135 */
    UR = paramdataPtr->parametervector[9]; /* 0.0006 */
    kcatur = paramdataPtr->parametervector[10]; /* 5.5 */
    Kurpiiump = paramdataPtr->parametervector[11]; /* 0.0023 */
    Kurump = paramdataPtr->parametervector[12]; /* 8.4 */
    Kglnur = paramdataPtr->parametervector[13]; /* 0.07 */
    a1 = paramdataPtr->parametervector[14]; /* 1e-22 */
    b1 = paramdataPtr->parametervector[15]; /* 0.5166 */
    c1 = paramdataPtr->parametervector[16]; /* 0.5974 */
    d1 = paramdataPtr->parametervector[17]; /* 0.0387 */
    Vad = paramdataPtr->parametervector[18]; /* 0.5 */
    Kadpiikg = paramdataPtr->parametervector[19]; /* 1.052e-05 */
    Kadgln = paramdataPtr->parametervector[20]; /* 0.9714 */
    Kadgs = paramdataPtr->parametervector[21]; /* 0.001703 */
    e1 = paramdataPtr->parametervector[22]; /* 1e-22 */
    f1 = paramdataPtr->parametervector[23]; /* 2.766 */
    g1 = paramdataPtr->parametervector[24]; /* 3.323 */
    h1 = paramdataPtr->parametervector[25]; /* 0.2148 */
    i1 = paramdataPtr->parametervector[26]; /* 1e-22 */
    j1 = paramdataPtr->parametervector[27]; /* 1e-22 */
    k1 = paramdataPtr->parametervector[28]; /* 1e-22 */
    l1 = paramdataPtr->parametervector[29]; /* 0.02316 */
    m1 = paramdataPtr->parametervector[30]; /* 0.8821 */
    n1 = paramdataPtr->parametervector[31]; /* 8.491 */
    o1 = paramdataPtr->parametervector[32]; /* 0.8791 */
    Vdead = paramdataPtr->parametervector[33]; /* 0.5 */
    Kdeadpiikg = paramdataPtr->parametervector[34]; /* 2.274e-06 */
    Kdeadgln = paramdataPtr->parametervector[35]; /* 0.04444 */
    Kdeadpiiu = paramdataPtr->parametervector[36]; /* 1.805e-05 */
    Kdeadgsa = paramdataPtr->parametervector[37]; /* 0.0002015 */
    Vgdh = paramdataPtr->parametervector[38]; /* 360 */
    Kgdhkg = paramdataPtr->parametervector[39]; /* 0.32 */
    Kgdhnh = paramdataPtr->parametervector[40]; /* 1.1 */
    Kgdhglu = paramdataPtr->parametervector[41]; /* 10 */
    Kgdhnadph = paramdataPtr->parametervector[42]; /* 0.04 */
    Kgdhnadp = paramdataPtr->parametervector[43]; /* 0.042 */
    Keqgdh = paramdataPtr->parametervector[44]; /* 1290 */
    Kgdhazglu = paramdataPtr->parametervector[45]; /* 2.5 */
    Vgog = paramdataPtr->parametervector[46]; /* 85 */
    Kgoggln = paramdataPtr->parametervector[47]; /* 0.175 */
    Kgogkg = paramdataPtr->parametervector[48]; /* 0.007 */
    Kgognadph = paramdataPtr->parametervector[49]; /* 0.0015 */
    Kgogglu = paramdataPtr->parametervector[50]; /* 11 */
    Kgognadp = paramdataPtr->parametervector[51]; /* 0.0037 */
    Kgogaz = paramdataPtr->parametervector[52]; /* 0.65 */
    Vgs = paramdataPtr->parametervector[53]; /* 600 */
    aamp = paramdataPtr->parametervector[54]; /* 10 */
    bamp = paramdataPtr->parametervector[55]; /* 2.3667 */
    camp = paramdataPtr->parametervector[56]; /* 0.1012 */
    damp = paramdataPtr->parametervector[57]; /* 10.8688 */
    n1amp = paramdataPtr->parametervector[58]; /* 1.1456 */
    n2amp = paramdataPtr->parametervector[59]; /* 19.2166 */
    Kgseq = paramdataPtr->parametervector[60]; /* 460 */
    Kgsatp = paramdataPtr->parametervector[61]; /* 0.35 */
    Kgsglu = paramdataPtr->parametervector[62]; /* 4.1 */
    Kgsnh = paramdataPtr->parametervector[63]; /* 0.1 */
    Kgsadp = paramdataPtr->parametervector[64]; /* 0.0585 */
    Kgspi = paramdataPtr->parametervector[65]; /* 3.7 */
    Kgsgln = paramdataPtr->parametervector[66]; /* 5.65 */
    Keq = paramdataPtr->parametervector[67]; /* 460 */
    Vgludem = paramdataPtr->parametervector[68]; /* 120 */
    Kgludemglu = paramdataPtr->parametervector[69]; /* 8 */
    Kgludemeq = paramdataPtr->parametervector[70]; /* 1e+10 */
    Kgludemazglu = paramdataPtr->parametervector[71]; /* 0.5 */
    Vglndem = paramdataPtr->parametervector[72]; /* 70 */
    Kglndemgln = paramdataPtr->parametervector[73]; /* 2 */
    Kglndemeq = paramdataPtr->parametervector[74]; /* 1e+10 */
    Kglndemazgln = paramdataPtr->parametervector[75]; /* 0.25 */
    Vazglndem = paramdataPtr->parametervector[76]; /* 20 */
    Kazglndemazgln = paramdataPtr->parametervector[77]; /* 1 */
    Kazglndemeq = paramdataPtr->parametervector[78]; /* 1e+10 */
    Kazglndemazinter = paramdataPtr->parametervector[79]; /* 0.5 */
    Vazgludem = paramdataPtr->parametervector[80]; /* 30 */
    Kazgludemazglu = paramdataPtr->parametervector[81]; /* 0.3 */
    Kazgludemeq = paramdataPtr->parametervector[82]; /* 1e+10 */
    Kazgludemazinter = paramdataPtr->parametervector[83]; /* 0.5 */
    Vadp = paramdataPtr->parametervector[84]; /* 100 */
    Kadp = paramdataPtr->parametervector[85]; /* 0.5 */
    ATPtot = paramdataPtr->parametervector[86]; /* 5.37 */
    GStot = paramdataPtr->parametervector[87]; /* 0.014 */
    PIItot = paramdataPtr->parametervector[88]; /* 0.003 */
    Kd1 = paramdataPtr->parametervector[89]; /* 0.005 */
    Kd2 = paramdataPtr->parametervector[90]; /* 0.15 */
    Kd3 = paramdataPtr->parametervector[91]; /* 0.15 */
    Kd1piiump = paramdataPtr->parametervector[92]; /* 0.025 */
    Kd2piiump = paramdataPtr->parametervector[93]; /* 0.15 */
    Kd3piiump = paramdataPtr->parametervector[94]; /* 0.15 */
    compartment = paramdataPtr->parametervector[95]; /* 1 */
    vAPP_GS = aamp*camp/((1.0+pow(12.0,n1amp)*pow(AMP/(bamp*GStot),n1amp))*(1.0+pow(12.0,n2amp)*pow(AMP/(damp*GStot),n2amp)))*Vgs;
    nAMP = 12.0*AMP/GStot;
    PIIKG1 = 3.0*PII*KG/Kd1/(1.0+3.0*KG/Kd1+3.0*pow(KG,2.0)/(Kd1*Kd2)+pow(KG,3.0)/(Kd1*Kd2*Kd3));
    PIIUMP3KG3 = PIIUMP3*pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump)/(1.0+3.0*KG/Kd1piiump+3.0*pow(KG,2.0)/(Kd1piiump*Kd2piiump)+pow(KG,3.0)/(Kd1piiump*Kd2piiump*Kd3piiump));
    PII = 0.003;
    UTP = 0.5;
    PIIUMP = 0.0;
    PPi = 0.05;
    GLN = 1.0;
    PIIUMP2 = 0.0;
    PIIUMP3 = 0.0;
    UMP = 0.01;
    GS = 0.014;
    AMP = 0.0;
    NH4 = 0.05;
    KG = 0.2;
    NADPH = 0.15;
    GLU = 1.0;
    NADP = 0.05;
    AZGLU = 1.0;
    ATP = 2.685;
    ADP = 2.685;
    AZglu = 0.1;
    AZGLN = 1.0;
    AZgln = 0.1;
    icVector[0] = PII;
    icVector[1] = UTP;
    icVector[2] = PIIUMP;
    icVector[3] = PPi;
    icVector[4] = GLN;
    icVector[5] = PIIUMP2;
    icVector[6] = PIIUMP3;
    icVector[7] = UMP;
    icVector[8] = GS;
    icVector[9] = AMP;
    icVector[10] = NH4;
    icVector[11] = KG;
    icVector[12] = NADPH;
    icVector[13] = GLU;
    icVector[14] = NADP;
    icVector[15] = AZGLU;
    icVector[16] = ATP;
    icVector[17] = ADP;
    icVector[18] = AZglu;
    icVector[19] = AZGLN;
    icVector[20] = AZgln;
}

