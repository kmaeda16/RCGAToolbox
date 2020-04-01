#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "Maeda2019_RefinedActive_Kim_20190413_1_mod_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

static double Function_for_vamtb(double Kamtbnh,double NH4int,double NH4surf,double Vamtb_app,double default0,double phi)
{
    return Vamtb_app*(NH4surf-NH4int/phi)/(Kamtbnh+NH4surf)/default0;
}

static double Function_for_vdiff(double NH3int,double NH3surf,double default0,double kdiff)
{
    return kdiff*(NH3surf-NH3int)/default0;
}

static double Function_for_vdead(double GSAMP,double Kdeadgsamp,double Vdead_app,double default0)
{
    return Vdead_app*GSAMP/(Kdeadgsamp+GSAMP)/default0;
}

static double Function_for_vglndemf(double GLNdemf,double default0,double mu)
{
    return mu*GLNdemf/default0;
}

static double Function_for_vglndemn(double GLNdemn,double default0,double mu)
{
    return mu*GLNdemn/default0;
}

static double Function_for_vgludemf(double GLUdemf,double default0,double mu)
{
    return mu*GLUdemf/default0;
}

static double Function_for_vad(double GS,double Kadgs,double Vad_app,double default0)
{
    return Vad_app*GS/(Kadgs+GS)/default0;
}

static double Function_for_vgdh(double GLU,double Kgdheq,double Kgdhglu,double Kgdhnadp,double Kgdhnadph,double Kgdhnh,double Kgdhog,double NADP,double NADPH,double NH4int,double OG,double Vgdh,double default0)
{
    return Vgdh/(Kgdhog*Kgdhnh*Kgdhnadph)*(OG*NH4int*NADPH-GLU*NADP/Kgdheq)/((1.0+NH4int/Kgdhnh)*(1.0+OG/Kgdhog+GLU/Kgdhglu)*(1.0+NADPH/Kgdhnadph+NADP/Kgdhnadp))/default0;
}

static double Function_for_vgludemn(double GLUdemn,double default0,double mu)
{
    return mu*GLUdemn/default0;
}

static double Function_for_vutglnk2(double GLN,double GlnKUMP,double GlnKUMP2,double GlnKUMP3,double GlnK_AmtBfree,double Kutgln,double Kutglnk,double Kutiglnk,double Kutiglnkump,double Kutippi,double Kututp,double PPi,double UTP,double UTase,double default0,double kcatutglnk)
{
    return kcatutglnk*UTase*GlnKUMP*UTP/((1.0+GLN/Kutgln)*(Kutiglnk*Kututp+Kututp*(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)+Kutglnk*UTP+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP+Kutglnk*UTP*(GlnKUMP+GlnKUMP2+GlnKUMP3)/Kutiglnkump+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP*PPi/Kutippi))/default0;
}

static double Function_for_vutglnk3(double GLN,double GlnKUMP,double GlnKUMP2,double GlnKUMP3,double GlnK_AmtBfree,double Kutgln,double Kutglnk,double Kutiglnk,double Kutiglnkump,double Kutippi,double Kututp,double PPi,double UTP,double UTase,double default0,double kcatutglnk)
{
    return kcatutglnk*UTase*GlnKUMP2*UTP/((1.0+GLN/Kutgln)*(Kutiglnk*Kututp+Kututp*(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)+Kutglnk*UTP+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP+Kutglnk*UTP*(GlnKUMP+GlnKUMP2+GlnKUMP3)/Kutiglnkump+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP*PPi/Kutippi))/default0;
}

static double Function_for_vurglnb1(double GLN,double GlnBUMP,double GlnBUMP2,double GlnBUMP3,double Kurgln,double Kurglnbump,double Kurump,double UMP,double UTase,double default0,double kcaturglnb)
{
    return kcaturglnb*UTase*GlnBUMP/((1.0+Kurgln/GLN)*(Kurglnbump+(1.0+UMP/Kurump)*(GlnBUMP+GlnBUMP2+GlnBUMP3)))/default0;
}

static double Function_for_vurglnb2(double GLN,double GlnBUMP,double GlnBUMP2,double GlnBUMP3,double Kurgln,double Kurglnbump,double Kurump,double UMP,double UTase,double default0,double kcaturglnb)
{
    return kcaturglnb*UTase*GlnBUMP2/((1.0+Kurgln/GLN)*(Kurglnbump+(1.0+UMP/Kurump)*(GlnBUMP+GlnBUMP2+GlnBUMP3)))/default0;
}

static double Function_for_vutglnb3(double GLN,double GlnB,double GlnBUMP,double GlnBUMP2,double GlnBUMP3,double Kutgln,double Kutglnb,double Kutiglnb,double Kutiglnbump,double Kutippi,double Kututp,double PPi,double UTP,double UTase,double default0,double kcatutglnb)
{
    return kcatutglnb*UTase*GlnBUMP2*UTP/((1.0+GLN/Kutgln)*(Kutiglnb*Kututp+Kututp*(GlnB+GlnBUMP+GlnBUMP2)+Kutglnb*UTP+(GlnB+GlnBUMP+GlnBUMP2)*UTP+Kutglnb*UTP*(GlnBUMP+GlnBUMP2+GlnBUMP3)/Kutiglnbump+(GlnB+GlnBUMP+GlnBUMP2)*UTP*PPi/Kutippi))/default0;
}

static double Function_for_vurglnk3(double GLN,double GlnKUMP,double GlnKUMP2,double GlnKUMP3,double Kurgln,double Kurglnkump,double Kurump,double UMP,double UTase,double default0,double kcaturglnk)
{
    return kcaturglnk*UTase*GlnKUMP3/((1.0+Kurgln/GLN)*(Kurglnkump+(1.0+UMP/Kurump)*(GlnKUMP+GlnKUMP2+GlnKUMP3)))/default0;
}

static double Function_for_vutglnk1(double GLN,double GlnKUMP,double GlnKUMP2,double GlnKUMP3,double GlnK_AmtBfree,double Kutgln,double Kutglnk,double Kutiglnk,double Kutiglnkump,double Kutippi,double Kututp,double PPi,double UTP,double UTase,double default0,double kcatutglnk)
{
    return kcatutglnk*UTase*GlnK_AmtBfree*UTP/((1.0+GLN/Kutgln)*(Kutiglnk*Kututp+Kututp*(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)+Kutglnk*UTP+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP+Kutglnk*UTP*(GlnKUMP+GlnKUMP2+GlnKUMP3)/Kutiglnkump+(GlnK_AmtBfree+GlnKUMP+GlnKUMP2)*UTP*PPi/Kutippi))/default0;
}

static double Function_for_vutglnb1(double GLN,double GlnB,double GlnBUMP,double GlnBUMP2,double GlnBUMP3,double Kutgln,double Kutglnb,double Kutiglnb,double Kutiglnbump,double Kutippi,double Kututp,double PPi,double UTP,double UTase,double default0,double kcatutglnb)
{
    return kcatutglnb*UTase*GlnB*UTP/((1.0+GLN/Kutgln)*(Kutiglnb*Kututp+Kututp*(GlnB+GlnBUMP+GlnBUMP2)+Kutglnb*UTP+(GlnB+GlnBUMP+GlnBUMP2)*UTP+Kutglnb*UTP*(GlnBUMP+GlnBUMP2+GlnBUMP3)/Kutiglnbump+(GlnB+GlnBUMP+GlnBUMP2)*UTP*PPi/Kutippi))/default0;
}

static double Function_for_vutglnb2(double GLN,double GlnB,double GlnBUMP,double GlnBUMP2,double GlnBUMP3,double Kutgln,double Kutglnb,double Kutiglnb,double Kutiglnbump,double Kutippi,double Kututp,double PPi,double UTP,double UTase,double default0,double kcatutglnb)
{
    return kcatutglnb*UTase*GlnBUMP*UTP/((1.0+GLN/Kutgln)*(Kutiglnb*Kututp+Kututp*(GlnB+GlnBUMP+GlnBUMP2)+Kutglnb*UTP+(GlnB+GlnBUMP+GlnBUMP2)*UTP+Kutglnb*UTP*(GlnBUMP+GlnBUMP2+GlnBUMP3)/Kutiglnbump+(GlnB+GlnBUMP+GlnBUMP2)*UTP*PPi/Kutippi))/default0;
}

static double Function_for_vurglnb3(double GLN,double GlnBUMP,double GlnBUMP2,double GlnBUMP3,double Kurgln,double Kurglnbump,double Kurump,double UMP,double UTase,double default0,double kcaturglnb)
{
    return kcaturglnb*UTase*GlnBUMP3/((1.0+Kurgln/GLN)*(Kurglnbump+(1.0+UMP/Kurump)*(GlnBUMP+GlnBUMP2+GlnBUMP3)))/default0;
}

static double Function_for_vgog(double GLN,double GLU,double Kgoggln,double Kgogglu,double Kgognadp,double Kgognadph,double Kgogog,double NADP,double NADPH,double OG,double Vgog,double default0)
{
    return Vgog*GLN*OG*NADPH/(Kgoggln*Kgogog*Kgognadph)/((1.0+GLN/Kgoggln+GLU/Kgogglu)*(1.0+OG/Kgogog+GLU/Kgogglu)*(1.0+NADPH/Kgognadph+NADP/Kgognadp))/default0;
}

static double Function_for_vurglnk1(double GLN,double GlnKUMP,double GlnKUMP2,double GlnKUMP3,double Kurgln,double Kurglnkump,double Kurump,double UMP,double UTase,double default0,double kcaturglnk)
{
    return kcaturglnk*UTase*GlnKUMP/((1.0+Kurgln/GLN)*(Kurglnkump+(1.0+UMP/Kurump)*(GlnKUMP+GlnKUMP2+GlnKUMP3)))/default0;
}

static double Function_for_vgs(double ADP,double ATP,double GLN,double GLU,double Kgsadp,double Kgsatp,double Kgseq,double Kgsgln,double Kgsglu,double Kgsnh,double Kgspi,double NH4int,double Pi,double Vgs_app,double default0)
{
    return Vgs_app/(Kgsatp*Kgsnh*Kgsglu)*(ATP*NH4int*GLU-ADP*GLN*Pi/Kgseq)/((1.0+ATP/Kgsatp+ADP/Kgsadp+Pi/Kgspi+ADP*Pi/(Kgsadp*Kgspi))*(1.0+NH4int/Kgsnh+GLN/Kgsgln+GLU/Kgsglu+GLN*NH4int/(Kgsgln*Kgsnh)+GLU*NH4int/(Kgsglu*Kgsnh)))/default0;
}

static double Function_for_vurglnk2(double GLN,double GlnKUMP,double GlnKUMP2,double GlnKUMP3,double Kurgln,double Kurglnkump,double Kurump,double UMP,double UTase,double default0,double kcaturglnk)
{
    return kcaturglnk*UTase*GlnKUMP2/((1.0+Kurgln/GLN)*(Kurglnkump+(1.0+UMP/Kurump)*(GlnKUMP+GlnKUMP2+GlnKUMP3)))/default0;
}

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double ADP,ATP,GLN,GLU,GS,GSAMP,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,GlnK,GlnKUMP,GlnKUMP2,GlnKUMP3,NADP,NADPH,NHxint,PPi,Pi,UMP;
    double UTP;
    double Acell,AmtB,Dpsi,EXTERNAL,F,GLNdemf,GLNdemn,GLUdemf,GLUdemn,Kadgln,Kadglnbog,Kadgs,Kamtbnh,Kdeadgln,Kdeadglnbog,Kdeadglnbump,Kdeadgsamp,Kgdheq,Kgdhglu,Kgdhnadp;
    double Kgdhnadph,Kgdhnh,Kgdhog,Kglnbog1,Kglnbog2,Kglnbog3,Kglnbumpog1,Kglnbumpog2,Kglnbumpog3,Kglnkamtb,Kglnkog1,Kglnkog2,Kglnkog3,Kgoggln,Kgogglu,Kgognadp,Kgognadph,Kgogog,Kgrowthgln,Kgrowthglu;
    double Kgsadp,Kgsatp,Kgseq,Kgsgln,Kgsglu,Kgsnh,Kgspi,Kurgln,Kurglnbump,Kurglnkump,Kurump,Kutgln,Kutglnb,Kutglnk,Kutiglnb,Kutiglnbump,Kutiglnk,Kutiglnkump,Kutippi,Kututp;
    double NHxext,Nintstar,OGbasal,Pcm,R,T,UTase,Vad,Vcell,Vdead,Vgdh,Vgog,a1,aamp,b1,bamp,c1,camp,d1,damp;
    double e1,f1,g1,h1,i1,j1,k1,kappa,kcatamtb,kcatgs,kcaturglnb,kcaturglnk,kcatutglnb,kcatutglnk,kdb,l1,m1,n1,n1amp,n2amp;
    double o1,pHext,pHint,pKa,tau0,default0;
    double GStotal,Vgs,Ka,Hint,Hext,NH4int,OG,GlnBOG1,GlnB_OGfree,GlnKOG2,GlnBUMP3OG3,GlnK_OGfree,AmtB_GlnKfree,GlnKOG3,GlnKOG1,GlnKAmtB,GlnK_AmtBfree,NHxsurf,NH4surf,NH3int;
    double NH4ext,Vamtb_app,Vamtb,NH3ext,NH3surf,theta_ad,Vad_app,theta_dead,Vdead_app,nAMP,theta_gs,Vgs_app,phi,kdiff,tau,mu;
    double vad,vamtb,vdead,vdiff,vgdh,vglndemf,vglndemn,vgludemf,vgludemn,vgog,vgs,vurglnb1,vurglnb2,vurglnb3,vurglnk1,vurglnk2,vurglnk3,vutglnb1,vutglnb2,vutglnb3;
    double vutglnk1,vutglnk2,vutglnk3;

    time = time_local;

    ADP = stateVector[0];
    ATP = stateVector[1];
    GLN = stateVector[2];
    GLU = stateVector[3];
    GS = stateVector[4];
    GSAMP = stateVector[5];
    GlnB = stateVector[6];
    GlnBUMP = stateVector[7];
    GlnBUMP2 = stateVector[8];
    GlnBUMP3 = stateVector[9];
    GlnK = stateVector[10];
    GlnKUMP = stateVector[11];
    GlnKUMP2 = stateVector[12];
    GlnKUMP3 = stateVector[13];
    NADP = stateVector[14];
    NADPH = stateVector[15];
    NHxint = stateVector[16];
    PPi = stateVector[17];
    Pi = stateVector[18];
    UMP = stateVector[19];
    UTP = stateVector[20];
    Acell = paramdataPtr->parametervector[0]; /* 9.18e-12 */
    AmtB = paramdataPtr->parametervector[1]; /* 0.00168395 */
    Dpsi = paramdataPtr->parametervector[2]; /* -0.15 */
    EXTERNAL = paramdataPtr->parametervector[3]; /* 0 */
    F = paramdataPtr->parametervector[4]; /* 96485 */
    GLNdemf = paramdataPtr->parametervector[5]; /* 84.6574 */
    GLNdemn = paramdataPtr->parametervector[6]; /* 645.163 */
    GLUdemf = paramdataPtr->parametervector[7]; /* 316.076 */
    GLUdemn = paramdataPtr->parametervector[8]; /* 2116.81 */
    Kadgln = paramdataPtr->parametervector[9]; /* 0.952511 */
    Kadglnbog = paramdataPtr->parametervector[10]; /* 1.19929e-05 */
    Kadgs = paramdataPtr->parametervector[11]; /* 7.48531e-05 */
    Kamtbnh = paramdataPtr->parametervector[12]; /* 0.005 */
    Kdeadgln = paramdataPtr->parametervector[13]; /* 0.202364 */
    Kdeadglnbog = paramdataPtr->parametervector[14]; /* 2.58892e-06 */
    Kdeadglnbump = paramdataPtr->parametervector[15]; /* 7.07454e-06 */
    Kdeadgsamp = paramdataPtr->parametervector[16]; /* 0.00044786 */
    Kgdheq = paramdataPtr->parametervector[17]; /* 1290 */
    Kgdhglu = paramdataPtr->parametervector[18]; /* 6.27836 */
    Kgdhnadp = paramdataPtr->parametervector[19]; /* 0.0348261 */
    Kgdhnadph = paramdataPtr->parametervector[20]; /* 0.0485857 */
    Kgdhnh = paramdataPtr->parametervector[21]; /* 1.1 */
    Kgdhog = paramdataPtr->parametervector[22]; /* 0.518879 */
    Kglnbog1 = paramdataPtr->parametervector[23]; /* 0.0049519 */
    Kglnbog2 = paramdataPtr->parametervector[24]; /* 0.137803 */
    Kglnbog3 = paramdataPtr->parametervector[25]; /* 0.148798 */
    Kglnbumpog1 = paramdataPtr->parametervector[26]; /* 0.0241288 */
    Kglnbumpog2 = paramdataPtr->parametervector[27]; /* 0.151613 */
    Kglnbumpog3 = paramdataPtr->parametervector[28]; /* 0.132972 */
    Kglnkamtb = paramdataPtr->parametervector[29]; /* 6.26171e-08 */
    Kglnkog1 = paramdataPtr->parametervector[30]; /* 9.5247 */
    Kglnkog2 = paramdataPtr->parametervector[31]; /* 5.41511 */
    Kglnkog3 = paramdataPtr->parametervector[32]; /* 5.19229 */
    Kgoggln = paramdataPtr->parametervector[33]; /* 0.287442 */
    Kgogglu = paramdataPtr->parametervector[34]; /* 6.92079 */
    Kgognadp = paramdataPtr->parametervector[35]; /* 0.00339109 */
    Kgognadph = paramdataPtr->parametervector[36]; /* 0.00164105 */
    Kgogog = paramdataPtr->parametervector[37]; /* 0.006879 */
    Kgrowthgln = paramdataPtr->parametervector[38]; /* 1.50715 */
    Kgrowthglu = paramdataPtr->parametervector[39]; /* 31.1697 */
    Kgsadp = paramdataPtr->parametervector[40]; /* 0.0730688 */
    Kgsatp = paramdataPtr->parametervector[41]; /* 0.264057 */
    Kgseq = paramdataPtr->parametervector[42]; /* 460 */
    Kgsgln = paramdataPtr->parametervector[43]; /* 5.81152 */
    Kgsglu = paramdataPtr->parametervector[44]; /* 4.12683 */
    Kgsnh = paramdataPtr->parametervector[45]; /* 0.1 */
    Kgspi = paramdataPtr->parametervector[46]; /* 4.51385 */
    Kurgln = paramdataPtr->parametervector[47]; /* 0.063601 */
    Kurglnbump = paramdataPtr->parametervector[48]; /* 0.0046321 */
    Kurglnkump = paramdataPtr->parametervector[49]; /* 0.00197953 */
    Kurump = paramdataPtr->parametervector[50]; /* 8.33975 */
    Kutgln = paramdataPtr->parametervector[51]; /* 0.0640159 */
    Kutglnb = paramdataPtr->parametervector[52]; /* 0.00292567 */
    Kutglnk = paramdataPtr->parametervector[53]; /* 0.00456417 */
    Kutiglnb = paramdataPtr->parametervector[54]; /* 0.00185694 */
    Kutiglnbump = paramdataPtr->parametervector[55]; /* 0.00363233 */
    Kutiglnk = paramdataPtr->parametervector[56]; /* 0.00176846 */
    Kutiglnkump = paramdataPtr->parametervector[57]; /* 0.0024576 */
    Kutippi = paramdataPtr->parametervector[58]; /* 0.107936 */
    Kututp = paramdataPtr->parametervector[59]; /* 0.0417302 */
    NHxext = paramdataPtr->parametervector[60]; /* 0.00411274 */
    Nintstar = paramdataPtr->parametervector[61]; /* 0.033 */
    OGbasal = paramdataPtr->parametervector[62]; /* 0.550038 */
    Pcm = paramdataPtr->parametervector[63]; /* 0.0733905 */
    R = paramdataPtr->parametervector[64]; /* 8.314 */
    T = paramdataPtr->parametervector[65]; /* 310 */
    UTase = paramdataPtr->parametervector[66]; /* 0.0006 */
    Vad = paramdataPtr->parametervector[67]; /* 0.418331 */
    Vcell = paramdataPtr->parametervector[68]; /* 2.15e-18 */
    Vdead = paramdataPtr->parametervector[69]; /* 0.659368 */
    Vgdh = paramdataPtr->parametervector[70]; /* 352.64 */
    Vgog = paramdataPtr->parametervector[71]; /* 83.1933 */
    a1 = paramdataPtr->parametervector[72]; /* 1e-22 */
    aamp = paramdataPtr->parametervector[73]; /* 10 */
    b1 = paramdataPtr->parametervector[74]; /* 0.753263 */
    bamp = paramdataPtr->parametervector[75]; /* 2.3667 */
    c1 = paramdataPtr->parametervector[76]; /* 0.295923 */
    camp = paramdataPtr->parametervector[77]; /* 0.1012 */
    d1 = paramdataPtr->parametervector[78]; /* 0.0158453 */
    damp = paramdataPtr->parametervector[79]; /* 10.8688 */
    e1 = paramdataPtr->parametervector[80]; /* 1e-22 */
    f1 = paramdataPtr->parametervector[81]; /* 0.662946 */
    g1 = paramdataPtr->parametervector[82]; /* 14.4267 */
    h1 = paramdataPtr->parametervector[83]; /* 0.20749 */
    i1 = paramdataPtr->parametervector[84]; /* 1e-22 */
    j1 = paramdataPtr->parametervector[85]; /* 1e-22 */
    k1 = paramdataPtr->parametervector[86]; /* 1e-22 */
    kappa = paramdataPtr->parametervector[87]; /* 7.88981 */
    kcatamtb = paramdataPtr->parametervector[88]; /* 795355 */
    kcatgs = paramdataPtr->parametervector[89]; /* 53049.4 */
    kcaturglnb = paramdataPtr->parametervector[90]; /* 2.81579 */
    kcaturglnk = paramdataPtr->parametervector[91]; /* 17.3128 */
    kcatutglnb = paramdataPtr->parametervector[92]; /* 132.248 */
    kcatutglnk = paramdataPtr->parametervector[93]; /* 74.2685 */
    kdb = paramdataPtr->parametervector[94]; /* 13.9379 */
    l1 = paramdataPtr->parametervector[95]; /* 0.017405 */
    m1 = paramdataPtr->parametervector[96]; /* 0.87943 */
    n1 = paramdataPtr->parametervector[97]; /* 9.96306 */
    n1amp = paramdataPtr->parametervector[98]; /* 1.1456 */
    n2amp = paramdataPtr->parametervector[99]; /* 19.2166 */
    o1 = paramdataPtr->parametervector[100]; /* 1.29171 */
    pHext = paramdataPtr->parametervector[101]; /* 7.4 */
    pHint = paramdataPtr->parametervector[102]; /* 7.6 */
    pKa = paramdataPtr->parametervector[103]; /* 8.95 */
    tau0 = paramdataPtr->parametervector[104]; /* 45.8312 */
    default0 = paramdataPtr->parametervector[105]; /* 1 */
    GStotal = GS+GSAMP;
    Vgs = kcatgs*GStotal;
    Ka = pow(10.0,-pKa)*1000.0;
    Hint = pow(10.0,-pHint)*1000.0;
    Hext = pow(10.0,-pHext)*1000.0;
    NH4int = NHxint*Hint/(Ka+Hint);
    OG = kappa*(1.0-NH4int/Nintstar)+OGbasal;
    GlnBOG1 = 3.0*GlnB*OG/Kglnbog1/(1.0+3.0*OG/Kglnbog1+3.0*pow(OG,2.0)/(Kglnbog1*Kglnbog2)+pow(OG,3.0)/(Kglnbog1*Kglnbog2*Kglnbog3));
    GlnB_OGfree = GlnB/(1.0+3.0*OG/Kglnbog1+3.0*pow(OG,2.0)/(Kglnbog1*Kglnbog2)+pow(OG,3.0)/(Kglnbog1*Kglnbog2*Kglnbog3));
    GlnKOG2 = 3.0*GlnK*pow(OG,2.0)/(Kglnkog1*Kglnkog2)/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    GlnBUMP3OG3 = GlnBUMP3*pow(OG,3.0)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3)/(1.0+3.0*OG/Kglnbumpog1+3.0*pow(OG,2.0)/(Kglnbumpog1*Kglnbumpog2)+pow(OG,3.0)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3));
    GlnK_OGfree = GlnK/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    AmtB_GlnKfree = 0.5*(-GlnK_OGfree+AmtB-Kglnkamtb+pow((GlnK_OGfree-AmtB+Kglnkamtb)*(GlnK_OGfree-AmtB+Kglnkamtb)+4.0*Kglnkamtb*AmtB,0.5));
    GlnKOG3 = GlnK*pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3)/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    GlnKOG1 = 3.0*GlnK*OG/Kglnkog1/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    GlnKAmtB = AmtB-AmtB_GlnKfree;
    GlnK_AmtBfree = GlnK-GlnKAmtB;
    NHxsurf = NHxext;
    NH4surf = NHxsurf*Hext/(Ka+Hext);
    NH3int = NHxint*Ka/(Ka+Hint);
    NH4ext = NHxext*Hext/(Ka+Hext);
    Vamtb_app = kcatamtb*AmtB_GlnKfree;
    Vamtb = kcatamtb*AmtB;
    NH3ext = NHxext*Ka/(Ka+Hext);
    NH3surf = NHxsurf*Ka/(Ka+Hext);
    theta_ad = (a1*GlnBOG1/Kadglnbog+b1*GLN/Kadgln+c1*GlnBOG1*GLN/(Kadglnbog*Kadgln))/(1.0+GlnBOG1/Kadglnbog+GLN/Kadgln+GlnBOG1*GLN/(d1*Kadglnbog*Kadgln));
    Vad_app = theta_ad*Vad;
    theta_dead = (e1*GlnBOG1/Kdeadglnbog+f1*GLN/Kdeadgln+g1*GlnBUMP3OG3/Kdeadglnbump+h1*GlnBOG1*GLN/(Kdeadglnbog*Kdeadgln)+i1*GlnBOG1*GlnBUMP3OG3/(Kdeadglnbog*Kdeadglnbump)+j1*GLN*GlnBUMP3OG3/(Kdeadgln*Kdeadglnbump)+k1*GlnBOG1*GLN*GlnBUMP3OG3/(Kdeadglnbog*Kdeadgln*Kdeadglnbump))/(1.0+GlnBOG1/Kdeadglnbog+GLN/Kdeadgln+GlnBUMP3OG3/Kdeadglnbump+GlnBOG1*GLN/(l1*Kdeadglnbog*Kdeadgln)+GlnBOG1*GlnBUMP3OG3/(m1*Kdeadglnbog*Kdeadglnbump)+GLN*GlnBUMP3OG3/(n1*Kdeadgln*Kdeadglnbump)+GlnBOG1*GLN*GlnBUMP3OG3/(o1*Kdeadglnbog*Kdeadgln*Kdeadglnbump));
    Vdead_app = theta_dead*Vdead;
    nAMP = 12.0*GSAMP/GStotal;
    theta_gs = aamp/(1.0+pow(nAMP/bamp,n1amp))*camp/(1.0+pow(nAMP/damp,n2amp));
    Vgs_app = theta_gs*Vgs;
    phi = exp(-F*Dpsi/(R*T));
    kdiff = Pcm*Acell/Vcell;
    tau = tau0*(1.0+pow(Kgrowthglu/GLU,2.0)+pow(Kgrowthgln/GLN,2.0));
    mu = log(2.0)/tau;
    vad = default0*Function_for_vad(GS,Kadgs,Vad_app,default0);
    vamtb = default0*Function_for_vamtb(Kamtbnh,NH4int,NH4surf,Vamtb_app,default0,phi);
    vdead = default0*Function_for_vdead(GSAMP,Kdeadgsamp,Vdead_app,default0);
    vdiff = default0*Function_for_vdiff(NH3int,NH3surf,default0,kdiff);
    vgdh = default0*Function_for_vgdh(GLU,Kgdheq,Kgdhglu,Kgdhnadp,Kgdhnadph,Kgdhnh,Kgdhog,NADP,NADPH,NH4int,OG,Vgdh,default0);
    vglndemf = default0*Function_for_vglndemf(GLNdemf,default0,mu);
    vglndemn = default0*Function_for_vglndemn(GLNdemn,default0,mu);
    vgludemf = default0*Function_for_vgludemf(GLUdemf,default0,mu);
    vgludemn = default0*Function_for_vgludemn(GLUdemn,default0,mu);
    vgog = default0*Function_for_vgog(GLN,GLU,Kgoggln,Kgogglu,Kgognadp,Kgognadph,Kgogog,NADP,NADPH,OG,Vgog,default0);
    vgs = default0*Function_for_vgs(ADP,ATP,GLN,GLU,Kgsadp,Kgsatp,Kgseq,Kgsgln,Kgsglu,Kgsnh,Kgspi,NH4int,Pi,Vgs_app,default0);
    vurglnb1 = default0*Function_for_vurglnb1(GLN,GlnBUMP,GlnBUMP2,GlnBUMP3,Kurgln,Kurglnbump,Kurump,UMP,UTase,default0,kcaturglnb);
    vurglnb2 = default0*Function_for_vurglnb2(GLN,GlnBUMP,GlnBUMP2,GlnBUMP3,Kurgln,Kurglnbump,Kurump,UMP,UTase,default0,kcaturglnb);
    vurglnb3 = default0*Function_for_vurglnb3(GLN,GlnBUMP,GlnBUMP2,GlnBUMP3,Kurgln,Kurglnbump,Kurump,UMP,UTase,default0,kcaturglnb);
    vurglnk1 = default0*Function_for_vurglnk1(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,Kurgln,Kurglnkump,Kurump,UMP,UTase,default0,kcaturglnk);
    vurglnk2 = default0*Function_for_vurglnk2(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,Kurgln,Kurglnkump,Kurump,UMP,UTase,default0,kcaturglnk);
    vurglnk3 = default0*Function_for_vurglnk3(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,Kurgln,Kurglnkump,Kurump,UMP,UTase,default0,kcaturglnk);
    vutglnb1 = default0*Function_for_vutglnb1(GLN,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,Kutgln,Kutglnb,Kutiglnb,Kutiglnbump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnb);
    vutglnb2 = default0*Function_for_vutglnb2(GLN,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,Kutgln,Kutglnb,Kutiglnb,Kutiglnbump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnb);
    vutglnb3 = default0*Function_for_vutglnb3(GLN,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,Kutgln,Kutglnb,Kutiglnb,Kutiglnbump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnb);
    vutglnk1 = default0*Function_for_vutglnk1(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,GlnK_AmtBfree,Kutgln,Kutglnk,Kutiglnk,Kutiglnkump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnk);
    vutglnk2 = default0*Function_for_vutglnk2(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,GlnK_AmtBfree,Kutgln,Kutglnk,Kutiglnk,Kutiglnkump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnk);
    vutglnk3 = default0*Function_for_vutglnk3(GLN,GlnKUMP,GlnKUMP2,GlnKUMP3,GlnK_AmtBfree,Kutgln,Kutglnk,Kutiglnk,Kutiglnkump,Kutippi,Kututp,PPi,UTP,UTase,default0,kcatutglnk);
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = 0.0;
    	DDTvector[1] = 0.0;
    	DDTvector[2] = (-vglndemf-vglndemn-vgog+vgs)/default0;
    	DDTvector[3] = (+vgdh+vglndemn-vgludemf-vgludemn+2.0*vgog-vgs)/default0;
    	DDTvector[4] = (-vad+vdead)/default0;
    	DDTvector[5] = (+vad-vdead)/default0;
    	DDTvector[6] = (+vurglnb1-vutglnb1)/default0;
    	DDTvector[7] = (-vurglnb1+vurglnb2+vutglnb1-vutglnb2)/default0;
    	DDTvector[8] = (-vurglnb2+vurglnb3+vutglnb2-vutglnb3)/default0;
    	DDTvector[9] = (-vurglnb3+vutglnb3)/default0;
    	DDTvector[10] = (+vurglnk1-vutglnk1)/default0;
    	DDTvector[11] = (-vurglnk1+vurglnk2+vutglnk1-vutglnk2)/default0;
    	DDTvector[12] = (-vurglnk2+vurglnk3+vutglnk2-vutglnk3)/default0;
    	DDTvector[13] = (-vurglnk3+vutglnk3)/default0;
    	DDTvector[14] = 0.0;
    	DDTvector[15] = 0.0;
    	DDTvector[16] = (+vamtb+vdiff-vgdh-vgs)/default0;
    	DDTvector[17] = 0.0;
    	DDTvector[18] = 0.0;
    	DDTvector[19] = 0.0;
    	DDTvector[20] = 0.0;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = GStotal;
        variableVector[1] = Vgs;
        variableVector[2] = Ka;
        variableVector[3] = Hint;
        variableVector[4] = Hext;
        variableVector[5] = NH4int;
        variableVector[6] = OG;
        variableVector[7] = GlnBOG1;
        variableVector[8] = GlnB_OGfree;
        variableVector[9] = GlnKOG2;
        variableVector[10] = GlnBUMP3OG3;
        variableVector[11] = GlnK_OGfree;
        variableVector[12] = AmtB_GlnKfree;
        variableVector[13] = GlnKOG3;
        variableVector[14] = GlnKOG1;
        variableVector[15] = GlnKAmtB;
        variableVector[16] = GlnK_AmtBfree;
        variableVector[17] = NHxsurf;
        variableVector[18] = NH4surf;
        variableVector[19] = NH3int;
        variableVector[20] = NH4ext;
        variableVector[21] = Vamtb_app;
        variableVector[22] = Vamtb;
        variableVector[23] = NH3ext;
        variableVector[24] = NH3surf;
        variableVector[25] = theta_ad;
        variableVector[26] = Vad_app;
        variableVector[27] = theta_dead;
        variableVector[28] = Vdead_app;
        variableVector[29] = nAMP;
        variableVector[30] = theta_gs;
        variableVector[31] = Vgs_app;
        variableVector[32] = phi;
        variableVector[33] = kdiff;
        variableVector[34] = tau;
        variableVector[35] = mu;
        reactionVector[0] = vad;
        reactionVector[1] = vamtb;
        reactionVector[2] = vdead;
        reactionVector[3] = vdiff;
        reactionVector[4] = vgdh;
        reactionVector[5] = vglndemf;
        reactionVector[6] = vglndemn;
        reactionVector[7] = vgludemf;
        reactionVector[8] = vgludemn;
        reactionVector[9] = vgog;
        reactionVector[10] = vgs;
        reactionVector[11] = vurglnb1;
        reactionVector[12] = vurglnb2;
        reactionVector[13] = vurglnb3;
        reactionVector[14] = vurglnk1;
        reactionVector[15] = vurglnk2;
        reactionVector[16] = vurglnk3;
        reactionVector[17] = vutglnb1;
        reactionVector[18] = vutglnb2;
        reactionVector[19] = vutglnb3;
        reactionVector[20] = vutglnk1;
        reactionVector[21] = vutglnk2;
        reactionVector[22] = vutglnk3;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double ADP,ATP,GLN,GLU,GS,GSAMP,GlnB,GlnBUMP,GlnBUMP2,GlnBUMP3,GlnK,GlnKUMP,GlnKUMP2,GlnKUMP3,NADP,NADPH,NHxint,PPi,Pi,UMP;
    double UTP;
    double Acell,AmtB,Dpsi,EXTERNAL,F,GLNdemf,GLNdemn,GLUdemf,GLUdemn,Kadgln,Kadglnbog,Kadgs,Kamtbnh,Kdeadgln,Kdeadglnbog,Kdeadglnbump,Kdeadgsamp,Kgdheq,Kgdhglu,Kgdhnadp;
    double Kgdhnadph,Kgdhnh,Kgdhog,Kglnbog1,Kglnbog2,Kglnbog3,Kglnbumpog1,Kglnbumpog2,Kglnbumpog3,Kglnkamtb,Kglnkog1,Kglnkog2,Kglnkog3,Kgoggln,Kgogglu,Kgognadp,Kgognadph,Kgogog,Kgrowthgln,Kgrowthglu;
    double Kgsadp,Kgsatp,Kgseq,Kgsgln,Kgsglu,Kgsnh,Kgspi,Kurgln,Kurglnbump,Kurglnkump,Kurump,Kutgln,Kutglnb,Kutglnk,Kutiglnb,Kutiglnbump,Kutiglnk,Kutiglnkump,Kutippi,Kututp;
    double NHxext,Nintstar,OGbasal,Pcm,R,T,UTase,Vad,Vcell,Vdead,Vgdh,Vgog,a1,aamp,b1,bamp,c1,camp,d1,damp;
    double e1,f1,g1,h1,i1,j1,k1,kappa,kcatamtb,kcatgs,kcaturglnb,kcaturglnk,kcatutglnb,kcatutglnk,kdb,l1,m1,n1,n1amp,n2amp;
    double o1,pHext,pHint,pKa,tau0,default0;
    double GStotal,Vgs,Ka,Hint,Hext,NH4int,OG,GlnBOG1,GlnB_OGfree,GlnKOG2,GlnBUMP3OG3,GlnK_OGfree,AmtB_GlnKfree,GlnKOG3,GlnKOG1,GlnKAmtB,GlnK_AmtBfree,NHxsurf,NH4surf,NH3int;
    double NH4ext,Vamtb_app,Vamtb,NH3ext,NH3surf,theta_ad,Vad_app,theta_dead,Vdead_app,nAMP,theta_gs,Vgs_app,phi,kdiff,tau,mu;
    Acell = paramdataPtr->parametervector[0]; /* 9.18e-12 */
    AmtB = paramdataPtr->parametervector[1]; /* 0.00168395 */
    Dpsi = paramdataPtr->parametervector[2]; /* -0.15 */
    EXTERNAL = paramdataPtr->parametervector[3]; /* 0 */
    F = paramdataPtr->parametervector[4]; /* 96485 */
    GLNdemf = paramdataPtr->parametervector[5]; /* 84.6574 */
    GLNdemn = paramdataPtr->parametervector[6]; /* 645.163 */
    GLUdemf = paramdataPtr->parametervector[7]; /* 316.076 */
    GLUdemn = paramdataPtr->parametervector[8]; /* 2116.81 */
    Kadgln = paramdataPtr->parametervector[9]; /* 0.952511 */
    Kadglnbog = paramdataPtr->parametervector[10]; /* 1.19929e-05 */
    Kadgs = paramdataPtr->parametervector[11]; /* 7.48531e-05 */
    Kamtbnh = paramdataPtr->parametervector[12]; /* 0.005 */
    Kdeadgln = paramdataPtr->parametervector[13]; /* 0.202364 */
    Kdeadglnbog = paramdataPtr->parametervector[14]; /* 2.58892e-06 */
    Kdeadglnbump = paramdataPtr->parametervector[15]; /* 7.07454e-06 */
    Kdeadgsamp = paramdataPtr->parametervector[16]; /* 0.00044786 */
    Kgdheq = paramdataPtr->parametervector[17]; /* 1290 */
    Kgdhglu = paramdataPtr->parametervector[18]; /* 6.27836 */
    Kgdhnadp = paramdataPtr->parametervector[19]; /* 0.0348261 */
    Kgdhnadph = paramdataPtr->parametervector[20]; /* 0.0485857 */
    Kgdhnh = paramdataPtr->parametervector[21]; /* 1.1 */
    Kgdhog = paramdataPtr->parametervector[22]; /* 0.518879 */
    Kglnbog1 = paramdataPtr->parametervector[23]; /* 0.0049519 */
    Kglnbog2 = paramdataPtr->parametervector[24]; /* 0.137803 */
    Kglnbog3 = paramdataPtr->parametervector[25]; /* 0.148798 */
    Kglnbumpog1 = paramdataPtr->parametervector[26]; /* 0.0241288 */
    Kglnbumpog2 = paramdataPtr->parametervector[27]; /* 0.151613 */
    Kglnbumpog3 = paramdataPtr->parametervector[28]; /* 0.132972 */
    Kglnkamtb = paramdataPtr->parametervector[29]; /* 6.26171e-08 */
    Kglnkog1 = paramdataPtr->parametervector[30]; /* 9.5247 */
    Kglnkog2 = paramdataPtr->parametervector[31]; /* 5.41511 */
    Kglnkog3 = paramdataPtr->parametervector[32]; /* 5.19229 */
    Kgoggln = paramdataPtr->parametervector[33]; /* 0.287442 */
    Kgogglu = paramdataPtr->parametervector[34]; /* 6.92079 */
    Kgognadp = paramdataPtr->parametervector[35]; /* 0.00339109 */
    Kgognadph = paramdataPtr->parametervector[36]; /* 0.00164105 */
    Kgogog = paramdataPtr->parametervector[37]; /* 0.006879 */
    Kgrowthgln = paramdataPtr->parametervector[38]; /* 1.50715 */
    Kgrowthglu = paramdataPtr->parametervector[39]; /* 31.1697 */
    Kgsadp = paramdataPtr->parametervector[40]; /* 0.0730688 */
    Kgsatp = paramdataPtr->parametervector[41]; /* 0.264057 */
    Kgseq = paramdataPtr->parametervector[42]; /* 460 */
    Kgsgln = paramdataPtr->parametervector[43]; /* 5.81152 */
    Kgsglu = paramdataPtr->parametervector[44]; /* 4.12683 */
    Kgsnh = paramdataPtr->parametervector[45]; /* 0.1 */
    Kgspi = paramdataPtr->parametervector[46]; /* 4.51385 */
    Kurgln = paramdataPtr->parametervector[47]; /* 0.063601 */
    Kurglnbump = paramdataPtr->parametervector[48]; /* 0.0046321 */
    Kurglnkump = paramdataPtr->parametervector[49]; /* 0.00197953 */
    Kurump = paramdataPtr->parametervector[50]; /* 8.33975 */
    Kutgln = paramdataPtr->parametervector[51]; /* 0.0640159 */
    Kutglnb = paramdataPtr->parametervector[52]; /* 0.00292567 */
    Kutglnk = paramdataPtr->parametervector[53]; /* 0.00456417 */
    Kutiglnb = paramdataPtr->parametervector[54]; /* 0.00185694 */
    Kutiglnbump = paramdataPtr->parametervector[55]; /* 0.00363233 */
    Kutiglnk = paramdataPtr->parametervector[56]; /* 0.00176846 */
    Kutiglnkump = paramdataPtr->parametervector[57]; /* 0.0024576 */
    Kutippi = paramdataPtr->parametervector[58]; /* 0.107936 */
    Kututp = paramdataPtr->parametervector[59]; /* 0.0417302 */
    NHxext = paramdataPtr->parametervector[60]; /* 0.00411274 */
    Nintstar = paramdataPtr->parametervector[61]; /* 0.033 */
    OGbasal = paramdataPtr->parametervector[62]; /* 0.550038 */
    Pcm = paramdataPtr->parametervector[63]; /* 0.0733905 */
    R = paramdataPtr->parametervector[64]; /* 8.314 */
    T = paramdataPtr->parametervector[65]; /* 310 */
    UTase = paramdataPtr->parametervector[66]; /* 0.0006 */
    Vad = paramdataPtr->parametervector[67]; /* 0.418331 */
    Vcell = paramdataPtr->parametervector[68]; /* 2.15e-18 */
    Vdead = paramdataPtr->parametervector[69]; /* 0.659368 */
    Vgdh = paramdataPtr->parametervector[70]; /* 352.64 */
    Vgog = paramdataPtr->parametervector[71]; /* 83.1933 */
    a1 = paramdataPtr->parametervector[72]; /* 1e-22 */
    aamp = paramdataPtr->parametervector[73]; /* 10 */
    b1 = paramdataPtr->parametervector[74]; /* 0.753263 */
    bamp = paramdataPtr->parametervector[75]; /* 2.3667 */
    c1 = paramdataPtr->parametervector[76]; /* 0.295923 */
    camp = paramdataPtr->parametervector[77]; /* 0.1012 */
    d1 = paramdataPtr->parametervector[78]; /* 0.0158453 */
    damp = paramdataPtr->parametervector[79]; /* 10.8688 */
    e1 = paramdataPtr->parametervector[80]; /* 1e-22 */
    f1 = paramdataPtr->parametervector[81]; /* 0.662946 */
    g1 = paramdataPtr->parametervector[82]; /* 14.4267 */
    h1 = paramdataPtr->parametervector[83]; /* 0.20749 */
    i1 = paramdataPtr->parametervector[84]; /* 1e-22 */
    j1 = paramdataPtr->parametervector[85]; /* 1e-22 */
    k1 = paramdataPtr->parametervector[86]; /* 1e-22 */
    kappa = paramdataPtr->parametervector[87]; /* 7.88981 */
    kcatamtb = paramdataPtr->parametervector[88]; /* 795355 */
    kcatgs = paramdataPtr->parametervector[89]; /* 53049.4 */
    kcaturglnb = paramdataPtr->parametervector[90]; /* 2.81579 */
    kcaturglnk = paramdataPtr->parametervector[91]; /* 17.3128 */
    kcatutglnb = paramdataPtr->parametervector[92]; /* 132.248 */
    kcatutglnk = paramdataPtr->parametervector[93]; /* 74.2685 */
    kdb = paramdataPtr->parametervector[94]; /* 13.9379 */
    l1 = paramdataPtr->parametervector[95]; /* 0.017405 */
    m1 = paramdataPtr->parametervector[96]; /* 0.87943 */
    n1 = paramdataPtr->parametervector[97]; /* 9.96306 */
    n1amp = paramdataPtr->parametervector[98]; /* 1.1456 */
    n2amp = paramdataPtr->parametervector[99]; /* 19.2166 */
    o1 = paramdataPtr->parametervector[100]; /* 1.29171 */
    pHext = paramdataPtr->parametervector[101]; /* 7.4 */
    pHint = paramdataPtr->parametervector[102]; /* 7.6 */
    pKa = paramdataPtr->parametervector[103]; /* 8.95 */
    tau0 = paramdataPtr->parametervector[104]; /* 45.8312 */
    default0 = paramdataPtr->parametervector[105]; /* 1 */
    GStotal = GS+GSAMP;
    Vgs = kcatgs*GStotal;
    Ka = pow(10.0,-pKa)*1000.0;
    Hint = pow(10.0,-pHint)*1000.0;
    Hext = pow(10.0,-pHext)*1000.0;
    NH4int = NHxint*Hint/(Ka+Hint);
    OG = kappa*(1.0-NH4int/Nintstar)+OGbasal;
    GlnBOG1 = 3.0*GlnB*OG/Kglnbog1/(1.0+3.0*OG/Kglnbog1+3.0*pow(OG,2.0)/(Kglnbog1*Kglnbog2)+pow(OG,3.0)/(Kglnbog1*Kglnbog2*Kglnbog3));
    GlnB_OGfree = GlnB/(1.0+3.0*OG/Kglnbog1+3.0*pow(OG,2.0)/(Kglnbog1*Kglnbog2)+pow(OG,3.0)/(Kglnbog1*Kglnbog2*Kglnbog3));
    GlnKOG2 = 3.0*GlnK*pow(OG,2.0)/(Kglnkog1*Kglnkog2)/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    GlnBUMP3OG3 = GlnBUMP3*pow(OG,3.0)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3)/(1.0+3.0*OG/Kglnbumpog1+3.0*pow(OG,2.0)/(Kglnbumpog1*Kglnbumpog2)+pow(OG,3.0)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3));
    GlnK_OGfree = GlnK/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    AmtB_GlnKfree = 0.5*(-GlnK_OGfree+AmtB-Kglnkamtb+pow((GlnK_OGfree-AmtB+Kglnkamtb)*(GlnK_OGfree-AmtB+Kglnkamtb)+4.0*Kglnkamtb*AmtB,0.5));
    GlnKOG3 = GlnK*pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3)/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    GlnKOG1 = 3.0*GlnK*OG/Kglnkog1/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3));
    GlnKAmtB = AmtB-AmtB_GlnKfree;
    GlnK_AmtBfree = GlnK-GlnKAmtB;
    NHxsurf = NHxext;
    NH4surf = NHxsurf*Hext/(Ka+Hext);
    NH3int = NHxint*Ka/(Ka+Hint);
    NH4ext = NHxext*Hext/(Ka+Hext);
    Vamtb_app = kcatamtb*AmtB_GlnKfree;
    Vamtb = kcatamtb*AmtB;
    NH3ext = NHxext*Ka/(Ka+Hext);
    NH3surf = NHxsurf*Ka/(Ka+Hext);
    theta_ad = (a1*GlnBOG1/Kadglnbog+b1*GLN/Kadgln+c1*GlnBOG1*GLN/(Kadglnbog*Kadgln))/(1.0+GlnBOG1/Kadglnbog+GLN/Kadgln+GlnBOG1*GLN/(d1*Kadglnbog*Kadgln));
    Vad_app = theta_ad*Vad;
    theta_dead = (e1*GlnBOG1/Kdeadglnbog+f1*GLN/Kdeadgln+g1*GlnBUMP3OG3/Kdeadglnbump+h1*GlnBOG1*GLN/(Kdeadglnbog*Kdeadgln)+i1*GlnBOG1*GlnBUMP3OG3/(Kdeadglnbog*Kdeadglnbump)+j1*GLN*GlnBUMP3OG3/(Kdeadgln*Kdeadglnbump)+k1*GlnBOG1*GLN*GlnBUMP3OG3/(Kdeadglnbog*Kdeadgln*Kdeadglnbump))/(1.0+GlnBOG1/Kdeadglnbog+GLN/Kdeadgln+GlnBUMP3OG3/Kdeadglnbump+GlnBOG1*GLN/(l1*Kdeadglnbog*Kdeadgln)+GlnBOG1*GlnBUMP3OG3/(m1*Kdeadglnbog*Kdeadglnbump)+GLN*GlnBUMP3OG3/(n1*Kdeadgln*Kdeadglnbump)+GlnBOG1*GLN*GlnBUMP3OG3/(o1*Kdeadglnbog*Kdeadgln*Kdeadglnbump));
    Vdead_app = theta_dead*Vdead;
    nAMP = 12.0*GSAMP/GStotal;
    theta_gs = aamp/(1.0+pow(nAMP/bamp,n1amp))*camp/(1.0+pow(nAMP/damp,n2amp));
    Vgs_app = theta_gs*Vgs;
    phi = exp(-F*Dpsi/(R*T));
    kdiff = Pcm*Acell/Vcell;
    tau = tau0*(1.0+pow(Kgrowthglu/GLU,2.0)+pow(Kgrowthgln/GLN,2.0));
    mu = log(2.0)/tau;
    ADP = 0.56;
    ATP = 9.6;
    GLN = 2.0;
    GLU = 80.0;
    GS = 0.0108376;
    GSAMP = 0.0;
    GlnB = 0.00065;
    GlnBUMP = 0.0;
    GlnBUMP2 = 0.0;
    GlnBUMP3 = 0.0;
    GlnK = 0.00203338;
    GlnKUMP = 0.0;
    GlnKUMP2 = 0.0;
    GlnKUMP3 = 0.0;
    NADP = 0.076;
    NADPH = 0.12;
    NHxint = 0.02;
    PPi = 0.05;
    Pi = 10.0;
    UMP = 0.01;
    UTP = 8.3;
    icVector[0] = ADP;
    icVector[1] = ATP;
    icVector[2] = GLN;
    icVector[3] = GLU;
    icVector[4] = GS;
    icVector[5] = GSAMP;
    icVector[6] = GlnB;
    icVector[7] = GlnBUMP;
    icVector[8] = GlnBUMP2;
    icVector[9] = GlnBUMP3;
    icVector[10] = GlnK;
    icVector[11] = GlnKUMP;
    icVector[12] = GlnKUMP2;
    icVector[13] = GlnKUMP3;
    icVector[14] = NADP;
    icVector[15] = NADPH;
    icVector[16] = NHxint;
    icVector[17] = PPi;
    icVector[18] = Pi;
    icVector[19] = UMP;
    icVector[20] = UTP;
}

