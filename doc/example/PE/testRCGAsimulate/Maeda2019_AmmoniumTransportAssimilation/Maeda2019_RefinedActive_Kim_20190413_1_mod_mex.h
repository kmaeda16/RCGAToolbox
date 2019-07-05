#include "mex.h"

const int NRSTATES = 21;
const int NRPARAMETERS = 106;
const int NRVARIABLES = 36;
const int NRREACTIONS = 23;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[21] = {
	0.56,9.6,2,80,0.0108376,0,0.00065,0,0,0,0.00203338,0,0,0,0.076,0.12,0.02,0.05,10,0.01,
	8.3};
char *defaultICs_nonnum[1];

double defaultParam[106] = {
	9.18e-12,0.00168395,-0.15,0,96485,84.6574,645.163,316.076,2116.81,0.952511,1.19929e-05,7.48531e-05,0.005,0.202364,2.58892e-06,7.07454e-06,0.00044786,1290,6.27836,0.0348261,
	0.0485857,1.1,0.518879,0.0049519,0.137803,0.148798,0.0241288,0.151613,0.132972,6.26171e-08,9.5247,5.41511,5.19229,0.287442,6.92079,0.00339109,0.00164105,0.006879,1.50715,31.1697,
	0.0730688,0.264057,460,5.81152,4.12683,0.1,4.51385,0.063601,0.0046321,0.00197953,8.33975,0.0640159,0.00292567,0.00456417,0.00185694,0.00363233,0.00176846,0.0024576,0.107936,0.0417302,
	0.00411274,0.033,0.550038,0.0733905,8.314,310,0.0006,0.418331,2.15e-18,0.659368,352.64,83.1933,1e-22,10,0.753263,2.3667,0.295923,0.1012,0.0158453,10.8688,
	1e-22,0.662946,14.4267,0.20749,1e-22,1e-22,1e-22,7.88981,795355,53049.4,2.81579,17.3128,132.248,74.2685,13.9379,0.017405,0.87943,9.96306,1.1456,19.2166,
	1.29171,7.4,7.6,8.95,45.8312,1};
char *stateNames[21] = {
	"ADP","ATP","GLN","GLU","GS","GSAMP","GlnB","GlnBUMP","GlnBUMP2","GlnBUMP3","GlnK","GlnKUMP","GlnKUMP2","GlnKUMP3","NADP","NADPH","NHxint","PPi","Pi","UMP",
	"UTP"};
char *parameterNames[106] = {
	"Acell","AmtB","Dpsi","EXTERNAL","F","GLNdemf","GLNdemn","GLUdemf","GLUdemn","Kadgln","Kadglnbog","Kadgs","Kamtbnh","Kdeadgln","Kdeadglnbog","Kdeadglnbump","Kdeadgsamp","Kgdheq","Kgdhglu","Kgdhnadp",
	"Kgdhnadph","Kgdhnh","Kgdhog","Kglnbog1","Kglnbog2","Kglnbog3","Kglnbumpog1","Kglnbumpog2","Kglnbumpog3","Kglnkamtb","Kglnkog1","Kglnkog2","Kglnkog3","Kgoggln","Kgogglu","Kgognadp","Kgognadph","Kgogog","Kgrowthgln","Kgrowthglu",
	"Kgsadp","Kgsatp","Kgseq","Kgsgln","Kgsglu","Kgsnh","Kgspi","Kurgln","Kurglnbump","Kurglnkump","Kurump","Kutgln","Kutglnb","Kutglnk","Kutiglnb","Kutiglnbump","Kutiglnk","Kutiglnkump","Kutippi","Kututp",
	"NHxext","Nintstar","OGbasal","Pcm","R","T","UTase","Vad","Vcell","Vdead","Vgdh","Vgog","a1","aamp","b1","bamp","c1","camp","d1","damp",
	"e1","f1","g1","h1","i1","j1","k1","kappa","kcatamtb","kcatgs","kcaturglnb","kcaturglnk","kcatutglnb","kcatutglnk","kdb","l1","m1","n1","n1amp","n2amp",
	"o1","pHext","pHint","pKa","tau0","default0"};
char *variableNames[36] = {
	"GStotal","Vgs","Ka","Hint","Hext","NH4int","OG","GlnBOG1","GlnB_OGfree","GlnKOG2","GlnBUMP3OG3","GlnK_OGfree","AmtB_GlnKfree","GlnKOG3","GlnKOG1","GlnKAmtB","GlnK_AmtBfree","NHxsurf","NH4surf","NH3int",
	"NH4ext","Vamtb_app","Vamtb","NH3ext","NH3surf","theta_ad","Vad_app","theta_dead","Vdead_app","nAMP","theta_gs","Vgs_app","phi","kdiff","tau","mu"};
char *variableFormulas[36] = {
	"GS+GSAMP","kcatgs*GStotal","pow(10.0,-pKa)*1000.0","pow(10.0,-pHint)*1000.0","pow(10.0,-pHext)*1000.0","NHxint*Hint/(Ka+Hint)","kappa*(1.0-NH4int/Nintstar)+OGbasal","3.0*GlnB*OG/Kglnbog1/(1.0+3.0*OG/Kglnbog1+3.0*pow(OG,2.0)/(Kglnbog1*Kglnbog2)+pow(OG,3.0)/(Kglnbog1*Kglnbog2*Kglnbog3))","GlnB/(1.0+3.0*OG/Kglnbog1+3.0*pow(OG,2.0)/(Kglnbog1*Kglnbog2)+pow(OG,3.0)/(Kglnbog1*Kglnbog2*Kglnbog3))","3.0*GlnK*pow(OG,2.0)/(Kglnkog1*Kglnkog2)/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3))","GlnBUMP3*pow(OG,3.0)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3)/(1.0+3.0*OG/Kglnbumpog1+3.0*pow(OG,2.0)/(Kglnbumpog1*Kglnbumpog2)+pow(OG,3.0)/(Kglnbumpog1*Kglnbumpog2*Kglnbumpog3))","GlnK/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3))","0.5*(-GlnK_OGfree+AmtB-Kglnkamtb+pow((GlnK_OGfree-AmtB+Kglnkamtb)*(GlnK_OGfree-AmtB+Kglnkamtb)+4.0*Kglnkamtb*AmtB,0.5))","GlnK*pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3)/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3))","3.0*GlnK*OG/Kglnkog1/(1.0+3.0*OG/Kglnkog1+3.0*pow(OG,2.0)/(Kglnkog1*Kglnkog2)+pow(OG,3.0)/(Kglnkog1*Kglnkog2*Kglnkog3))","AmtB-AmtB_GlnKfree","GlnK-GlnKAmtB","NHxext","NHxsurf*Hext/(Ka+Hext)","NHxint*Ka/(Ka+Hint)",
	"NHxext*Hext/(Ka+Hext)","kcatamtb*AmtB_GlnKfree","kcatamtb*AmtB","NHxext*Ka/(Ka+Hext)","NHxsurf*Ka/(Ka+Hext)","(a1*GlnBOG1/Kadglnbog+b1*GLN/Kadgln+c1*GlnBOG1*GLN/(Kadglnbog*Kadgln))/(1.0+GlnBOG1/Kadglnbog+GLN/Kadgln+GlnBOG1*GLN/(d1*Kadglnbog*Kadgln))","theta_ad*Vad","(e1*GlnBOG1/Kdeadglnbog+f1*GLN/Kdeadgln+g1*GlnBUMP3OG3/Kdeadglnbump+h1*GlnBOG1*GLN/(Kdeadglnbog*Kdeadgln)+i1*GlnBOG1*GlnBUMP3OG3/(Kdeadglnbog*Kdeadglnbump)+j1*GLN*GlnBUMP3OG3/(Kdeadgln*Kdeadglnbump)+k1*GlnBOG1*GLN*GlnBUMP3OG3/(Kdeadglnbog*Kdeadgln*Kdeadglnbump))/(1.0+GlnBOG1/Kdeadglnbog+GLN/Kdeadgln+GlnBUMP3OG3/Kdeadglnbump+GlnBOG1*GLN/(l1*Kdeadglnbog*Kdeadgln)+GlnBOG1*GlnBUMP3OG3/(m1*Kdeadglnbog*Kdeadglnbump)+GLN*GlnBUMP3OG3/(n1*Kdeadgln*Kdeadglnbump)+GlnBOG1*GLN*GlnBUMP3OG3/(o1*Kdeadglnbog*Kdeadgln*Kdeadglnbump))","theta_dead*Vdead","12.0*GSAMP/GStotal","aamp/(1.0+pow(nAMP/bamp,n1amp))*camp/(1.0+pow(nAMP/damp,n2amp))","theta_gs*Vgs","exp(-F*Dpsi/(R*T))","Pcm*Acell/Vcell","tau0*(1.0+pow(Kgrowthglu/GLU,2.0)+pow(Kgrowthgln/GLN,2.0))","log(2.0)/tau"};
char *reactionNames[23] = {
	"vad","vamtb","vdead","vdiff","vgdh","vglndemf","vglndemn","vgludemf","vgludemn","vgog","vgs","vurglnb1","vurglnb2","vurglnb3","vurglnk1","vurglnk2","vurglnk3","vutglnb1","vutglnb2","vutglnb3",
	"vutglnk1","vutglnk2","vutglnk3"};
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
