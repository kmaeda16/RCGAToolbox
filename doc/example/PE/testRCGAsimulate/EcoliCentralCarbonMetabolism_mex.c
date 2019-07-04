#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "EcoliCentralCarbonMetabolism_mex.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

static double Function_for_vD_CS(double CS,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*CS/default0;
}

static double Function_for_vD_Fum(double Fum,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Fum/default0;
}

static double Function_for_vE_TktA(double GAP,double KTktA_eq,double R5P,double S7P,double X5P,double default0,double vTktA_max)
{
    return vTktA_max*(R5P*X5P-S7P*GAP/KTktA_eq)/default0;
}

static double Function_for_vBM_G6P(double G6P,double alphaACE,double alphaGLC,double default0,double kBM_ACE_G6P,double kBM_GLC_G6P)
{
    return (alphaGLC*kBM_GLC_G6P+alphaACE*kBM_ACE_G6P)*G6P/default0;
}

static double Function_for_vD_F6P(double F6P,double default0,double mu)
{
    return mu*F6P/default0;
}

static double Function_for_vD_E4P(double E4P,double default0,double mu)
{
    return mu*E4P/default0;
}

static double Function_for_vE_CS(double AcCoA,double CS,double KCS_AcCoA,double KCS_OAA,double KCS_OAA_AcCoA,double KCS_aKG,double OAA,double aKG,double default0,double kCS_cat)
{
    return CS*kCS_cat*OAA*AcCoA/((1.0+aKG/KCS_aKG)*KCS_OAA_AcCoA*KCS_AcCoA+KCS_AcCoA*OAA+(1.0+aKG/KCS_aKG)*KCS_OAA*AcCoA+OAA*AcCoA)/default0;
}

static double Function_for_vE_Pta(double AcCoA,double AcP,double CoA,double KPta_AcCoA_i,double KPta_AcP_i,double KPta_AcP_m,double KPta_CoA_i,double KPta_Pi_i,double KPta_Pi_m,double KPta_eq,double Pi,double default0,double vPta_max)
{
    return vPta_max*(1.0/(KPta_AcCoA_i*KPta_Pi_m))*(AcCoA*Pi-AcP*CoA/KPta_eq)/(1.0+AcCoA/KPta_AcCoA_i+Pi/KPta_Pi_i+AcP/KPta_AcP_i+CoA/KPta_CoA_i+AcCoA*Pi/(KPta_AcCoA_i*KPta_Pi_m)+AcP*CoA/(KPta_AcP_m*KPta_CoA_i))/default0;
}

static double Function_for_vG_maeB(double SS_Mez,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*SS_Mez/default0;
}

static double Function_for_vG_gltA(double CrpcAMP,double KgltA_Crp,double default0,double kexpr,double mu,double ngltA,double vgltA_Crp_bound,double vgltA_Crp_unbound)
{
    return mu*kexpr*((1.0-pow(CrpcAMP,ngltA)/(pow(CrpcAMP,ngltA)+pow(KgltA_Crp,ngltA)))*vgltA_Crp_unbound+pow(CrpcAMP,ngltA)/(pow(CrpcAMP,ngltA)+pow(KgltA_Crp,ngltA))*vgltA_Crp_bound)/default0;
}

static double Function_for_vE_Fba(double FBP,double Fba,double GAP,double KFba_DHAP,double KFba_FBP,double KFba_GAP,double KFba_GAP_inh,double KFba_eq,double VFba_blf,double default0,double kFba_cat)
{
    return Fba*kFba_cat*(FBP-pow(GAP,2.0)/KFba_eq)/(KFba_FBP+FBP+KFba_GAP*GAP/(KFba_eq*VFba_blf)+KFba_DHAP*GAP/(KFba_eq*VFba_blf)+FBP*GAP/KFba_GAP_inh+pow(GAP,2.0)/(KFba_eq*VFba_blf))/default0;
}

static double Function_for_vG_mdh(double CrpcAMP,double Kmdh_Crp,double default0,double kexpr,double mu,double vmdh_Crp_bound,double vmdh_Crp_unbound)
{
    return mu*kexpr*((1.0-CrpcAMP/(CrpcAMP+Kmdh_Crp))*vmdh_Crp_unbound+CrpcAMP/(CrpcAMP+Kmdh_Crp)*vmdh_Crp_bound)/default0;
}

static double Function_for_vG_icdA(double Cra,double KicdA_Cra,double default0,double kexpr,double mu,double vicdA_Cra_bound,double vicdA_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KicdA_Cra))*vicdA_Cra_unbound+Cra/(Cra+KicdA_Cra)*vicdA_Cra_bound)/default0;
}

static double Function_for_vE_SDH(double FUM,double KSDH_SUC_m,double KSDH_eq,double SUC,double default0,double vSDH1_max,double vSDH2_max)
{
    return vSDH1_max*vSDH2_max*(SUC-FUM/KSDH_eq)/(KSDH_SUC_m*vSDH2_max+vSDH2_max*SUC+vSDH1_max*FUM/KSDH_eq)/default0;
}

static double Function_for_vE_Pgl(double H,double KPgl_6PGL_m,double KPgl_6PG_m,double KPgl_eq,double KPgl_h1,double KPgl_h2,double default0,double sixPG,double sixPGL,double vPgl_max)
{
    return vPgl_max*(sixPGL-sixPG/KPgl_eq)/((1.0+H/KPgl_h1+KPgl_h2/H)*(KPgl_6PGL_m+sixPGL+KPgl_6PGL_m/KPgl_6PG_m*sixPG))/default0;
}

static double Function_for_vE_Tal(double E4P,double F6P,double GAP,double KTal_eq,double S7P,double default0,double vTal_max)
{
    return vTal_max*(GAP*S7P-E4P*F6P/KTal_eq)/default0;
}

static double Function_for_vE_Pgi(double F6P,double G6P,double KPgi_F6P,double KPgi_F6P_6pginh,double KPgi_G6P,double KPgi_G6P_6pginh,double KPgi_eq,double default0,double sixPG,double vPgi_max)
{
    return vPgi_max*(G6P-F6P/KPgi_eq)/(KPgi_G6P*(1.0+F6P/(KPgi_F6P*(1.0+sixPG/KPgi_F6P_6pginh))+sixPG/KPgi_G6P_6pginh)+G6P)/default0;
}

static double Function_for_vD_FUM(double FUM,double default0,double mu)
{
    return mu*FUM/default0;
}

static double Function_for_vD_G6P(double G6P,double default0,double mu)
{
    return mu*G6P/default0;
}

static double Function_for_vD_Fbp(double Fbp,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Fbp/default0;
}

static double Function_for_vD_GLC(double GLC,double default0,double mu)
{
    return mu*GLC/default0;
}

static double Function_for_vE_R5PI(double KR5PI_eq,double R5P,double RU5P,double default0,double vR5PI_max)
{
    return vR5PI_max*(RU5P-R5P/KR5PI_eq)/default0;
}

static double Function_for_vG_aceK(double Cra,double CrpcAMP,double Factor_aceK,double GOX,double IclRtotal,double KaceBAK_Cra,double KaceBAK_Crp,double KaceBAK_DNA,double KaceBAK_GOX,double KaceBAK_PYR,double KaceBAK_PYRprime,double LaceBAK,double PYR,double aceBAK_DNA,double default0,double kaceBAK_cat_IclR,double kexpr,double mu,double vaceBAK_Cra_bound,double vaceBAK_Cra_unbound,double vaceBAK_Crp_bound,double vaceBAK_Crp_unbound)
{
    return Factor_aceK*mu*kexpr*((1.0-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1.0-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1.0-aceBAK_DNA/KaceBAK_DNA*(1.0+PYR/KaceBAK_PYRprime)/(1.0+1.0/LaceBAK*(GOX/KaceBAK_GOX)*(1.0+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)/default0;
}

static double Function_for_vG_glk(double Cra,double Kglk_Cra,double default0,double kexpr,double mu,double vglk_Cra_bound,double vglk_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+Kglk_Cra))*vglk_Cra_unbound+Cra/(Cra+Kglk_Cra)*vglk_Cra_bound)/default0;
}

static double Function_for_vD_GAPDH(double GAPDH,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*GAPDH/default0;
}

static double Function_for_vD_GLCex(double D,double GLCex,double default0)
{
    return D*GLCex/default0;
}

static double Function_for_vD_GAP(double GAP,double default0,double mu)
{
    return mu*GAP/default0;
}

static double Function_for_vD_Fba(double Fba,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Fba/default0;
}

static double Function_for_vD_FBP(double FBP,double default0,double mu)
{
    return mu*FBP/default0;
}

static double Function_for_vE_Mez(double AcCoA,double KMez_AcCoA,double KMez_MAL,double KMez_cAMP,double LMez,double MAL,double Mez,double cAMP,double default0,double kMez_cat,double nMez)
{
    return Mez*kMez_cat*MAL/KMez_MAL*pow(1.0+MAL/KMez_MAL,nMez-1.0)/(pow(1.0+MAL/KMez_MAL,nMez)+LMez*pow(1.0+AcCoA/KMez_AcCoA+cAMP/KMez_cAMP,nMez))/default0;
}

static double Function_for_vD_Ppc(double Ppc,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Ppc/default0;
}

static double Function_for_vD_aKGDH(double aKGDH,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*aKGDH/default0;
}

static double Function_for_vE_Ppc(double FBP,double KPpc_FBP,double KPpc_PEP,double LPpc,double PEP,double Ppc,double default0,double kPpc_cat,double nPpc)
{
    return Ppc*kPpc_cat*PEP/KPpc_PEP*pow(1.0+PEP/KPpc_PEP,nPpc-1.0)/(pow(1.0+PEP/KPpc_PEP,nPpc)+LPpc/pow(1.0+FBP/KPpc_FBP,nPpc))/default0;
}

static double Function_for_vE_G6PDH(double G6P,double KG6PDH_G6P,double KG6PDH_NADP,double KG6PDH_NADPH_g6pinh,double KG6PDH_NADPH_nadpinh,double NADP,double NADPH,double default0,double vG6PDH_max)
{
    return vG6PDH_max*G6P*NADP/((G6P+KG6PDH_G6P)*(1.0+NADPH/KG6PDH_NADPH_g6pinh)*(KG6PDH_NADP*(1.0+NADPH/KG6PDH_NADPH_nadpinh)+NADP))/default0;
}

static double Function_for_vE_Glk(double ATP,double G6P,double GLC,double Glk,double KGlk_ATP_m,double KGlk_G6P_i,double KGlk_GLC_m,double default0,double kGlk_cat)
{
    return Glk*kGlk_cat*(GLC/KGlk_GLC_m)*(ATP/(KGlk_ATP_m*(1.0+G6P/KGlk_G6P_i)))/(1.0+GLC/KGlk_GLC_m+ATP/(KGlk_ATP_m*(1.0+G6P/KGlk_G6P_i))+GLC*ATP/(KGlk_GLC_m*KGlk_ATP_m*(1.0+G6P/KGlk_G6P_i))+G6P/KGlk_G6P_i)/default0;
}

static double Function_for_vE_MDH(double KMDH_MAL_I,double KMDH_MAL_m,double KMDH_NADH_I,double KMDH_NADH_m,double KMDH_NAD_I,double KMDH_NAD_II,double KMDH_NAD_m,double KMDH_OAA_I,double KMDH_OAA_II,double KMDH_OAA_m,double KMDH_eq,double MAL,double NAD,double NADH,double OAA,double default0,double vMDH1_max,double vMDH2_max)
{
    return vMDH1_max*vMDH2_max*(NAD*MAL-NADH*OAA/KMDH_eq)/(KMDH_NAD_I*KMDH_MAL_m*vMDH2_max+KMDH_MAL_m*vMDH2_max*NAD+KMDH_NAD_m*vMDH2_max*MAL+vMDH2_max*NAD*MAL+KMDH_OAA_m*vMDH1_max*NADH/KMDH_eq+KMDH_NADH_m*vMDH1_max*OAA/KMDH_eq+vMDH1_max*NADH*OAA/KMDH_eq+vMDH1_max*KMDH_OAA_m*NAD*NADH/(KMDH_eq*KMDH_NAD_I)+vMDH2_max*KMDH_NAD_m*MAL*OAA/KMDH_OAA_I+vMDH2_max*NAD*MAL*NADH/KMDH_NADH_I+vMDH1_max*MAL*NADH*OAA/(KMDH_eq*KMDH_MAL_I)+vMDH2_max*NAD*MAL*OAA/KMDH_OAA_II+vMDH1_max*NAD*NADH*OAA/(KMDH_NAD_II*KMDH_eq)+KMDH_NAD_I*vMDH2_max*NAD*MAL*NADH*OAA/(KMDH_NAD_II*KMDH_OAA_m*KMDH_NADH_I))/default0;
}

static double Function_for_vD_RU5P(double RU5P,double default0,double mu)
{
    return mu*RU5P/default0;
}

static double Function_for_vE_GAPDH(double GAP,double GAPDH,double KGAPDH_GAP,double KGAPDH_NAD,double KGAPDH_NADH,double KGAPDH_PGP,double KGAPDH_eq,double NAD,double NADH,double PEP,double default0,double kGAPDH_cat)
{
    return GAPDH*kGAPDH_cat*(GAP*NAD-PEP*NADH/KGAPDH_eq)/((KGAPDH_GAP*(1.0+PEP/KGAPDH_PGP)+GAP)*(KGAPDH_NAD*(1.0+NADH/KGAPDH_NADH)+NAD))/default0;
}

static double Function_for_vD_OAA(double OAA,double default0,double mu)
{
    return mu*OAA/default0;
}

static double Function_for_vE_Fbp(double FBP,double Fbp,double KFbp_FBP,double KFbp_PEP,double LFbp,double PEP,double default0,double kFbp_cat,double nFbp)
{
    return Fbp*kFbp_cat*FBP/KFbp_FBP*pow(1.0+FBP/KFbp_FBP,nFbp-1.0)/(pow(1.0+FBP/KFbp_FBP,nFbp)+LFbp/pow(1.0+PEP/KFbp_PEP,nFbp))/default0;
}

static double Function_for_vE_MS(double AcCoA,double GOX,double KMS_AcCoA,double KMS_GOX,double KMS_GOX_AcCoA,double MS,double default0,double kMS_cat)
{
    return MS*kMS_cat*GOX*AcCoA/(KMS_GOX_AcCoA*KMS_AcCoA+KMS_AcCoA*GOX+KMS_GOX*AcCoA+GOX*AcCoA)/default0;
}

static double Function_for_vD_R5P(double R5P,double default0,double mu)
{
    return mu*R5P/default0;
}

static double Function_for_vD_Pfk(double Pfk,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Pfk/default0;
}

static double Function_for_vD_PYR(double PYR,double default0,double mu)
{
    return mu*PYR/default0;
}

static double Function_for_vD_X(double D,double X,double default0)
{
    return D*X/default0;
}

static double Function_for_vE_Icl(double GAP,double ICIT,double Icl,double KIcl_3PG,double KIcl_ICIT,double KIcl_PEP,double KIcl_aKG,double LIcl,double PEP,double aKG,double default0,double kIcl_cat,double nIcl)
{
    return Icl*kIcl_cat*ICIT/KIcl_ICIT*pow(1.0+ICIT/KIcl_ICIT,nIcl-1.0)/(pow(1.0+ICIT/KIcl_ICIT,nIcl)+LIcl*pow(1.0+PEP/KIcl_PEP+GAP/KIcl_3PG+aKG/KIcl_aKG,nIcl))/default0;
}

static double Function_for_vE_Pps(double KPps_PEP,double KPps_PYR,double LPps,double PEP,double PYR,double Pps,double default0,double kPps_cat,double nPps)
{
    return Pps*kPps_cat*PYR/KPps_PYR*pow(1.0+PYR/KPps_PYR,nPps-1.0)/(pow(1.0+PYR/KPps_PYR,nPps)+LPps*pow(1.0+PEP/KPps_PEP,nPps))/default0;
}

static double Function_for_vD_SDH(double SDH,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*SDH/default0;
}

static double Function_for_vD_Pck(double Pck,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Pck/default0;
}

static double Function_for_vD_PEP(double PEP,double default0,double mu)
{
    return mu*PEP/default0;
}

static double Function_for_vD_PDH(double PDH,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*PDH/default0;
}

static double Function_for_vD_Pyk(double Pyk,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Pyk/default0;
}

static double Function_for_vE_Ack_medium(double ACEex,double ADP,double ATP,double AcP,double KAck_ACE_m,double KAck_ADP_m,double KAck_ATP_m,double KAck_AcP_m,double KAck_eq,double X,double default0,double rho,double vAck_max)
{
    return vAck_max*(1.0/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1.0+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1.0+ADP/KAck_ADP_m+ATP/KAck_ATP_m))*X/rho/default0;
}

static double Function_for_vE_Acs(double ACEex,double Acs,double KAcs_ACE,double default0,double kAcs_cat)
{
    return Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE)/default0;
}

static double Function_for_vE_AceKki(double AceK,double GAP,double GOX,double ICDH,double ICIT,double KAceK_3PG,double KAceK_GOX,double KAceK_ICDH,double KAceK_ICIT,double KAceK_OAA,double KAceK_PEP,double KAceK_PYR,double KAceK_aKG,double LAceK,double OAA,double PEP,double PYR,double aKG,double default0,double kAceKki_cat,double nAceK)
{
    return AceK*kAceKki_cat*ICDH/KAceK_ICDH*pow(1.0+ICDH/KAceK_ICDH,nAceK-1.0)/(pow(1.0+ICDH/KAceK_ICDH,nAceK)+LAceK*pow(1.0+ICIT/KAceK_ICIT+GOX/KAceK_GOX+OAA/KAceK_OAA+aKG/KAceK_aKG+PEP/KAceK_PEP+GAP/KAceK_3PG+PYR/KAceK_PYR,nAceK))/default0;
}

static double Function_for_vD_SUC(double SUC,double default0,double mu)
{
    return mu*SUC/default0;
}

static double Function_for_vD_Pps(double Pps,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Pps/default0;
}

static double Function_for_vD_S7P(double S7P,double default0,double mu)
{
    return mu*S7P/default0;
}

static double Function_for_vD_aKG(double aKG,double default0,double mu)
{
    return mu*aKG/default0;
}

static double Function_for_vE_Acs_medium(double ACEex,double Acs,double KAcs_ACE,double X,double default0,double kAcs_cat,double rho)
{
    return Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE)*X/rho/default0;
}

static double Function_for_vE_6PGDH(double ATP,double K6PGDH_6PG,double K6PGDH_ATP_inh,double K6PGDH_NADP,double K6PGDH_NADPH_inh,double NADP,double NADPH,double default0,double sixPG,double v6PGDH_max)
{
    return v6PGDH_max*sixPG*NADP/((sixPG+K6PGDH_6PG)*(NADP+K6PGDH_NADP*(1.0+NADPH/K6PGDH_NADPH_inh)*(1.0+ATP/K6PGDH_ATP_inh)))/default0;
}

static double Function_for_vE_Ack(double ACEex,double ADP,double ATP,double AcP,double KAck_ACE_m,double KAck_ADP_m,double KAck_ATP_m,double KAck_AcP_m,double KAck_eq,double default0,double vAck_max)
{
    return vAck_max*(1.0/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1.0+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1.0+ADP/KAck_ADP_m+ATP/KAck_ATP_m))/default0;
}

static double Function_for_vE_Cya(double EIIAP,double KCya_EIIAP,double default0,double vCya_max)
{
    return vCya_max*EIIAP/(EIIAP+KCya_EIIAP)/default0;
}

static double Function_for_vG_aceA(double Cra,double CrpcAMP,double GOX,double IclRtotal,double KaceBAK_Cra,double KaceBAK_Crp,double KaceBAK_DNA,double KaceBAK_GOX,double KaceBAK_PYR,double KaceBAK_PYRprime,double LaceBAK,double PYR,double aceBAK_DNA,double default0,double kaceBAK_cat_IclR,double kexpr,double mu,double vaceBAK_Cra_bound,double vaceBAK_Cra_unbound,double vaceBAK_Crp_bound,double vaceBAK_Crp_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1.0-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1.0-aceBAK_DNA/KaceBAK_DNA*(1.0+PYR/KaceBAK_PYRprime)/(1.0+1.0/LaceBAK*(GOX/KaceBAK_GOX)*(1.0+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)/default0;
}

static double Function_for_vG_fbaA(double Cra,double CrpcAMP,double KfbaA_Cra,double KfbaA_Crp,double default0,double kexpr,double mu,double vfbaA_Cra_bound,double vfbaA_Cra_unbound,double vfbaA_Crp_bound,double vfbaA_Crp_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KfbaA_Cra))*vfbaA_Cra_unbound+Cra/(Cra+KfbaA_Cra)*vfbaA_Cra_bound+(1.0-CrpcAMP/(CrpcAMP+KfbaA_Crp))*vfbaA_Crp_unbound+CrpcAMP/(CrpcAMP+KfbaA_Crp)*vfbaA_Crp_bound)/default0;
}

static double Function_for_vG_acs(double CrpcAMP,double Kacs_Crp,double default0,double kexpr,double mu,double nacs,double vacs_Crp_bound,double vacs_Crp_unbound)
{
    return mu*kexpr*((1.0-pow(CrpcAMP,nacs)/(pow(CrpcAMP,nacs)+pow(Kacs_Crp,nacs)))*vacs_Crp_unbound+pow(CrpcAMP,nacs)/(pow(CrpcAMP,nacs)+pow(Kacs_Crp,nacs))*vacs_Crp_bound)/default0;
}

static double Function_for_vD_X5P(double X5P,double default0,double mu)
{
    return mu*X5P/default0;
}

static double Function_for_vD_GLCfeed(double D,double GLCfeed,double default0)
{
    return D*GLCfeed/default0;
}

static double Function_for_vE_Edd(double KDPG,double KEdd_6PG_m,double KEdd_KDPG_m,double KEdd_eq,double QEdd_pH,double default0,double sixPG,double vEdd_max)
{
    return vEdd_max*QEdd_pH*(sixPG-KDPG/KEdd_eq)/(KEdd_6PG_m+sixPG+KEdd_6PG_m*KDPG/KEdd_KDPG_m)/default0;
}

static double Function_for_vG_fbp(double Cra,double Kfbp_Cra,double default0,double kexpr,double mu,double vfbp_Cra_bound,double vfbp_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+Kfbp_Cra))*vfbp_Cra_unbound+Cra/(Cra+Kfbp_Cra)*vfbp_Cra_bound)/default0;
}

static double Function_for_vG_fumABC(double CrpcAMP,double KfumABC_Crp,double default0,double kexpr,double mu,double nfumABC,double vfumABC_Crp_bound,double vfumABC_Crp_unbound)
{
    return mu*kexpr*((1.0-pow(CrpcAMP,nfumABC)/(pow(CrpcAMP,nfumABC)+pow(KfumABC_Crp,nfumABC)))*vfumABC_Crp_unbound+pow(CrpcAMP,nfumABC)/(pow(CrpcAMP,nfumABC)+pow(KfumABC_Crp,nfumABC))*vfumABC_Crp_bound)/default0;
}

static double Function_for_vE_Eda(double GAP,double KDPG,double KEda_GAP_m,double KEda_KDPG_m,double KEda_PYR_m,double KEda_eq,double PYR,double QEda_pH,double default0,double vEda_max)
{
    return vEda_max*QEda_pH*(KDPG-GAP*PYR/KEda_eq)/(KEda_KDPG_m+KDPG+KEda_KDPG_m*(PYR/KEda_PYR_m+GAP/KEda_GAP_m+PYR*GAP/(KEda_PYR_m*KEda_GAP_m)))/default0;
}

static double Function_for_vG_gapA(double Cra,double CrpcAMP,double KgapA_Cra,double KgapA_Crp,double default0,double kexpr,double mu,double vgapA_Cra_bound,double vgapA_Cra_unbound,double vgapA_Crp_bound,double vgapA_Crp_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KgapA_Cra))*vgapA_Cra_unbound+Cra/(Cra+KgapA_Cra)*vgapA_Cra_bound+(1.0-CrpcAMP/(CrpcAMP+KgapA_Crp))*vgapA_Crp_unbound+CrpcAMP/(CrpcAMP+KgapA_Crp)*vgapA_Crp_bound)/default0;
}

static double Function_for_vE_AceKph(double AceK,double GAP,double ICDHP,double KAceK_3PG,double KAceK_ICDHP,double KAceK_OAA,double KAceK_PEP,double KAceK_PYR,double KAceK_aKG,double LAceK,double OAA,double PEP,double PYR,double aKG,double default0,double kAceKph_cat,double nAceK)
{
    return AceK*kAceKph_cat*ICDHP/KAceK_ICDHP*pow(1.0+ICDHP/KAceK_ICDHP,nAceK-1.0)/(pow(1.0+ICDHP/KAceK_ICDHP,nAceK)+LAceK/pow(1.0+OAA/KAceK_OAA+aKG/KAceK_aKG+PEP/KAceK_PEP+GAP/KAceK_3PG+PYR/KAceK_PYR,nAceK))/default0;
}

static double Function_for_vD_cAMP(double cAMP,double default0,double mu)
{
    return mu*cAMP/default0;
}

static double Function_for_vG_aceB(double Cra,double CrpcAMP,double Factor_aceB,double GOX,double IclRtotal,double KaceBAK_Cra,double KaceBAK_Crp,double KaceBAK_DNA,double KaceBAK_GOX,double KaceBAK_PYR,double KaceBAK_PYRprime,double LaceBAK,double PYR,double aceBAK_DNA,double default0,double kaceBAK_cat_IclR,double kexpr,double mu,double vaceBAK_Cra_bound,double vaceBAK_Cra_unbound,double vaceBAK_Crp_bound,double vaceBAK_Crp_unbound)
{
    return Factor_aceB*mu*kexpr*((1.0-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1.0-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1.0-aceBAK_DNA/KaceBAK_DNA*(1.0+PYR/KaceBAK_PYRprime)/(1.0+1.0/LaceBAK*(GOX/KaceBAK_GOX)*(1.0+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)/default0;
}

static double Function_for_vPTS1(double EIIA,double EIIAP,double PEP,double PYR,double default0,double kPTS1,double kmPTS1)
{
    return (kPTS1*PEP*EIIA-kmPTS1*PYR*EIIAP)/default0;
}

static double Function_for_vE_Pyk(double ADP,double AMP,double ATP,double FBP,double KPyk_ADP,double KPyk_AMP,double KPyk_ATP,double KPyk_FBP,double KPyk_PEP,double LPyk,double PEP,double Pyk,double default0,double kPyk_cat,double nPyk)
{
    return Pyk*kPyk_cat*PEP*pow(PEP/KPyk_PEP+1.0,nPyk-1.0)*ADP/(KPyk_PEP*(LPyk*pow((1.0+ATP/KPyk_ATP)/(FBP/KPyk_FBP+AMP/KPyk_AMP+1.0),nPyk)+pow(PEP/KPyk_PEP+1.0,nPyk))*(ADP+KPyk_ADP))/default0;
}

static double Function_for_vE_Ru5P(double KRu5P_eq,double RU5P,double X5P,double default0,double vRu5P_max)
{
    return vRu5P_max*(RU5P-X5P/KRu5P_eq)/default0;
}

static double Function_for_vE_Pfk(double A,double ADP,double ATP,double B,double F6P,double KPfk_ADP_c,double KPfk_ATP_s,double KPfk_F6P_s,double LPfk,double Pfk,double default0,double kPfk_cat,double nPfk)
{
    return Pfk*kPfk_cat*ATP*F6P/((ATP+KPfk_ATP_s*(1.0+ADP/KPfk_ADP_c))*(F6P+KPfk_F6P_s*A/B)*(1.0+LPfk/pow(1.0+F6P*B/(KPfk_F6P_s*A),nPfk)))/default0;
}

static double Function_for_vG_sdhCDAB(double CrpcAMP,double KsdhCDAB_Crp,double default0,double kexpr,double mu,double nsdhCDAB,double vsdhCDAB_Crp_bound,double vsdhCDAB_Crp_unbound)
{
    return mu*kexpr*((1.0-pow(CrpcAMP,nsdhCDAB)/(pow(CrpcAMP,nsdhCDAB)+pow(KsdhCDAB_Crp,nsdhCDAB)))*vsdhCDAB_Crp_unbound+pow(CrpcAMP,nsdhCDAB)/(pow(CrpcAMP,nsdhCDAB)+pow(KsdhCDAB_Crp,nsdhCDAB))*vsdhCDAB_Crp_bound)/default0;
}

static double Function_for_vBM_GAP(double GAP,double alphaACE,double alphaGLC,double default0,double kBM_ACE_GAP,double kBM_GLC_GAP)
{
    return (alphaGLC*kBM_GLC_GAP+alphaACE*kBM_ACE_GAP)*GAP/default0;
}

static double Function_for_vD_ICDHP(double ICDHP,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*ICDHP/default0;
}

static double Function_for_vD_MDH(double MDH,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*MDH/default0;
}

static double Function_for_vPTS4(double EIIAP,double GLCex,double KPTS_EIIA,double KPTS_GLC,double default0,double vPTS4_max)
{
    return vPTS4_max*EIIAP*GLCex/((KPTS_EIIA+EIIAP)*(KPTS_GLC+GLCex))/default0;
}

static double Function_for_vgrowth(double X,double default0,double mu)
{
    return mu*X/default0;
}

static double Function_for_vPTS4_medium(double EIIAP,double GLCex,double KPTS_EIIA,double KPTS_GLC,double X,double default0,double rho,double vPTS4_max)
{
    return vPTS4_max*EIIAP*GLCex/((KPTS_EIIA+EIIAP)*(KPTS_GLC+GLCex))*X/rho/default0;
}

static double Function_for_vNonPTS_medium(double EIIA,double GLCex,double KNonPTS_I,double KNonPTS_S,double X,double default0,double rho,double vNonPTS_max)
{
    return vNonPTS_max*GLCex/(KNonPTS_S+(1.0+EIIA/KNonPTS_I)*GLCex)*X/rho/default0;
}

static double Function_for_vNonPTS(double EIIA,double GLCex,double KNonPTS_I,double KNonPTS_S,double default0,double vNonPTS_max)
{
    return vNonPTS_max*GLCex/(KNonPTS_S+(1.0+EIIA/KNonPTS_I)*GLCex)/default0;
}

static double Function_for_vD_GOX(double GOX,double default0,double mu)
{
    return mu*GOX/default0;
}

static double Function_for_vE_cAMPdegr(double KcAMPdegr_cAMP,double cAMP,double default0,double vcAMPdegr_max)
{
    return vcAMPdegr_max*cAMP/(cAMP+KcAMPdegr_cAMP)/default0;
}

static double Function_for_vBM_E4P(double E4P,double alphaACE,double alphaGLC,double default0,double kBM_ACE_E4P,double kBM_GLC_E4P)
{
    return (alphaGLC*kBM_GLC_E4P+alphaACE*kBM_ACE_E4P)*E4P/default0;
}

static double Function_for_vE_TktB(double E4P,double F6P,double GAP,double KTktB_eq,double X5P,double default0,double vTktB_max)
{
    return vTktB_max*(X5P*E4P-F6P*GAP/KTktB_eq)/default0;
}

static double Function_for_vG_sucAB(double CrpcAMP,double KsucAB_Crp,double default0,double kexpr,double mu,double nsucAB,double vsucAB_Crp_bound,double vsucAB_Crp_unbound)
{
    return mu*kexpr*((1.0-pow(CrpcAMP,nsucAB)/(pow(CrpcAMP,nsucAB)+pow(KsucAB_Crp,nsucAB)))*vsucAB_Crp_unbound+pow(CrpcAMP,nsucAB)/(pow(CrpcAMP,nsucAB)+pow(KsucAB_Crp,nsucAB))*vsucAB_Crp_bound)/default0;
}

static double Function_for_vE_Pck(double ADP,double ATP,double KPck_ADP_i,double KPck_ATP_I,double KPck_ATP_i,double KPck_OAA,double KPck_OAA_I,double KPck_PEP,double KPck_PEP_i,double OAA,double PEP,double Pck,double default0,double kPck_cat)
{
    return Pck*kPck_cat*OAA*ATP/ADP/(KPck_OAA*ATP/ADP+OAA*ATP/ADP+KPck_ATP_i*KPck_OAA/KPck_ADP_i+KPck_ATP_i*KPck_OAA/(KPck_PEP*KPck_ADP_i)*PEP+KPck_ATP_i*KPck_OAA/(KPck_PEP_i*KPck_ATP_I)*ATP/ADP*PEP+KPck_ATP_i*KPck_OAA/(KPck_ADP_i*KPck_OAA_I)*OAA)/default0;
}

static double Function_for_vG_pfkA(double Cra,double KpfkA_Cra,double default0,double kexpr,double mu,double vpfkA_Cra_bound,double vpfkA_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KpfkA_Cra))*vpfkA_Cra_unbound+Cra/(Cra+KpfkA_Cra)*vpfkA_Cra_bound)/default0;
}

static double Function_for_vG_pckA(double Cra,double KpckA_Cra,double default0,double kexpr,double mu,double vpckA_Cra_bound,double vpckA_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KpckA_Cra))*vpckA_Cra_unbound+Cra/(Cra+KpckA_Cra)*vpckA_Cra_bound)/default0;
}

static double Function_for_vE_aKGDH(double CoA,double KaKGDH_CoA_m,double KaKGDH_NADH_I,double KaKGDH_NAD_m,double KaKGDH_SUC_I,double KaKGDH_Z,double KaKGDH_aKG_I,double KaKGDH_aKG_m,double NAD,double NADH,double SUC,double aKG,double aKGDH,double default0,double kaKGDH_cat)
{
    return aKGDH*kaKGDH_cat*aKG*CoA*NAD/(KaKGDH_NAD_m*aKG*CoA+KaKGDH_CoA_m*aKG*NAD+KaKGDH_aKG_m*CoA*NAD+aKG*CoA*NAD+KaKGDH_aKG_m*KaKGDH_Z*SUC*NADH/KaKGDH_SUC_I+KaKGDH_NAD_m*aKG*CoA*NADH/KaKGDH_NADH_I+KaKGDH_CoA_m*aKG*NAD*SUC/KaKGDH_SUC_I+KaKGDH_aKG_m*KaKGDH_Z*aKG*SUC*NADH/(KaKGDH_aKG_I*KaKGDH_SUC_I))/default0;
}

static double Function_for_vD_ICDH(double ICDH,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*ICDH/default0;
}

static double Function_for_vD_MS(double MS,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*MS/default0;
}

static double Function_for_vD_Icl(double Icl,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Icl/default0;
}

static double Function_for_vG_pdh(double Kpdh_PdhR,double PdhR,double default0,double kexpr,double mu,double vpdh_PdhR_bound,double vpdh_PdhR_unbound)
{
    return mu*kexpr*((1.0-PdhR/(PdhR+Kpdh_PdhR))*vpdh_PdhR_unbound+PdhR/(PdhR+Kpdh_PdhR)*vpdh_PdhR_bound)/default0;
}

static double Function_for_vD_Mez(double Mez,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Mez/default0;
}

static double Function_for_vG_ppc(double SS_Ppc,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*SS_Ppc/default0;
}

static double Function_for_vD_Glk(double Glk,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Glk/default0;
}

static double Function_for_vG_ppsA(double Cra,double KppsA_Cra,double default0,double kexpr,double mu,double vppsA_Cra_bound,double vppsA_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KppsA_Cra))*vppsA_Cra_unbound+Cra/(Cra+KppsA_Cra)*vppsA_Cra_bound)/default0;
}

static double Function_for_vBM_AcCoA(double AcCoA,double alphaACE,double alphaGLC,double default0,double kBM_ACE_AcCoA,double kBM_GLC_AcCoA)
{
    return (alphaGLC*kBM_GLC_AcCoA+alphaACE*kBM_ACE_AcCoA)*AcCoA/default0;
}

static double Function_for_vD_KDPG(double KDPG,double default0,double mu)
{
    return mu*KDPG/default0;
}

static double Function_for_vG_pykF(double Cra,double KpykF_Cra,double default0,double kexpr,double mu,double vpykF_Cra_bound,double vpykF_Cra_unbound)
{
    return mu*kexpr*((1.0-Cra/(Cra+KpykF_Cra))*vpykF_Cra_unbound+Cra/(Cra+KpykF_Cra)*vpykF_Cra_bound)/default0;
}

static double Function_for_vBM_PEP(double PEP,double alphaACE,double alphaGLC,double default0,double kBM_ACE_PEP,double kBM_GLC_PEP)
{
    return (alphaGLC*kBM_GLC_PEP+alphaACE*kBM_ACE_PEP)*PEP/default0;
}

static double Function_for_vD_MAL(double MAL,double default0,double mu)
{
    return mu*MAL/default0;
}

static double Function_for_vD_6PG(double default0,double mu,double sixPG)
{
    return mu*sixPG/default0;
}

static double Function_for_vD_AcCoA(double AcCoA,double default0,double mu)
{
    return mu*AcCoA/default0;
}

static double Function_for_vD_AcP(double AcP,double default0,double mu)
{
    return mu*AcP/default0;
}

static double Function_for_vD_Acs(double Acs,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*Acs/default0;
}

static double Function_for_vD_6PGL(double default0,double mu,double sixPGL)
{
    return mu*sixPGL/default0;
}

static double Function_for_vBM_aKG(double aKG,double alphaACE,double alphaGLC,double default0,double kBM_ACE_aKG,double kBM_GLC_aKG)
{
    return (alphaGLC*kBM_GLC_aKG+alphaACE*kBM_ACE_aKG)*aKG/default0;
}

static double Function_for_vBM_OAA(double OAA,double alphaACE,double alphaGLC,double default0,double kBM_ACE_OAA,double kBM_GLC_OAA)
{
    return (alphaGLC*kBM_GLC_OAA+alphaACE*kBM_ACE_OAA)*OAA/default0;
}

static double Function_for_vD_ACEex(double ACEex,double D,double default0)
{
    return D*ACEex/default0;
}

static double Function_for_vD_ICIT(double ICIT,double default0,double mu)
{
    return mu*ICIT/default0;
}

static double Function_for_vBM_SUC(double SUC,double alphaACE,double alphaGLC,double default0,double kBM_ACE_SUC,double kBM_GLC_SUC)
{
    return (alphaGLC*kBM_GLC_SUC+alphaACE*kBM_ACE_SUC)*SUC/default0;
}

static double Function_for_vE_PDH(double AcCoA,double CoA,double KPDH_AcCoA_m,double KPDH_CoA_m,double KPDH_NADH_m,double KPDH_NAD_m,double KPDH_PYR_m,double KPDH_i,double NAD,double NADH,double PDH,double PYR,double default0,double kPDH_cat)
{
    return PDH*kPDH_cat*(1.0/(1.0+KPDH_i*NADH/NAD))*(PYR/KPDH_PYR_m)*(NAD/KPDH_NAD_m)*(CoA/KPDH_CoA_m)/((1.0+PYR/KPDH_PYR_m)*(1.0+NAD/KPDH_NAD_m+NADH/KPDH_NADH_m)*(1.0+CoA/KPDH_CoA_m+AcCoA/KPDH_AcCoA_m))/default0;
}

static double Function_for_vE_ICDH(double ICDH,double ICIT,double KICDH_ICIT,double KICDH_PEP,double LICDH,double PEP,double default0,double kICDH_cat,double nICDH)
{
    return ICDH*kICDH_cat*ICIT/KICDH_ICIT*pow(1.0+ICIT/KICDH_ICIT,nICDH-1.0)/(pow(1.0+ICIT/KICDH_ICIT,nICDH)+LICDH*pow(1.0+PEP/KICDH_PEP,nICDH))/default0;
}

static double Function_for_vBM_PYR(double PYR,double alphaACE,double alphaGLC,double default0,double kBM_ACE_PYR,double kBM_GLC_PYR)
{
    return (alphaGLC*kBM_GLC_PYR+alphaACE*kBM_ACE_PYR)*PYR/default0;
}

static double Function_for_vBM_R5P(double R5P,double alphaACE,double alphaGLC,double default0,double kBM_ACE_R5P,double kBM_GLC_R5P)
{
    return (alphaGLC*kBM_GLC_R5P+alphaACE*kBM_ACE_R5P)*R5P/default0;
}

static double Function_for_vE_Fum(double FUM,double KFum_FUM_m,double KFum_eq,double MAL,double default0,double vFum1_max,double vFum2_max)
{
    return vFum1_max*vFum2_max*(FUM-MAL/KFum_eq)/(KFum_FUM_m*vFum2_max+vFum2_max*FUM+vFum1_max*MAL/KFum_eq)/default0;
}

static double Function_for_vD_AceK(double AceK,double default0,double kdegr,double mu)
{
    return (mu+kdegr)*AceK/default0;
}

static double Function_for_vBM_F6P(double F6P,double alphaACE,double alphaGLC,double default0,double kBM_ACE_F6P,double kBM_GLC_F6P)
{
    return (alphaGLC*kBM_GLC_F6P+alphaACE*kBM_ACE_F6P)*F6P/default0;
}

static double Function_for_vBM_FUM(double FUM,double alphaACE,double alphaGLC,double default0,double kBM_ACE_FUM,double kBM_GLC_FUM)
{
    return (alphaGLC*kBM_GLC_FUM+alphaACE*kBM_ACE_FUM)*FUM/default0;
}

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double ACEex,AcCoA,AcP,AceK,Acs,CS,E4P,EIIAP,F6P,FBP,FUM,Fba,Fbp,Fum,G6P,GAP,GAPDH,GLC,GLCex,GLCfeed;
    double GOX,Glk,ICDH,ICDHP,ICIT,Icl,KDPG,MAL,MDH,MS,Mez,OAA,PDH,PEP,PYR,Pck,Pfk,Ppc,Pps,Pyk;
    double R5P,RU5P,S7P,SDH,SUC,X,X5P,aKG,aKGDH,cAMP,sixPG,sixPGL;
    double ADP,AMP,ATP,CoA,Cratotal,Crptotal,D,EIIAtotal,EXTERNAL,Factor_aceB,Factor_aceK,IclRtotal,K6PGDH_6PG,K6PGDH_ATP_inh,K6PGDH_NADP,K6PGDH_NADPH_inh,KAceK_3PG,KAceK_GOX,KAceK_ICDH,KAceK_ICDHP;
    double KAceK_ICIT,KAceK_OAA,KAceK_PEP,KAceK_PYR,KAceK_aKG,KAck_ACE_m,KAck_ADP_m,KAck_ATP_m,KAck_AcP_m,KAck_eq,KAcs_ACE,KCS_AcCoA,KCS_OAA,KCS_OAA_AcCoA,KCS_aKG,KCraFBP,KCrpcAMP,KCya_EIIAP,KEda_GAP_m,KEda_KDPG_m;
    double KEda_PYR_m,KEda_eq,KEdd_6PG_m,KEdd_KDPG_m,KEdd_eq,KFba_DHAP,KFba_FBP,KFba_GAP,KFba_GAP_inh,KFba_eq,KFbp_FBP,KFbp_PEP,KFum_FUM_m,KFum_eq,KG6PDH_G6P,KG6PDH_NADP,KG6PDH_NADPH_g6pinh,KG6PDH_NADPH_nadpinh,KGAPDH_GAP,KGAPDH_NAD;
    double KGAPDH_NADH,KGAPDH_PGP,KGAPDH_eq,KGlk_ATP_m,KGlk_G6P_i,KGlk_GLC_m,KICDH_ICIT,KICDH_PEP,KIcl_3PG,KIcl_ICIT,KIcl_PEP,KIcl_aKG,KMDH_MAL_I,KMDH_MAL_m,KMDH_NADH_I,KMDH_NADH_m,KMDH_NAD_I,KMDH_NAD_II,KMDH_NAD_m,KMDH_OAA_I;
    double KMDH_OAA_II,KMDH_OAA_m,KMDH_eq,KMS_AcCoA,KMS_GOX,KMS_GOX_AcCoA,KMez_AcCoA,KMez_MAL,KMez_cAMP,KNonPTS_I,KNonPTS_S,KPDH_AcCoA_m,KPDH_CoA_m,KPDH_NADH_m,KPDH_NAD_m,KPDH_PYR_m,KPDH_i,KPTS_EIIA,KPTS_GLC,KPck_ADP_i;
    double KPck_ATP_I,KPck_ATP_i,KPck_OAA,KPck_OAA_I,KPck_PEP,KPck_PEP_i,KPdhRPYR,KPfk_ADP_a,KPfk_ADP_b,KPfk_ADP_c,KPfk_AMP_a,KPfk_AMP_b,KPfk_ATP_s,KPfk_F6P_s,KPfk_PEP,KPgi_F6P,KPgi_F6P_6pginh,KPgi_G6P,KPgi_G6P_6pginh,KPgi_eq;
    double KPgl_6PGL_m,KPgl_6PG_m,KPgl_eq,KPgl_h1,KPgl_h2,KPpc_FBP,KPpc_PEP,KPps_PEP,KPps_PYR,KPta_AcCoA_i,KPta_AcP_i,KPta_AcP_m,KPta_CoA_i,KPta_Pi_i,KPta_Pi_m,KPta_eq,KPyk_ADP,KPyk_AMP,KPyk_ATP,KPyk_FBP;
    double KPyk_PEP,KR5PI_eq,KRu5P_eq,KSDH_SUC_m,KSDH_eq,KTal_eq,KTktA_eq,KTktB_eq,KaKGDH_CoA_m,KaKGDH_NADH_I,KaKGDH_NAD_m,KaKGDH_SUC_I,KaKGDH_Z,KaKGDH_aKG_I,KaKGDH_aKG_m,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR;
    double KaceBAK_PYRprime,Kacs_Crp,KcAMPdegr_cAMP,KfbaA_Cra,KfbaA_Crp,Kfbp_Cra,KfumABC_Crp,KgapA_Cra,KgapA_Crp,Kglk_Cra,KgltA_Crp,KicdA_Cra,Kmdh_Crp,KpckA_Cra,Kpdh_PdhR,KpfkA_Cra,KppsA_Cra,KpykF_Cra,KsdhCDAB_Crp,KsucAB_Crp;
    double LAceK,LFbp,LICDH,LIcl,LMez,LPfk,LPpc,LPps,LPyk,LaceBAK,NAD,NADH,NADP,NADPH,POratio,POratio_prime,PdhRtotal,Pi,SS_Mez_ACE,SS_Mez_GLC;
    double SS_Ppc_ACE,SS_Ppc_GLC,VFba_blf,aceBAK_DNA,kATP,kAceKki_cat,kAceKph_cat,kAcs_cat,kBM_ACE_AcCoA,kBM_ACE_E4P,kBM_ACE_F6P,kBM_ACE_FUM,kBM_ACE_G6P,kBM_ACE_GAP,kBM_ACE_OAA,kBM_ACE_PEP,kBM_ACE_PYR,kBM_ACE_R5P,kBM_ACE_SUC,kBM_ACE_aKG;
    double kBM_GLC_AcCoA,kBM_GLC_E4P,kBM_GLC_F6P,kBM_GLC_FUM,kBM_GLC_G6P,kBM_GLC_GAP,kBM_GLC_OAA,kBM_GLC_PEP,kBM_GLC_PYR,kBM_GLC_R5P,kBM_GLC_SUC,kBM_GLC_aKG,kCS_cat,kFba_cat,kFbp_cat,kFum1_cat,kFum2_cat,kGAPDH_cat,kGlk_cat,kICDH_cat;
    double kIcl_cat,kMDH1_cat,kMDH2_cat,kMS_cat,kMez_cat,kPDH_cat,kPTS1,kPck_cat,kPfk_cat,kPpc_cat,kPps_cat,kPyk_cat,kSDH1_cat,kSDH2_cat,kaKGDH_cat,kaceBAK_cat_IclR,kdegr,kexpr,kmPTS1,nAceK;
    double nCraFBP,nCrpcAMP,nFbp,nICDH,nIcl,nMez,nPdhRPYR,nPfk,nPpc,nPps,nPyk,nacs,nfumABC,ngltA,nsdhCDAB,nsucAB,pH,pH_Eda_m,pH_Edd_m,pK_Eda;
    double pK_Edd,rho,v6PGDH_max,vAck_max,vCya_max,vEda_max,vEdd_max,vG6PDH_max,vNonPTS_max,vPTS4_max,vPgi_max,vPgl_max,vPta_max,vR5PI_max,vRu5P_max,vTal_max,vTktA_max,vTktB_max,vaceBAK_Cra_bound,vaceBAK_Cra_unbound;
    double vaceBAK_Crp_bound,vaceBAK_Crp_unbound,vacs_Crp_bound,vacs_Crp_unbound,vcAMPdegr_max,vfbaA_Cra_bound,vfbaA_Cra_unbound,vfbaA_Crp_bound,vfbaA_Crp_unbound,vfbp_Cra_bound,vfbp_Cra_unbound,vfumABC_Crp_bound,vfumABC_Crp_unbound,vgapA_Cra_bound,vgapA_Cra_unbound,vgapA_Crp_bound,vgapA_Crp_unbound,vglk_Cra_bound,vglk_Cra_unbound,vgltA_Crp_bound;
    double vgltA_Crp_unbound,vicdA_Cra_bound,vicdA_Cra_unbound,vmdh_Crp_bound,vmdh_Crp_unbound,vpckA_Cra_bound,vpckA_Cra_unbound,vpdh_PdhR_bound,vpdh_PdhR_unbound,vpfkA_Cra_bound,vpfkA_Cra_unbound,vppsA_Cra_bound,vppsA_Cra_unbound,vpykF_Cra_bound,vpykF_Cra_unbound,vsdhCDAB_Crp_bound,vsdhCDAB_Crp_unbound,vsucAB_Crp_bound,vsucAB_Crp_unbound,default0;
    double CraFBP,B,A,Cra,CrpcAMP,Crp,H,EIIA,alphaGLC,alphaACE,SS_Ppc,SS_Mez,OP_NADH,OP_FADH2,PdhRPYR,PdhR,QEda_pH,QEdd_pH,vMDH1_max,vFum2_max;
    double vSDH2_max,vSDH1_max,vATP,mu,vMDH2_max,vFum1_max;
    double vBM_AcCoA,vBM_E4P,vBM_F6P,vBM_FUM,vBM_G6P,vBM_GAP,vBM_OAA,vBM_PEP,vBM_PYR,vBM_R5P,vBM_SUC,vBM_aKG,vD_6PG,vD_6PGL,vD_ACEex,vD_AcCoA,vD_AcP,vD_AceK,vD_Acs,vD_CS;
    double vD_E4P,vD_F6P,vD_FBP,vD_FUM,vD_Fba,vD_Fbp,vD_Fum,vD_G6P,vD_GAP,vD_GAPDH,vD_GLC,vD_GLCex,vD_GLCfeed,vD_GOX,vD_Glk,vD_ICDH,vD_ICDHP,vD_ICIT,vD_Icl,vD_KDPG;
    double vD_MAL,vD_MDH,vD_MS,vD_Mez,vD_OAA,vD_PDH,vD_PEP,vD_PYR,vD_Pck,vD_Pfk,vD_Ppc,vD_Pps,vD_Pyk,vD_R5P,vD_RU5P,vD_S7P,vD_SDH,vD_SUC,vD_X,vD_X5P;
    double vD_aKG,vD_aKGDH,vD_cAMP,vE_6PGDH,vE_AceKki,vE_AceKph,vE_Ack,vE_Ack_medium,vE_Acs,vE_Acs_medium,vE_CS,vE_Cya,vE_Eda,vE_Edd,vE_Fba,vE_Fbp,vE_Fum,vE_G6PDH,vE_GAPDH,vE_Glk;
    double vE_ICDH,vE_Icl,vE_MDH,vE_MS,vE_Mez,vE_PDH,vE_Pck,vE_Pfk,vE_Pgi,vE_Pgl,vE_Ppc,vE_Pps,vE_Pta,vE_Pyk,vE_R5PI,vE_Ru5P,vE_SDH,vE_Tal,vE_TktA,vE_TktB;
    double vE_aKGDH,vE_cAMPdegr,vG_aceA,vG_aceB,vG_aceK,vG_acs,vG_fbaA,vG_fbp,vG_fumABC,vG_gapA,vG_glk,vG_gltA,vG_icdA,vG_maeB,vG_mdh,vG_pckA,vG_pdh,vG_pfkA,vG_ppc,vG_ppsA;
    double vG_pykF,vG_sdhCDAB,vG_sucAB,vNonPTS,vNonPTS_medium,vPTS1,vPTS4,vPTS4_medium,vgrowth;

    time = time_local;

    ACEex = stateVector[0];
    AcCoA = stateVector[1];
    AcP = stateVector[2];
    AceK = stateVector[3];
    Acs = stateVector[4];
    CS = stateVector[5];
    E4P = stateVector[6];
    EIIAP = stateVector[7];
    F6P = stateVector[8];
    FBP = stateVector[9];
    FUM = stateVector[10];
    Fba = stateVector[11];
    Fbp = stateVector[12];
    Fum = stateVector[13];
    G6P = stateVector[14];
    GAP = stateVector[15];
    GAPDH = stateVector[16];
    GLC = stateVector[17];
    GLCex = stateVector[18];
    GLCfeed = stateVector[19];
    GOX = stateVector[20];
    Glk = stateVector[21];
    ICDH = stateVector[22];
    ICDHP = stateVector[23];
    ICIT = stateVector[24];
    Icl = stateVector[25];
    KDPG = stateVector[26];
    MAL = stateVector[27];
    MDH = stateVector[28];
    MS = stateVector[29];
    Mez = stateVector[30];
    OAA = stateVector[31];
    PDH = stateVector[32];
    PEP = stateVector[33];
    PYR = stateVector[34];
    Pck = stateVector[35];
    Pfk = stateVector[36];
    Ppc = stateVector[37];
    Pps = stateVector[38];
    Pyk = stateVector[39];
    R5P = stateVector[40];
    RU5P = stateVector[41];
    S7P = stateVector[42];
    SDH = stateVector[43];
    SUC = stateVector[44];
    X = stateVector[45];
    X5P = stateVector[46];
    aKG = stateVector[47];
    aKGDH = stateVector[48];
    cAMP = stateVector[49];
    sixPG = stateVector[50];
    sixPGL = stateVector[51];
    ADP = paramdataPtr->parametervector[0]; /* 0.56 */
    AMP = paramdataPtr->parametervector[1]; /* 0.28 */
    ATP = paramdataPtr->parametervector[2]; /* 9.6 */
    CoA = paramdataPtr->parametervector[3]; /* 1.4 */
    Cratotal = paramdataPtr->parametervector[4]; /* 0.0003 */
    Crptotal = paramdataPtr->parametervector[5]; /* 0.0115 */
    D = paramdataPtr->parametervector[6]; /* 0 */
    EIIAtotal = paramdataPtr->parametervector[7]; /* 0.0769 */
    EXTERNAL = paramdataPtr->parametervector[8]; /* 0 */
    Factor_aceB = paramdataPtr->parametervector[9]; /* 0.509128 */
    Factor_aceK = paramdataPtr->parametervector[10]; /* 0.0293109 */
    IclRtotal = paramdataPtr->parametervector[11]; /* 8.3e-05 */
    K6PGDH_6PG = paramdataPtr->parametervector[12]; /* 0.0999956 */
    K6PGDH_ATP_inh = paramdataPtr->parametervector[13]; /* 3.03764 */
    K6PGDH_NADP = paramdataPtr->parametervector[14]; /* 0.00774325 */
    K6PGDH_NADPH_inh = paramdataPtr->parametervector[15]; /* 0.0180077 */
    KAceK_3PG = paramdataPtr->parametervector[16]; /* 1.0927 */
    KAceK_GOX = paramdataPtr->parametervector[17]; /* 0.415957 */
    KAceK_ICDH = paramdataPtr->parametervector[18]; /* 0.731925 */
    KAceK_ICDHP = paramdataPtr->parametervector[19]; /* 8.19069 */
    KAceK_ICIT = paramdataPtr->parametervector[20]; /* 0.0658771 */
    KAceK_OAA = paramdataPtr->parametervector[21]; /* 0.0345864 */
    KAceK_PEP = paramdataPtr->parametervector[22]; /* 0.232436 */
    KAceK_PYR = paramdataPtr->parametervector[23]; /* 0.0311645 */
    KAceK_aKG = paramdataPtr->parametervector[24]; /* 0.723859 */
    KAck_ACE_m = paramdataPtr->parametervector[25]; /* 6.33002 */
    KAck_ADP_m = paramdataPtr->parametervector[26]; /* 0.296006 */
    KAck_ATP_m = paramdataPtr->parametervector[27]; /* 0.0881512 */
    KAck_AcP_m = paramdataPtr->parametervector[28]; /* 0.0742685 */
    KAck_eq = paramdataPtr->parametervector[29]; /* 222.087 */
    KAcs_ACE = paramdataPtr->parametervector[30]; /* 0.0233726 */
    KCS_AcCoA = paramdataPtr->parametervector[31]; /* 0.172622 */
    KCS_OAA = paramdataPtr->parametervector[32]; /* 0.0127217 */
    KCS_OAA_AcCoA = paramdataPtr->parametervector[33]; /* 0.0177911 */
    KCS_aKG = paramdataPtr->parametervector[34]; /* 0.271462 */
    KCraFBP = paramdataPtr->parametervector[35]; /* 0.0322704 */
    KCrpcAMP = paramdataPtr->parametervector[36]; /* 0.368586 */
    KCya_EIIAP = paramdataPtr->parametervector[37]; /* 0.00242536 */
    KEda_GAP_m = paramdataPtr->parametervector[38]; /* 1.37717 */
    KEda_KDPG_m = paramdataPtr->parametervector[39]; /* 0.350113 */
    KEda_PYR_m = paramdataPtr->parametervector[40]; /* 11.6552 */
    KEda_eq = paramdataPtr->parametervector[41]; /* 0.500041 */
    KEdd_6PG_m = paramdataPtr->parametervector[42]; /* 0.347444 */
    KEdd_KDPG_m = paramdataPtr->parametervector[43]; /* 0.648821 */
    KEdd_eq = paramdataPtr->parametervector[44]; /* 1000.01 */
    KFba_DHAP = paramdataPtr->parametervector[45]; /* 0.0879915 */
    KFba_FBP = paramdataPtr->parametervector[46]; /* 0.13302 */
    KFba_GAP = paramdataPtr->parametervector[47]; /* 0.0879908 */
    KFba_GAP_inh = paramdataPtr->parametervector[48]; /* 0.771367 */
    KFba_eq = paramdataPtr->parametervector[49]; /* 0.33082 */
    KFbp_FBP = paramdataPtr->parametervector[50]; /* 0.00255028 */
    KFbp_PEP = paramdataPtr->parametervector[51]; /* 0.225362 */
    KFum_FUM_m = paramdataPtr->parametervector[52]; /* 0.0554577 */
    KFum_eq = paramdataPtr->parametervector[53]; /* 18.9482 */
    KG6PDH_G6P = paramdataPtr->parametervector[54]; /* 0.0699847 */
    KG6PDH_NADP = paramdataPtr->parametervector[55]; /* 0.00607658 */
    KG6PDH_NADPH_g6pinh = paramdataPtr->parametervector[56]; /* 0.192067 */
    KG6PDH_NADPH_nadpinh = paramdataPtr->parametervector[57]; /* 0.0251647 */
    KGAPDH_GAP = paramdataPtr->parametervector[58]; /* 0.104745 */
    KGAPDH_NAD = paramdataPtr->parametervector[59]; /* 0.449961 */
    KGAPDH_NADH = paramdataPtr->parametervector[60]; /* 0.0200034 */
    KGAPDH_PGP = paramdataPtr->parametervector[61]; /* 0.0995525 */
    KGAPDH_eq = paramdataPtr->parametervector[62]; /* 0.217878 */
    KGlk_ATP_m = paramdataPtr->parametervector[63]; /* 0.799718 */
    KGlk_G6P_i = paramdataPtr->parametervector[64]; /* 14.9942 */
    KGlk_GLC_m = paramdataPtr->parametervector[65]; /* 0.219961 */
    KICDH_ICIT = paramdataPtr->parametervector[66]; /* 9.55731e-05 */
    KICDH_PEP = paramdataPtr->parametervector[67]; /* 0.209919 */
    KIcl_3PG = paramdataPtr->parametervector[68]; /* 0.492465 */
    KIcl_ICIT = paramdataPtr->parametervector[69]; /* 0.0257752 */
    KIcl_PEP = paramdataPtr->parametervector[70]; /* 0.025628 */
    KIcl_aKG = paramdataPtr->parametervector[71]; /* 0.202044 */
    KMDH_MAL_I = paramdataPtr->parametervector[72]; /* 1.32753 */
    KMDH_MAL_m = paramdataPtr->parametervector[73]; /* 1.33028 */
    KMDH_NADH_I = paramdataPtr->parametervector[74]; /* 0.0242784 */
    KMDH_NADH_m = paramdataPtr->parametervector[75]; /* 0.0168477 */
    KMDH_NAD_I = paramdataPtr->parametervector[76]; /* 0.309937 */
    KMDH_NAD_II = paramdataPtr->parametervector[77]; /* 0.330748 */
    KMDH_NAD_m = paramdataPtr->parametervector[78]; /* 0.099998 */
    KMDH_OAA_I = paramdataPtr->parametervector[79]; /* 0.25207 */
    KMDH_OAA_II = paramdataPtr->parametervector[80]; /* 0.396105 */
    KMDH_OAA_m = paramdataPtr->parametervector[81]; /* 0.269873 */
    KMDH_eq = paramdataPtr->parametervector[82]; /* 0.557468 */
    KMS_AcCoA = paramdataPtr->parametervector[83]; /* 0.381655 */
    KMS_GOX = paramdataPtr->parametervector[84]; /* 0.361087 */
    KMS_GOX_AcCoA = paramdataPtr->parametervector[85]; /* 0.311119 */
    KMez_AcCoA = paramdataPtr->parametervector[86]; /* 2.99104 */
    KMez_MAL = paramdataPtr->parametervector[87]; /* 0.00293864 */
    KMez_cAMP = paramdataPtr->parametervector[88]; /* 2.40366 */
    KNonPTS_I = paramdataPtr->parametervector[89]; /* 0.009166 */
    KNonPTS_S = paramdataPtr->parametervector[90]; /* 2.34579 */
    KPDH_AcCoA_m = paramdataPtr->parametervector[91]; /* 0.0080012 */
    KPDH_CoA_m = paramdataPtr->parametervector[92]; /* 0.00698447 */
    KPDH_NADH_m = paramdataPtr->parametervector[93]; /* 0.248943 */
    KPDH_NAD_m = paramdataPtr->parametervector[94]; /* 0.399784 */
    KPDH_PYR_m = paramdataPtr->parametervector[95]; /* 1.20561 */
    KPDH_i = paramdataPtr->parametervector[96]; /* 21.7028 */
    KPTS_EIIA = paramdataPtr->parametervector[97]; /* 0.00660253 */
    KPTS_GLC = paramdataPtr->parametervector[98]; /* 0.0127905 */
    KPck_ADP_i = paramdataPtr->parametervector[99]; /* 0.0388343 */
    KPck_ATP_I = paramdataPtr->parametervector[100]; /* 0.0396311 */
    KPck_ATP_i = paramdataPtr->parametervector[101]; /* 0.0415227 */
    KPck_OAA = paramdataPtr->parametervector[102]; /* 0.669902 */
    KPck_OAA_I = paramdataPtr->parametervector[103]; /* 0.532385 */
    KPck_PEP = paramdataPtr->parametervector[104]; /* 0.0701171 */
    KPck_PEP_i = paramdataPtr->parametervector[105]; /* 0.059769 */
    KPdhRPYR = paramdataPtr->parametervector[106]; /* 0.0856307 */
    KPfk_ADP_a = paramdataPtr->parametervector[107]; /* 238.935 */
    KPfk_ADP_b = paramdataPtr->parametervector[108]; /* 0.250106 */
    KPfk_ADP_c = paramdataPtr->parametervector[109]; /* 0.359964 */
    KPfk_AMP_a = paramdataPtr->parametervector[110]; /* 6.05077 */
    KPfk_AMP_b = paramdataPtr->parametervector[111]; /* 0.0207231 */
    KPfk_ATP_s = paramdataPtr->parametervector[112]; /* 0.160031 */
    KPfk_F6P_s = paramdataPtr->parametervector[113]; /* 0.0231651 */
    KPfk_PEP = paramdataPtr->parametervector[114]; /* 1.93645 */
    KPgi_F6P = paramdataPtr->parametervector[115]; /* 0.199983 */
    KPgi_F6P_6pginh = paramdataPtr->parametervector[116]; /* 0.19997 */
    KPgi_G6P = paramdataPtr->parametervector[117]; /* 2.45964 */
    KPgi_G6P_6pginh = paramdataPtr->parametervector[118]; /* 0.200129 */
    KPgi_eq = paramdataPtr->parametervector[119]; /* 0.717321 */
    KPgl_6PGL_m = paramdataPtr->parametervector[120]; /* 0.0229252 */
    KPgl_6PG_m = paramdataPtr->parametervector[121]; /* 9.96181 */
    KPgl_eq = paramdataPtr->parametervector[122]; /* 42.7993 */
    KPgl_h1 = paramdataPtr->parametervector[123]; /* 0.00225241 */
    KPgl_h2 = paramdataPtr->parametervector[124]; /* 9.71558e-06 */
    KPpc_FBP = paramdataPtr->parametervector[125]; /* 0.126971 */
    KPpc_PEP = paramdataPtr->parametervector[126]; /* 0.0385597 */
    KPps_PEP = paramdataPtr->parametervector[127]; /* 0.000604094 */
    KPps_PYR = paramdataPtr->parametervector[128]; /* 0.0006405 */
    KPta_AcCoA_i = paramdataPtr->parametervector[129]; /* 0.088354 */
    KPta_AcP_i = paramdataPtr->parametervector[130]; /* 0.0678373 */
    KPta_AcP_m = paramdataPtr->parametervector[131]; /* 0.523402 */
    KPta_CoA_i = paramdataPtr->parametervector[132]; /* 0.0244299 */
    KPta_Pi_i = paramdataPtr->parametervector[133]; /* 1.61647 */
    KPta_Pi_m = paramdataPtr->parametervector[134]; /* 0.699177 */
    KPta_eq = paramdataPtr->parametervector[135]; /* 0.0326788 */
    KPyk_ADP = paramdataPtr->parametervector[136]; /* 0.213774 */
    KPyk_AMP = paramdataPtr->parametervector[137]; /* 0.209251 */
    KPyk_ATP = paramdataPtr->parametervector[138]; /* 22.4085 */
    KPyk_FBP = paramdataPtr->parametervector[139]; /* 0.190009 */
    KPyk_PEP = paramdataPtr->parametervector[140]; /* 0.370704 */
    KR5PI_eq = paramdataPtr->parametervector[141]; /* 0.484213 */
    KRu5P_eq = paramdataPtr->parametervector[142]; /* 1.62946 */
    KSDH_SUC_m = paramdataPtr->parametervector[143]; /* 0.220018 */
    KSDH_eq = paramdataPtr->parametervector[144]; /* 8.67636 */
    KTal_eq = paramdataPtr->parametervector[145]; /* 1.17897 */
    KTktA_eq = paramdataPtr->parametervector[146]; /* 1.2001 */
    KTktB_eq = paramdataPtr->parametervector[147]; /* 10 */
    KaKGDH_CoA_m = paramdataPtr->parametervector[148]; /* 0.00391074 */
    KaKGDH_NADH_I = paramdataPtr->parametervector[149]; /* 0.0180044 */
    KaKGDH_NAD_m = paramdataPtr->parametervector[150]; /* 0.0700101 */
    KaKGDH_SUC_I = paramdataPtr->parametervector[151]; /* 0.533941 */
    KaKGDH_Z = paramdataPtr->parametervector[152]; /* 1.41777 */
    KaKGDH_aKG_I = paramdataPtr->parametervector[153]; /* 1.78176 */
    KaKGDH_aKG_m = paramdataPtr->parametervector[154]; /* 0.99931 */
    KaceBAK_Cra = paramdataPtr->parametervector[155]; /* 0.00018881 */
    KaceBAK_Crp = paramdataPtr->parametervector[156]; /* 4.56257 */
    KaceBAK_DNA = paramdataPtr->parametervector[157]; /* 1.13147e-06 */
    KaceBAK_GOX = paramdataPtr->parametervector[158]; /* 0.00566376 */
    KaceBAK_PYR = paramdataPtr->parametervector[159]; /* 0.249242 */
    KaceBAK_PYRprime = paramdataPtr->parametervector[160]; /* 0.00858996 */
    Kacs_Crp = paramdataPtr->parametervector[161]; /* 0.00846341 */
    KcAMPdegr_cAMP = paramdataPtr->parametervector[162]; /* 0.0662446 */
    KfbaA_Cra = paramdataPtr->parametervector[163]; /* 0.00373129 */
    KfbaA_Crp = paramdataPtr->parametervector[164]; /* 0.03915 */
    Kfbp_Cra = paramdataPtr->parametervector[165]; /* 0.000122915 */
    KfumABC_Crp = paramdataPtr->parametervector[166]; /* 0.122193 */
    KgapA_Cra = paramdataPtr->parametervector[167]; /* 0.00494581 */
    KgapA_Crp = paramdataPtr->parametervector[168]; /* 0.0296522 */
    Kglk_Cra = paramdataPtr->parametervector[169]; /* 2.91919e-08 */
    KgltA_Crp = paramdataPtr->parametervector[170]; /* 0.0312059 */
    KicdA_Cra = paramdataPtr->parametervector[171]; /* 5.3399e-05 */
    Kmdh_Crp = paramdataPtr->parametervector[172]; /* 0.255211 */
    KpckA_Cra = paramdataPtr->parametervector[173]; /* 0.000275018 */
    Kpdh_PdhR = paramdataPtr->parametervector[174]; /* 1.24734e-05 */
    KpfkA_Cra = paramdataPtr->parametervector[175]; /* 2.41644e-08 */
    KppsA_Cra = paramdataPtr->parametervector[176]; /* 0.00062288 */
    KpykF_Cra = paramdataPtr->parametervector[177]; /* 0.000117808 */
    KsdhCDAB_Crp = paramdataPtr->parametervector[178]; /* 0.045294 */
    KsucAB_Crp = paramdataPtr->parametervector[179]; /* 0.117486 */
    LAceK = paramdataPtr->parametervector[180]; /* 7.59785e+07 */
    LFbp = paramdataPtr->parametervector[181]; /* 1.11817e+07 */
    LICDH = paramdataPtr->parametervector[182]; /* 126.46 */
    LIcl = paramdataPtr->parametervector[183]; /* 33725.8 */
    LMez = paramdataPtr->parametervector[184]; /* 197779 */
    LPfk = paramdataPtr->parametervector[185]; /* 1.26952e+06 */
    LPpc = paramdataPtr->parametervector[186]; /* 3.8947e+06 */
    LPps = paramdataPtr->parametervector[187]; /* 4.03223e-80 */
    LPyk = paramdataPtr->parametervector[188]; /* 1000.15 */
    LaceBAK = paramdataPtr->parametervector[189]; /* 1243.52 */
    NAD = paramdataPtr->parametervector[190]; /* 2.6 */
    NADH = paramdataPtr->parametervector[191]; /* 0.083 */
    NADP = paramdataPtr->parametervector[192]; /* 0.0021 */
    NADPH = paramdataPtr->parametervector[193]; /* 0.12 */
    POratio = paramdataPtr->parametervector[194]; /* 3.30385 */
    POratio_prime = paramdataPtr->parametervector[195]; /* 1.50583 */
    PdhRtotal = paramdataPtr->parametervector[196]; /* 6.66e-05 */
    Pi = paramdataPtr->parametervector[197]; /* 10 */
    SS_Mez_ACE = paramdataPtr->parametervector[198]; /* 0.0334723 */
    SS_Mez_GLC = paramdataPtr->parametervector[199]; /* 0.0133738 */
    SS_Ppc_ACE = paramdataPtr->parametervector[200]; /* 0.000852337 */
    SS_Ppc_GLC = paramdataPtr->parametervector[201]; /* 0.00323703 */
    VFba_blf = paramdataPtr->parametervector[202]; /* 2.00107 */
    aceBAK_DNA = paramdataPtr->parametervector[203]; /* 5.15e-07 */
    kATP = paramdataPtr->parametervector[204]; /* 1.30699e-05 */
    kAceKki_cat = paramdataPtr->parametervector[205]; /* 9.19676e+15 */
    kAceKph_cat = paramdataPtr->parametervector[206]; /* 8.85104e+12 */
    kAcs_cat = paramdataPtr->parametervector[207]; /* 116474 */
    kBM_ACE_AcCoA = paramdataPtr->parametervector[208]; /* 164.317 */
    kBM_ACE_E4P = paramdataPtr->parametervector[209]; /* 284.114 */
    kBM_ACE_F6P = paramdataPtr->parametervector[210]; /* 321.291 */
    kBM_ACE_FUM = paramdataPtr->parametervector[211]; /* 143.731 */
    kBM_ACE_G6P = paramdataPtr->parametervector[212]; /* 118.486 */
    kBM_ACE_GAP = paramdataPtr->parametervector[213]; /* 734.492 */
    kBM_ACE_OAA = paramdataPtr->parametervector[214]; /* 7011.8 */
    kBM_ACE_PEP = paramdataPtr->parametervector[215]; /* 202.218 */
    kBM_ACE_PYR = paramdataPtr->parametervector[216]; /* 12876.3 */
    kBM_ACE_R5P = paramdataPtr->parametervector[217]; /* 240.315 */
    kBM_ACE_SUC = paramdataPtr->parametervector[218]; /* 150.808 */
    kBM_ACE_aKG = paramdataPtr->parametervector[219]; /* 324.664 */
    kBM_GLC_AcCoA = paramdataPtr->parametervector[220]; /* 2656.7 */
    kBM_GLC_E4P = paramdataPtr->parametervector[221]; /* 603.434 */
    kBM_GLC_F6P = paramdataPtr->parametervector[222]; /* 966.423 */
    kBM_GLC_FUM = paramdataPtr->parametervector[223]; /* 1091.82 */
    kBM_GLC_G6P = paramdataPtr->parametervector[224]; /* 52.0836 */
    kBM_GLC_GAP = paramdataPtr->parametervector[225]; /* 375.616 */
    kBM_GLC_OAA = paramdataPtr->parametervector[226]; /* 19606.2 */
    kBM_GLC_PEP = paramdataPtr->parametervector[227]; /* 708.44 */
    kBM_GLC_PYR = paramdataPtr->parametervector[228]; /* 707.651 */
    kBM_GLC_R5P = paramdataPtr->parametervector[229]; /* 247.339 */
    kBM_GLC_SUC = paramdataPtr->parametervector[230]; /* 2467.94 */
    kBM_GLC_aKG = paramdataPtr->parametervector[231]; /* 4673.85 */
    kCS_cat = paramdataPtr->parametervector[232]; /* 792390 */
    kFba_cat = paramdataPtr->parametervector[233]; /* 1.60675e+06 */
    kFbp_cat = paramdataPtr->parametervector[234]; /* 2.41475e+06 */
    kFum1_cat = paramdataPtr->parametervector[235]; /* 681983 */
    kFum2_cat = paramdataPtr->parametervector[236]; /* 681983 */
    kGAPDH_cat = paramdataPtr->parametervector[237]; /* 8.79789e+07 */
    kGlk_cat = paramdataPtr->parametervector[238]; /* 1.98397e+06 */
    kICDH_cat = paramdataPtr->parametervector[239]; /* 211215 */
    kIcl_cat = paramdataPtr->parametervector[240]; /* 1.09342e+06 */
    kMDH1_cat = paramdataPtr->parametervector[241]; /* 328277 */
    kMDH2_cat = paramdataPtr->parametervector[242]; /* 328277 */
    kMS_cat = paramdataPtr->parametervector[243]; /* 13791.5 */
    kMez_cat = paramdataPtr->parametervector[244]; /* 1.65343e+06 */
    kPDH_cat = paramdataPtr->parametervector[245]; /* 6.59379e+07 */
    kPTS1 = paramdataPtr->parametervector[246]; /* 42555.9 */
    kPck_cat = paramdataPtr->parametervector[247]; /* 2.60105e+06 */
    kPfk_cat = paramdataPtr->parametervector[248]; /* 1.49637e+10 */
    kPpc_cat = paramdataPtr->parametervector[249]; /* 1.0471e+07 */
    kPps_cat = paramdataPtr->parametervector[250]; /* 1070.78 */
    kPyk_cat = paramdataPtr->parametervector[251]; /* 17530.3 */
    kSDH1_cat = paramdataPtr->parametervector[252]; /* 22311.9 */
    kSDH2_cat = paramdataPtr->parametervector[253]; /* 22311.9 */
    kaKGDH_cat = paramdataPtr->parametervector[254]; /* 2.04051e+08 */
    kaceBAK_cat_IclR = paramdataPtr->parametervector[255]; /* 3.36933 */
    kdegr = paramdataPtr->parametervector[256]; /* 0.265324 */
    kexpr = paramdataPtr->parametervector[257]; /* 4.85939 */
    kmPTS1 = paramdataPtr->parametervector[258]; /* 12994.3 */
    nAceK = paramdataPtr->parametervector[259]; /* 2 */
    nCraFBP = paramdataPtr->parametervector[260]; /* 2 */
    nCrpcAMP = paramdataPtr->parametervector[261]; /* 1 */
    nFbp = paramdataPtr->parametervector[262]; /* 4 */
    nICDH = paramdataPtr->parametervector[263]; /* 2 */
    nIcl = paramdataPtr->parametervector[264]; /* 4 */
    nMez = paramdataPtr->parametervector[265]; /* 1.92774 */
    nPdhRPYR = paramdataPtr->parametervector[266]; /* 1 */
    nPfk = paramdataPtr->parametervector[267]; /* 4 */
    nPpc = paramdataPtr->parametervector[268]; /* 3 */
    nPps = paramdataPtr->parametervector[269]; /* 2 */
    nPyk = paramdataPtr->parametervector[270]; /* 4 */
    nacs = paramdataPtr->parametervector[271]; /* 2.31 */
    nfumABC = paramdataPtr->parametervector[272]; /* 0.74 */
    ngltA = paramdataPtr->parametervector[273]; /* 1.07 */
    nsdhCDAB = paramdataPtr->parametervector[274]; /* 0.74 */
    nsucAB = paramdataPtr->parametervector[275]; /* 0.74 */
    pH = paramdataPtr->parametervector[276]; /* 7.5 */
    pH_Eda_m = paramdataPtr->parametervector[277]; /* 7.49978 */
    pH_Edd_m = paramdataPtr->parametervector[278]; /* 7.15931 */
    pK_Eda = paramdataPtr->parametervector[279]; /* 13.1433 */
    pK_Edd = paramdataPtr->parametervector[280]; /* 3.42317 */
    rho = paramdataPtr->parametervector[281]; /* 564 */
    v6PGDH_max = paramdataPtr->parametervector[282]; /* 193585 */
    vAck_max = paramdataPtr->parametervector[283]; /* 533045 */
    vCya_max = paramdataPtr->parametervector[284]; /* 13.2427 */
    vEda_max = paramdataPtr->parametervector[285]; /* 582.202 */
    vEdd_max = paramdataPtr->parametervector[286]; /* 944.737 */
    vG6PDH_max = paramdataPtr->parametervector[287]; /* 53780.5 */
    vNonPTS_max = paramdataPtr->parametervector[288]; /* 4679.47 */
    vPTS4_max = paramdataPtr->parametervector[289]; /* 10119.9 */
    vPgi_max = paramdataPtr->parametervector[290]; /* 3.62077e+06 */
    vPgl_max = paramdataPtr->parametervector[291]; /* 45257.5 */
    vPta_max = paramdataPtr->parametervector[292]; /* 5915.4 */
    vR5PI_max = paramdataPtr->parametervector[293]; /* 42445.8 */
    vRu5P_max = paramdataPtr->parametervector[294]; /* 23336.1 */
    vTal_max = paramdataPtr->parametervector[295]; /* 42057.1 */
    vTktA_max = paramdataPtr->parametervector[296]; /* 7656.02 */
    vTktB_max = paramdataPtr->parametervector[297]; /* 285535 */
    vaceBAK_Cra_bound = paramdataPtr->parametervector[298]; /* 0.0596534 */
    vaceBAK_Cra_unbound = paramdataPtr->parametervector[299]; /* 7.53474e-05 */
    vaceBAK_Crp_bound = paramdataPtr->parametervector[300]; /* 1.59563e-05 */
    vaceBAK_Crp_unbound = paramdataPtr->parametervector[301]; /* 0.00262305 */
    vacs_Crp_bound = paramdataPtr->parametervector[302]; /* 0.000592461 */
    vacs_Crp_unbound = paramdataPtr->parametervector[303]; /* 0 */
    vcAMPdegr_max = paramdataPtr->parametervector[304]; /* 9.91179 */
    vfbaA_Cra_bound = paramdataPtr->parametervector[305]; /* 0 */
    vfbaA_Cra_unbound = paramdataPtr->parametervector[306]; /* 0.0173187 */
    vfbaA_Crp_bound = paramdataPtr->parametervector[307]; /* 0.0150909 */
    vfbaA_Crp_unbound = paramdataPtr->parametervector[308]; /* 0 */
    vfbp_Cra_bound = paramdataPtr->parametervector[309]; /* 0.000701946 */
    vfbp_Cra_unbound = paramdataPtr->parametervector[310]; /* 0 */
    vfumABC_Crp_bound = paramdataPtr->parametervector[311]; /* 0.0540414 */
    vfumABC_Crp_unbound = paramdataPtr->parametervector[312]; /* 0 */
    vgapA_Cra_bound = paramdataPtr->parametervector[313]; /* 0 */
    vgapA_Cra_unbound = paramdataPtr->parametervector[314]; /* 0.0156125 */
    vgapA_Crp_bound = paramdataPtr->parametervector[315]; /* 0.0242728 */
    vgapA_Crp_unbound = paramdataPtr->parametervector[316]; /* 0 */
    vglk_Cra_bound = paramdataPtr->parametervector[317]; /* 0.00139121 */
    vglk_Cra_unbound = paramdataPtr->parametervector[318]; /* 0.184677 */
    vgltA_Crp_bound = paramdataPtr->parametervector[319]; /* 0.0505062 */
    vgltA_Crp_unbound = paramdataPtr->parametervector[320]; /* 0 */
    vicdA_Cra_bound = paramdataPtr->parametervector[321]; /* 0.0159016 */
    vicdA_Cra_unbound = paramdataPtr->parametervector[322]; /* 0.00401129 */
    vmdh_Crp_bound = paramdataPtr->parametervector[323]; /* 0.158286 */
    vmdh_Crp_unbound = paramdataPtr->parametervector[324]; /* 0 */
    vpckA_Cra_bound = paramdataPtr->parametervector[325]; /* 0.0110092 */
    vpckA_Cra_unbound = paramdataPtr->parametervector[326]; /* 0 */
    vpdh_PdhR_bound = paramdataPtr->parametervector[327]; /* 9.45249e-06 */
    vpdh_PdhR_unbound = paramdataPtr->parametervector[328]; /* 0.00156527 */
    vpfkA_Cra_bound = paramdataPtr->parametervector[329]; /* 0.00272437 */
    vpfkA_Cra_unbound = paramdataPtr->parametervector[330]; /* 0.0754715 */
    vppsA_Cra_bound = paramdataPtr->parametervector[331]; /* 0.113449 */
    vppsA_Cra_unbound = paramdataPtr->parametervector[332]; /* 0 */
    vpykF_Cra_bound = paramdataPtr->parametervector[333]; /* 8.98315e-05 */
    vpykF_Cra_unbound = paramdataPtr->parametervector[334]; /* 0.0038128 */
    vsdhCDAB_Crp_bound = paramdataPtr->parametervector[335]; /* 0.357936 */
    vsdhCDAB_Crp_unbound = paramdataPtr->parametervector[336]; /* 0 */
    vsucAB_Crp_bound = paramdataPtr->parametervector[337]; /* 0.0282708 */
    vsucAB_Crp_unbound = paramdataPtr->parametervector[338]; /* 0 */
    default0 = paramdataPtr->parametervector[339]; /* 1 */
    CraFBP = Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP));
    B = 1.0+ADP/KPfk_ADP_a+AMP/KPfk_AMP_a;
    A = 1.0+PEP/KPfk_PEP+ADP/KPfk_ADP_b+AMP/KPfk_AMP_b;
    Cra = Cratotal-CraFBP;
    CrpcAMP = Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP));
    Crp = Crptotal-CrpcAMP;
    H = pow(10.0,-pH)*1000.0;
    EIIA = EIIAtotal-EIIAP;
    alphaGLC = GLCex/(GLCex+KPTS_GLC);
    alphaACE = ACEex/(ACEex+KAcs_ACE)*(1.0-alphaGLC);
    SS_Ppc = alphaGLC*SS_Ppc_GLC+alphaACE*SS_Ppc_ACE;
    SS_Mez = alphaGLC*SS_Mez_GLC+alphaACE*SS_Mez_ACE;
    OP_NADH = (vE_GAPDH+vE_PDH+vE_aKGDH+vE_MDH)*POratio;
    OP_FADH2 = vE_SDH*POratio_prime;
    PdhRPYR = PdhRtotal*pow(PYR,nPdhRPYR)/(pow(PYR,nPdhRPYR)+pow(KPdhRPYR,nPdhRPYR));
    PdhR = PdhRtotal-PdhRPYR;
    QEda_pH = 1.0+2.0*pow(10.0,pH_Eda_m-pK_Eda)/(1.0+pow(10.0,pH-pK_Eda)+pow(10.0,2.0*pH_Eda_m-pH-pK_Eda));
    QEdd_pH = 1.0+2.0*pow(10.0,pH_Edd_m-pK_Edd)/(1.0+pow(10.0,pH-pK_Edd)+pow(10.0,2.0*pH_Edd_m-pH-pK_Edd));
    vMDH1_max = MDH*kMDH1_cat;
    vFum2_max = Fum*kFum2_cat;
    vSDH2_max = SDH*kSDH2_cat;
    vSDH1_max = SDH*kSDH1_cat;
    vATP = OP_NADH+OP_FADH2-vE_Glk-vE_Pfk+vE_GAPDH+vE_Pyk-vE_Pps+vE_Ack-vE_Acs+vE_aKGDH-vE_Pck-vE_AceKki-vE_Cya;
    mu = kATP*vATP;
    vMDH2_max = MDH*kMDH2_cat;
    vFum1_max = Fum*kFum1_cat;
    vBM_AcCoA = default0*Function_for_vBM_AcCoA(AcCoA,alphaACE,alphaGLC,default0,kBM_ACE_AcCoA,kBM_GLC_AcCoA);
    vBM_E4P = default0*Function_for_vBM_E4P(E4P,alphaACE,alphaGLC,default0,kBM_ACE_E4P,kBM_GLC_E4P);
    vBM_F6P = default0*Function_for_vBM_F6P(F6P,alphaACE,alphaGLC,default0,kBM_ACE_F6P,kBM_GLC_F6P);
    vBM_FUM = default0*Function_for_vBM_FUM(FUM,alphaACE,alphaGLC,default0,kBM_ACE_FUM,kBM_GLC_FUM);
    vBM_G6P = default0*Function_for_vBM_G6P(G6P,alphaACE,alphaGLC,default0,kBM_ACE_G6P,kBM_GLC_G6P);
    vBM_GAP = default0*Function_for_vBM_GAP(GAP,alphaACE,alphaGLC,default0,kBM_ACE_GAP,kBM_GLC_GAP);
    vBM_OAA = default0*Function_for_vBM_OAA(OAA,alphaACE,alphaGLC,default0,kBM_ACE_OAA,kBM_GLC_OAA);
    vBM_PEP = default0*Function_for_vBM_PEP(PEP,alphaACE,alphaGLC,default0,kBM_ACE_PEP,kBM_GLC_PEP);
    vBM_PYR = default0*Function_for_vBM_PYR(PYR,alphaACE,alphaGLC,default0,kBM_ACE_PYR,kBM_GLC_PYR);
    vBM_R5P = default0*Function_for_vBM_R5P(R5P,alphaACE,alphaGLC,default0,kBM_ACE_R5P,kBM_GLC_R5P);
    vBM_SUC = default0*Function_for_vBM_SUC(SUC,alphaACE,alphaGLC,default0,kBM_ACE_SUC,kBM_GLC_SUC);
    vBM_aKG = default0*Function_for_vBM_aKG(aKG,alphaACE,alphaGLC,default0,kBM_ACE_aKG,kBM_GLC_aKG);
    vD_6PG = default0*Function_for_vD_6PG(default0,mu,sixPG);
    vD_6PGL = default0*Function_for_vD_6PGL(default0,mu,sixPGL);
    vD_ACEex = default0*Function_for_vD_ACEex(ACEex,D,default0);
    vD_AcCoA = default0*Function_for_vD_AcCoA(AcCoA,default0,mu);
    vD_AcP = default0*Function_for_vD_AcP(AcP,default0,mu);
    vD_AceK = default0*Function_for_vD_AceK(AceK,default0,kdegr,mu);
    vD_Acs = default0*Function_for_vD_Acs(Acs,default0,kdegr,mu);
    vD_CS = default0*Function_for_vD_CS(CS,default0,kdegr,mu);
    vD_E4P = default0*Function_for_vD_E4P(E4P,default0,mu);
    vD_F6P = default0*Function_for_vD_F6P(F6P,default0,mu);
    vD_FBP = default0*Function_for_vD_FBP(FBP,default0,mu);
    vD_FUM = default0*Function_for_vD_FUM(FUM,default0,mu);
    vD_Fba = default0*Function_for_vD_Fba(Fba,default0,kdegr,mu);
    vD_Fbp = default0*Function_for_vD_Fbp(Fbp,default0,kdegr,mu);
    vD_Fum = default0*Function_for_vD_Fum(Fum,default0,kdegr,mu);
    vD_G6P = default0*Function_for_vD_G6P(G6P,default0,mu);
    vD_GAP = default0*Function_for_vD_GAP(GAP,default0,mu);
    vD_GAPDH = default0*Function_for_vD_GAPDH(GAPDH,default0,kdegr,mu);
    vD_GLC = default0*Function_for_vD_GLC(GLC,default0,mu);
    vD_GLCex = default0*Function_for_vD_GLCex(D,GLCex,default0);
    vD_GLCfeed = default0*Function_for_vD_GLCfeed(D,GLCfeed,default0);
    vD_GOX = default0*Function_for_vD_GOX(GOX,default0,mu);
    vD_Glk = default0*Function_for_vD_Glk(Glk,default0,kdegr,mu);
    vD_ICDH = default0*Function_for_vD_ICDH(ICDH,default0,kdegr,mu);
    vD_ICDHP = default0*Function_for_vD_ICDHP(ICDHP,default0,kdegr,mu);
    vD_ICIT = default0*Function_for_vD_ICIT(ICIT,default0,mu);
    vD_Icl = default0*Function_for_vD_Icl(Icl,default0,kdegr,mu);
    vD_KDPG = default0*Function_for_vD_KDPG(KDPG,default0,mu);
    vD_MAL = default0*Function_for_vD_MAL(MAL,default0,mu);
    vD_MDH = default0*Function_for_vD_MDH(MDH,default0,kdegr,mu);
    vD_MS = default0*Function_for_vD_MS(MS,default0,kdegr,mu);
    vD_Mez = default0*Function_for_vD_Mez(Mez,default0,kdegr,mu);
    vD_OAA = default0*Function_for_vD_OAA(OAA,default0,mu);
    vD_PDH = default0*Function_for_vD_PDH(PDH,default0,kdegr,mu);
    vD_PEP = default0*Function_for_vD_PEP(PEP,default0,mu);
    vD_PYR = default0*Function_for_vD_PYR(PYR,default0,mu);
    vD_Pck = default0*Function_for_vD_Pck(Pck,default0,kdegr,mu);
    vD_Pfk = default0*Function_for_vD_Pfk(Pfk,default0,kdegr,mu);
    vD_Ppc = default0*Function_for_vD_Ppc(Ppc,default0,kdegr,mu);
    vD_Pps = default0*Function_for_vD_Pps(Pps,default0,kdegr,mu);
    vD_Pyk = default0*Function_for_vD_Pyk(Pyk,default0,kdegr,mu);
    vD_R5P = default0*Function_for_vD_R5P(R5P,default0,mu);
    vD_RU5P = default0*Function_for_vD_RU5P(RU5P,default0,mu);
    vD_S7P = default0*Function_for_vD_S7P(S7P,default0,mu);
    vD_SDH = default0*Function_for_vD_SDH(SDH,default0,kdegr,mu);
    vD_SUC = default0*Function_for_vD_SUC(SUC,default0,mu);
    vD_X = default0*Function_for_vD_X(D,X,default0);
    vD_X5P = default0*Function_for_vD_X5P(X5P,default0,mu);
    vD_aKG = default0*Function_for_vD_aKG(aKG,default0,mu);
    vD_aKGDH = default0*Function_for_vD_aKGDH(aKGDH,default0,kdegr,mu);
    vD_cAMP = default0*Function_for_vD_cAMP(cAMP,default0,mu);
    vE_6PGDH = default0*Function_for_vE_6PGDH(ATP,K6PGDH_6PG,K6PGDH_ATP_inh,K6PGDH_NADP,K6PGDH_NADPH_inh,NADP,NADPH,default0,sixPG,v6PGDH_max);
    vE_AceKki = default0*Function_for_vE_AceKki(AceK,GAP,GOX,ICDH,ICIT,KAceK_3PG,KAceK_GOX,KAceK_ICDH,KAceK_ICIT,KAceK_OAA,KAceK_PEP,KAceK_PYR,KAceK_aKG,LAceK,OAA,PEP,PYR,aKG,default0,kAceKki_cat,nAceK);
    vE_AceKph = default0*Function_for_vE_AceKph(AceK,GAP,ICDHP,KAceK_3PG,KAceK_ICDHP,KAceK_OAA,KAceK_PEP,KAceK_PYR,KAceK_aKG,LAceK,OAA,PEP,PYR,aKG,default0,kAceKph_cat,nAceK);
    vE_Ack = default0*Function_for_vE_Ack(ACEex,ADP,ATP,AcP,KAck_ACE_m,KAck_ADP_m,KAck_ATP_m,KAck_AcP_m,KAck_eq,default0,vAck_max);
    vE_Ack_medium = default0*Function_for_vE_Ack_medium(ACEex,ADP,ATP,AcP,KAck_ACE_m,KAck_ADP_m,KAck_ATP_m,KAck_AcP_m,KAck_eq,X,default0,rho,vAck_max);
    vE_Acs = default0*Function_for_vE_Acs(ACEex,Acs,KAcs_ACE,default0,kAcs_cat);
    vE_Acs_medium = default0*Function_for_vE_Acs_medium(ACEex,Acs,KAcs_ACE,X,default0,kAcs_cat,rho);
    vE_CS = default0*Function_for_vE_CS(AcCoA,CS,KCS_AcCoA,KCS_OAA,KCS_OAA_AcCoA,KCS_aKG,OAA,aKG,default0,kCS_cat);
    vE_Cya = default0*Function_for_vE_Cya(EIIAP,KCya_EIIAP,default0,vCya_max);
    vE_Eda = default0*Function_for_vE_Eda(GAP,KDPG,KEda_GAP_m,KEda_KDPG_m,KEda_PYR_m,KEda_eq,PYR,QEda_pH,default0,vEda_max);
    vE_Edd = default0*Function_for_vE_Edd(KDPG,KEdd_6PG_m,KEdd_KDPG_m,KEdd_eq,QEdd_pH,default0,sixPG,vEdd_max);
    vE_Fba = default0*Function_for_vE_Fba(FBP,Fba,GAP,KFba_DHAP,KFba_FBP,KFba_GAP,KFba_GAP_inh,KFba_eq,VFba_blf,default0,kFba_cat);
    vE_Fbp = default0*Function_for_vE_Fbp(FBP,Fbp,KFbp_FBP,KFbp_PEP,LFbp,PEP,default0,kFbp_cat,nFbp);
    vE_Fum = default0*Function_for_vE_Fum(FUM,KFum_FUM_m,KFum_eq,MAL,default0,vFum1_max,vFum2_max);
    vE_G6PDH = default0*Function_for_vE_G6PDH(G6P,KG6PDH_G6P,KG6PDH_NADP,KG6PDH_NADPH_g6pinh,KG6PDH_NADPH_nadpinh,NADP,NADPH,default0,vG6PDH_max);
    vE_GAPDH = default0*Function_for_vE_GAPDH(GAP,GAPDH,KGAPDH_GAP,KGAPDH_NAD,KGAPDH_NADH,KGAPDH_PGP,KGAPDH_eq,NAD,NADH,PEP,default0,kGAPDH_cat);
    vE_Glk = default0*Function_for_vE_Glk(ATP,G6P,GLC,Glk,KGlk_ATP_m,KGlk_G6P_i,KGlk_GLC_m,default0,kGlk_cat);
    vE_ICDH = default0*Function_for_vE_ICDH(ICDH,ICIT,KICDH_ICIT,KICDH_PEP,LICDH,PEP,default0,kICDH_cat,nICDH);
    vE_Icl = default0*Function_for_vE_Icl(GAP,ICIT,Icl,KIcl_3PG,KIcl_ICIT,KIcl_PEP,KIcl_aKG,LIcl,PEP,aKG,default0,kIcl_cat,nIcl);
    vE_MDH = default0*Function_for_vE_MDH(KMDH_MAL_I,KMDH_MAL_m,KMDH_NADH_I,KMDH_NADH_m,KMDH_NAD_I,KMDH_NAD_II,KMDH_NAD_m,KMDH_OAA_I,KMDH_OAA_II,KMDH_OAA_m,KMDH_eq,MAL,NAD,NADH,OAA,default0,vMDH1_max,vMDH2_max);
    vE_MS = default0*Function_for_vE_MS(AcCoA,GOX,KMS_AcCoA,KMS_GOX,KMS_GOX_AcCoA,MS,default0,kMS_cat);
    vE_Mez = default0*Function_for_vE_Mez(AcCoA,KMez_AcCoA,KMez_MAL,KMez_cAMP,LMez,MAL,Mez,cAMP,default0,kMez_cat,nMez);
    vE_PDH = default0*Function_for_vE_PDH(AcCoA,CoA,KPDH_AcCoA_m,KPDH_CoA_m,KPDH_NADH_m,KPDH_NAD_m,KPDH_PYR_m,KPDH_i,NAD,NADH,PDH,PYR,default0,kPDH_cat);
    vE_Pck = default0*Function_for_vE_Pck(ADP,ATP,KPck_ADP_i,KPck_ATP_I,KPck_ATP_i,KPck_OAA,KPck_OAA_I,KPck_PEP,KPck_PEP_i,OAA,PEP,Pck,default0,kPck_cat);
    vE_Pfk = default0*Function_for_vE_Pfk(A,ADP,ATP,B,F6P,KPfk_ADP_c,KPfk_ATP_s,KPfk_F6P_s,LPfk,Pfk,default0,kPfk_cat,nPfk);
    vE_Pgi = default0*Function_for_vE_Pgi(F6P,G6P,KPgi_F6P,KPgi_F6P_6pginh,KPgi_G6P,KPgi_G6P_6pginh,KPgi_eq,default0,sixPG,vPgi_max);
    vE_Pgl = default0*Function_for_vE_Pgl(H,KPgl_6PGL_m,KPgl_6PG_m,KPgl_eq,KPgl_h1,KPgl_h2,default0,sixPG,sixPGL,vPgl_max);
    vE_Ppc = default0*Function_for_vE_Ppc(FBP,KPpc_FBP,KPpc_PEP,LPpc,PEP,Ppc,default0,kPpc_cat,nPpc);
    vE_Pps = default0*Function_for_vE_Pps(KPps_PEP,KPps_PYR,LPps,PEP,PYR,Pps,default0,kPps_cat,nPps);
    vE_Pta = default0*Function_for_vE_Pta(AcCoA,AcP,CoA,KPta_AcCoA_i,KPta_AcP_i,KPta_AcP_m,KPta_CoA_i,KPta_Pi_i,KPta_Pi_m,KPta_eq,Pi,default0,vPta_max);
    vE_Pyk = default0*Function_for_vE_Pyk(ADP,AMP,ATP,FBP,KPyk_ADP,KPyk_AMP,KPyk_ATP,KPyk_FBP,KPyk_PEP,LPyk,PEP,Pyk,default0,kPyk_cat,nPyk);
    vE_R5PI = default0*Function_for_vE_R5PI(KR5PI_eq,R5P,RU5P,default0,vR5PI_max);
    vE_Ru5P = default0*Function_for_vE_Ru5P(KRu5P_eq,RU5P,X5P,default0,vRu5P_max);
    vE_SDH = default0*Function_for_vE_SDH(FUM,KSDH_SUC_m,KSDH_eq,SUC,default0,vSDH1_max,vSDH2_max);
    vE_Tal = default0*Function_for_vE_Tal(E4P,F6P,GAP,KTal_eq,S7P,default0,vTal_max);
    vE_TktA = default0*Function_for_vE_TktA(GAP,KTktA_eq,R5P,S7P,X5P,default0,vTktA_max);
    vE_TktB = default0*Function_for_vE_TktB(E4P,F6P,GAP,KTktB_eq,X5P,default0,vTktB_max);
    vE_aKGDH = default0*Function_for_vE_aKGDH(CoA,KaKGDH_CoA_m,KaKGDH_NADH_I,KaKGDH_NAD_m,KaKGDH_SUC_I,KaKGDH_Z,KaKGDH_aKG_I,KaKGDH_aKG_m,NAD,NADH,SUC,aKG,aKGDH,default0,kaKGDH_cat);
    vE_cAMPdegr = default0*Function_for_vE_cAMPdegr(KcAMPdegr_cAMP,cAMP,default0,vcAMPdegr_max);
    vG_aceA = default0*Function_for_vG_aceA(Cra,CrpcAMP,GOX,IclRtotal,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR,KaceBAK_PYRprime,LaceBAK,PYR,aceBAK_DNA,default0,kaceBAK_cat_IclR,kexpr,mu,vaceBAK_Cra_bound,vaceBAK_Cra_unbound,vaceBAK_Crp_bound,vaceBAK_Crp_unbound);
    vG_aceB = default0*Function_for_vG_aceB(Cra,CrpcAMP,Factor_aceB,GOX,IclRtotal,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR,KaceBAK_PYRprime,LaceBAK,PYR,aceBAK_DNA,default0,kaceBAK_cat_IclR,kexpr,mu,vaceBAK_Cra_bound,vaceBAK_Cra_unbound,vaceBAK_Crp_bound,vaceBAK_Crp_unbound);
    vG_aceK = default0*Function_for_vG_aceK(Cra,CrpcAMP,Factor_aceK,GOX,IclRtotal,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR,KaceBAK_PYRprime,LaceBAK,PYR,aceBAK_DNA,default0,kaceBAK_cat_IclR,kexpr,mu,vaceBAK_Cra_bound,vaceBAK_Cra_unbound,vaceBAK_Crp_bound,vaceBAK_Crp_unbound);
    vG_acs = default0*Function_for_vG_acs(CrpcAMP,Kacs_Crp,default0,kexpr,mu,nacs,vacs_Crp_bound,vacs_Crp_unbound);
    vG_fbaA = default0*Function_for_vG_fbaA(Cra,CrpcAMP,KfbaA_Cra,KfbaA_Crp,default0,kexpr,mu,vfbaA_Cra_bound,vfbaA_Cra_unbound,vfbaA_Crp_bound,vfbaA_Crp_unbound);
    vG_fbp = default0*Function_for_vG_fbp(Cra,Kfbp_Cra,default0,kexpr,mu,vfbp_Cra_bound,vfbp_Cra_unbound);
    vG_fumABC = default0*Function_for_vG_fumABC(CrpcAMP,KfumABC_Crp,default0,kexpr,mu,nfumABC,vfumABC_Crp_bound,vfumABC_Crp_unbound);
    vG_gapA = default0*Function_for_vG_gapA(Cra,CrpcAMP,KgapA_Cra,KgapA_Crp,default0,kexpr,mu,vgapA_Cra_bound,vgapA_Cra_unbound,vgapA_Crp_bound,vgapA_Crp_unbound);
    vG_glk = default0*Function_for_vG_glk(Cra,Kglk_Cra,default0,kexpr,mu,vglk_Cra_bound,vglk_Cra_unbound);
    vG_gltA = default0*Function_for_vG_gltA(CrpcAMP,KgltA_Crp,default0,kexpr,mu,ngltA,vgltA_Crp_bound,vgltA_Crp_unbound);
    vG_icdA = default0*Function_for_vG_icdA(Cra,KicdA_Cra,default0,kexpr,mu,vicdA_Cra_bound,vicdA_Cra_unbound);
    vG_maeB = default0*Function_for_vG_maeB(SS_Mez,default0,kdegr,mu);
    vG_mdh = default0*Function_for_vG_mdh(CrpcAMP,Kmdh_Crp,default0,kexpr,mu,vmdh_Crp_bound,vmdh_Crp_unbound);
    vG_pckA = default0*Function_for_vG_pckA(Cra,KpckA_Cra,default0,kexpr,mu,vpckA_Cra_bound,vpckA_Cra_unbound);
    vG_pdh = default0*Function_for_vG_pdh(Kpdh_PdhR,PdhR,default0,kexpr,mu,vpdh_PdhR_bound,vpdh_PdhR_unbound);
    vG_pfkA = default0*Function_for_vG_pfkA(Cra,KpfkA_Cra,default0,kexpr,mu,vpfkA_Cra_bound,vpfkA_Cra_unbound);
    vG_ppc = default0*Function_for_vG_ppc(SS_Ppc,default0,kdegr,mu);
    vG_ppsA = default0*Function_for_vG_ppsA(Cra,KppsA_Cra,default0,kexpr,mu,vppsA_Cra_bound,vppsA_Cra_unbound);
    vG_pykF = default0*Function_for_vG_pykF(Cra,KpykF_Cra,default0,kexpr,mu,vpykF_Cra_bound,vpykF_Cra_unbound);
    vG_sdhCDAB = default0*Function_for_vG_sdhCDAB(CrpcAMP,KsdhCDAB_Crp,default0,kexpr,mu,nsdhCDAB,vsdhCDAB_Crp_bound,vsdhCDAB_Crp_unbound);
    vG_sucAB = default0*Function_for_vG_sucAB(CrpcAMP,KsucAB_Crp,default0,kexpr,mu,nsucAB,vsucAB_Crp_bound,vsucAB_Crp_unbound);
    vNonPTS = default0*Function_for_vNonPTS(EIIA,GLCex,KNonPTS_I,KNonPTS_S,default0,vNonPTS_max);
    vNonPTS_medium = default0*Function_for_vNonPTS_medium(EIIA,GLCex,KNonPTS_I,KNonPTS_S,X,default0,rho,vNonPTS_max);
    vPTS1 = default0*Function_for_vPTS1(EIIA,EIIAP,PEP,PYR,default0,kPTS1,kmPTS1);
    vPTS4 = default0*Function_for_vPTS4(EIIAP,GLCex,KPTS_EIIA,KPTS_GLC,default0,vPTS4_max);
    vPTS4_medium = default0*Function_for_vPTS4_medium(EIIAP,GLCex,KPTS_EIIA,KPTS_GLC,X,default0,rho,vPTS4_max);
    vgrowth = default0*Function_for_vgrowth(X,default0,mu);
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = (-vD_ACEex+vE_Ack_medium-vE_Acs_medium)/default0;
    	DDTvector[1] = (-vBM_AcCoA-vD_AcCoA+vE_Acs-vE_CS-vE_MS+vE_PDH-vE_Pta)/default0;
    	DDTvector[2] = (-vD_AcP-vE_Ack+vE_Pta)/default0;
    	DDTvector[3] = (-vD_AceK+vG_aceK)/default0;
    	DDTvector[4] = (-vD_Acs+vG_acs)/default0;
    	DDTvector[5] = (-vD_CS+vG_gltA)/default0;
    	DDTvector[6] = (-vBM_E4P-vD_E4P+vE_Tal-vE_TktB)/default0;
    	DDTvector[7] = (+vPTS1-vPTS4)/default0;
    	DDTvector[8] = (-vBM_F6P-vD_F6P+vE_Fbp-vE_Pfk+vE_Pgi+vE_Tal+vE_TktB)/default0;
    	DDTvector[9] = (-vD_FBP-vE_Fba-vE_Fbp+vE_Pfk)/default0;
    	DDTvector[10] = (-vBM_FUM-vD_FUM-vE_Fum+vE_SDH)/default0;
    	DDTvector[11] = (-vD_Fba+vG_fbaA)/default0;
    	DDTvector[12] = (-vD_Fbp+vG_fbp)/default0;
    	DDTvector[13] = (-vD_Fum+vG_fumABC)/default0;
    	DDTvector[14] = (-vBM_G6P-vD_G6P-vE_G6PDH+vE_Glk-vE_Pgi+vPTS4)/default0;
    	DDTvector[15] = (-vBM_GAP-vD_GAP+vE_Eda+2.0*vE_Fba-vE_GAPDH-vE_Tal+vE_TktA+vE_TktB)/default0;
    	DDTvector[16] = (-vD_GAPDH+vG_gapA)/default0;
    	DDTvector[17] = (-vD_GLC-vE_Glk+vNonPTS)/default0;
    	DDTvector[18] = (-vD_GLCex-vNonPTS_medium-vPTS4_medium)/default0;
    	DDTvector[19] = (-vD_GLCfeed)/default0;
    	DDTvector[20] = (-vD_GOX+vE_Icl-vE_MS)/default0;
    	DDTvector[21] = (-vD_Glk+vG_glk)/default0;
    	DDTvector[22] = (-vD_ICDH-vE_AceKki+vE_AceKph+vG_icdA)/default0;
    	DDTvector[23] = (-vD_ICDHP+vE_AceKki-vE_AceKph)/default0;
    	DDTvector[24] = (-vD_ICIT+vE_CS-vE_ICDH-vE_Icl)/default0;
    	DDTvector[25] = (-vD_Icl+vG_aceA)/default0;
    	DDTvector[26] = (-vD_KDPG-vE_Eda+vE_Edd)/default0;
    	DDTvector[27] = (-vD_MAL+vE_Fum-vE_MDH+vE_MS-vE_Mez)/default0;
    	DDTvector[28] = (-vD_MDH+vG_mdh)/default0;
    	DDTvector[29] = (-vD_MS+vG_aceB)/default0;
    	DDTvector[30] = (-vD_Mez+vG_maeB)/default0;
    	DDTvector[31] = (-vBM_OAA-vD_OAA-vE_CS+vE_MDH-vE_Pck+vE_Ppc)/default0;
    	DDTvector[32] = (-vD_PDH+vG_pdh)/default0;
    	DDTvector[33] = (-vBM_PEP-vD_PEP+vE_GAPDH+vE_Pck-vE_Ppc+vE_Pps-vE_Pyk-vPTS1)/default0;
    	DDTvector[34] = (-vBM_PYR-vD_PYR+vE_Eda+vE_Mez-vE_PDH-vE_Pps+vE_Pyk+vPTS1)/default0;
    	DDTvector[35] = (-vD_Pck+vG_pckA)/default0;
    	DDTvector[36] = (-vD_Pfk+vG_pfkA)/default0;
    	DDTvector[37] = (-vD_Ppc+vG_ppc)/default0;
    	DDTvector[38] = (-vD_Pps+vG_ppsA)/default0;
    	DDTvector[39] = (-vD_Pyk+vG_pykF)/default0;
    	DDTvector[40] = (-vBM_R5P-vD_R5P+vE_R5PI-vE_TktA)/default0;
    	DDTvector[41] = (-vD_RU5P+vE_6PGDH-vE_R5PI-vE_Ru5P)/default0;
    	DDTvector[42] = (-vD_S7P-vE_Tal+vE_TktA)/default0;
    	DDTvector[43] = (-vD_SDH+vG_sdhCDAB)/default0;
    	DDTvector[44] = (-vBM_SUC-vD_SUC+vE_Icl-vE_SDH+vE_aKGDH)/default0;
    	DDTvector[45] = (-vD_X+vgrowth)/default0;
    	DDTvector[46] = (-vD_X5P+vE_Ru5P-vE_TktA-vE_TktB)/default0;
    	DDTvector[47] = (-vBM_aKG-vD_aKG+vE_ICDH-vE_aKGDH)/default0;
    	DDTvector[48] = (-vD_aKGDH+vG_sucAB)/default0;
    	DDTvector[49] = (-vD_cAMP+vE_Cya-vE_cAMPdegr)/default0;
    	DDTvector[50] = (-vD_6PG-vE_6PGDH-vE_Edd+vE_Pgl)/default0;
    	DDTvector[51] = (-vD_6PGL+vE_G6PDH-vE_Pgl)/default0;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = CraFBP;
        variableVector[1] = B;
        variableVector[2] = A;
        variableVector[3] = Cra;
        variableVector[4] = CrpcAMP;
        variableVector[5] = Crp;
        variableVector[6] = H;
        variableVector[7] = EIIA;
        variableVector[8] = alphaGLC;
        variableVector[9] = alphaACE;
        variableVector[10] = SS_Ppc;
        variableVector[11] = SS_Mez;
        variableVector[12] = OP_NADH;
        variableVector[13] = OP_FADH2;
        variableVector[14] = PdhRPYR;
        variableVector[15] = PdhR;
        variableVector[16] = QEda_pH;
        variableVector[17] = QEdd_pH;
        variableVector[18] = vMDH1_max;
        variableVector[19] = vFum2_max;
        variableVector[20] = vSDH2_max;
        variableVector[21] = vSDH1_max;
        variableVector[22] = vATP;
        variableVector[23] = mu;
        variableVector[24] = vMDH2_max;
        variableVector[25] = vFum1_max;
        reactionVector[0] = vBM_AcCoA;
        reactionVector[1] = vBM_E4P;
        reactionVector[2] = vBM_F6P;
        reactionVector[3] = vBM_FUM;
        reactionVector[4] = vBM_G6P;
        reactionVector[5] = vBM_GAP;
        reactionVector[6] = vBM_OAA;
        reactionVector[7] = vBM_PEP;
        reactionVector[8] = vBM_PYR;
        reactionVector[9] = vBM_R5P;
        reactionVector[10] = vBM_SUC;
        reactionVector[11] = vBM_aKG;
        reactionVector[12] = vD_6PG;
        reactionVector[13] = vD_6PGL;
        reactionVector[14] = vD_ACEex;
        reactionVector[15] = vD_AcCoA;
        reactionVector[16] = vD_AcP;
        reactionVector[17] = vD_AceK;
        reactionVector[18] = vD_Acs;
        reactionVector[19] = vD_CS;
        reactionVector[20] = vD_E4P;
        reactionVector[21] = vD_F6P;
        reactionVector[22] = vD_FBP;
        reactionVector[23] = vD_FUM;
        reactionVector[24] = vD_Fba;
        reactionVector[25] = vD_Fbp;
        reactionVector[26] = vD_Fum;
        reactionVector[27] = vD_G6P;
        reactionVector[28] = vD_GAP;
        reactionVector[29] = vD_GAPDH;
        reactionVector[30] = vD_GLC;
        reactionVector[31] = vD_GLCex;
        reactionVector[32] = vD_GLCfeed;
        reactionVector[33] = vD_GOX;
        reactionVector[34] = vD_Glk;
        reactionVector[35] = vD_ICDH;
        reactionVector[36] = vD_ICDHP;
        reactionVector[37] = vD_ICIT;
        reactionVector[38] = vD_Icl;
        reactionVector[39] = vD_KDPG;
        reactionVector[40] = vD_MAL;
        reactionVector[41] = vD_MDH;
        reactionVector[42] = vD_MS;
        reactionVector[43] = vD_Mez;
        reactionVector[44] = vD_OAA;
        reactionVector[45] = vD_PDH;
        reactionVector[46] = vD_PEP;
        reactionVector[47] = vD_PYR;
        reactionVector[48] = vD_Pck;
        reactionVector[49] = vD_Pfk;
        reactionVector[50] = vD_Ppc;
        reactionVector[51] = vD_Pps;
        reactionVector[52] = vD_Pyk;
        reactionVector[53] = vD_R5P;
        reactionVector[54] = vD_RU5P;
        reactionVector[55] = vD_S7P;
        reactionVector[56] = vD_SDH;
        reactionVector[57] = vD_SUC;
        reactionVector[58] = vD_X;
        reactionVector[59] = vD_X5P;
        reactionVector[60] = vD_aKG;
        reactionVector[61] = vD_aKGDH;
        reactionVector[62] = vD_cAMP;
        reactionVector[63] = vE_6PGDH;
        reactionVector[64] = vE_AceKki;
        reactionVector[65] = vE_AceKph;
        reactionVector[66] = vE_Ack;
        reactionVector[67] = vE_Ack_medium;
        reactionVector[68] = vE_Acs;
        reactionVector[69] = vE_Acs_medium;
        reactionVector[70] = vE_CS;
        reactionVector[71] = vE_Cya;
        reactionVector[72] = vE_Eda;
        reactionVector[73] = vE_Edd;
        reactionVector[74] = vE_Fba;
        reactionVector[75] = vE_Fbp;
        reactionVector[76] = vE_Fum;
        reactionVector[77] = vE_G6PDH;
        reactionVector[78] = vE_GAPDH;
        reactionVector[79] = vE_Glk;
        reactionVector[80] = vE_ICDH;
        reactionVector[81] = vE_Icl;
        reactionVector[82] = vE_MDH;
        reactionVector[83] = vE_MS;
        reactionVector[84] = vE_Mez;
        reactionVector[85] = vE_PDH;
        reactionVector[86] = vE_Pck;
        reactionVector[87] = vE_Pfk;
        reactionVector[88] = vE_Pgi;
        reactionVector[89] = vE_Pgl;
        reactionVector[90] = vE_Ppc;
        reactionVector[91] = vE_Pps;
        reactionVector[92] = vE_Pta;
        reactionVector[93] = vE_Pyk;
        reactionVector[94] = vE_R5PI;
        reactionVector[95] = vE_Ru5P;
        reactionVector[96] = vE_SDH;
        reactionVector[97] = vE_Tal;
        reactionVector[98] = vE_TktA;
        reactionVector[99] = vE_TktB;
        reactionVector[100] = vE_aKGDH;
        reactionVector[101] = vE_cAMPdegr;
        reactionVector[102] = vG_aceA;
        reactionVector[103] = vG_aceB;
        reactionVector[104] = vG_aceK;
        reactionVector[105] = vG_acs;
        reactionVector[106] = vG_fbaA;
        reactionVector[107] = vG_fbp;
        reactionVector[108] = vG_fumABC;
        reactionVector[109] = vG_gapA;
        reactionVector[110] = vG_glk;
        reactionVector[111] = vG_gltA;
        reactionVector[112] = vG_icdA;
        reactionVector[113] = vG_maeB;
        reactionVector[114] = vG_mdh;
        reactionVector[115] = vG_pckA;
        reactionVector[116] = vG_pdh;
        reactionVector[117] = vG_pfkA;
        reactionVector[118] = vG_ppc;
        reactionVector[119] = vG_ppsA;
        reactionVector[120] = vG_pykF;
        reactionVector[121] = vG_sdhCDAB;
        reactionVector[122] = vG_sucAB;
        reactionVector[123] = vNonPTS;
        reactionVector[124] = vNonPTS_medium;
        reactionVector[125] = vPTS1;
        reactionVector[126] = vPTS4;
        reactionVector[127] = vPTS4_medium;
        reactionVector[128] = vgrowth;
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double ACEex,AcCoA,AcP,AceK,Acs,CS,E4P,EIIAP,F6P,FBP,FUM,Fba,Fbp,Fum,G6P,GAP,GAPDH,GLC,GLCex,GLCfeed;
    double GOX,Glk,ICDH,ICDHP,ICIT,Icl,KDPG,MAL,MDH,MS,Mez,OAA,PDH,PEP,PYR,Pck,Pfk,Ppc,Pps,Pyk;
    double R5P,RU5P,S7P,SDH,SUC,X,X5P,aKG,aKGDH,cAMP,sixPG,sixPGL;
    double ADP,AMP,ATP,CoA,Cratotal,Crptotal,D,EIIAtotal,EXTERNAL,Factor_aceB,Factor_aceK,IclRtotal,K6PGDH_6PG,K6PGDH_ATP_inh,K6PGDH_NADP,K6PGDH_NADPH_inh,KAceK_3PG,KAceK_GOX,KAceK_ICDH,KAceK_ICDHP;
    double KAceK_ICIT,KAceK_OAA,KAceK_PEP,KAceK_PYR,KAceK_aKG,KAck_ACE_m,KAck_ADP_m,KAck_ATP_m,KAck_AcP_m,KAck_eq,KAcs_ACE,KCS_AcCoA,KCS_OAA,KCS_OAA_AcCoA,KCS_aKG,KCraFBP,KCrpcAMP,KCya_EIIAP,KEda_GAP_m,KEda_KDPG_m;
    double KEda_PYR_m,KEda_eq,KEdd_6PG_m,KEdd_KDPG_m,KEdd_eq,KFba_DHAP,KFba_FBP,KFba_GAP,KFba_GAP_inh,KFba_eq,KFbp_FBP,KFbp_PEP,KFum_FUM_m,KFum_eq,KG6PDH_G6P,KG6PDH_NADP,KG6PDH_NADPH_g6pinh,KG6PDH_NADPH_nadpinh,KGAPDH_GAP,KGAPDH_NAD;
    double KGAPDH_NADH,KGAPDH_PGP,KGAPDH_eq,KGlk_ATP_m,KGlk_G6P_i,KGlk_GLC_m,KICDH_ICIT,KICDH_PEP,KIcl_3PG,KIcl_ICIT,KIcl_PEP,KIcl_aKG,KMDH_MAL_I,KMDH_MAL_m,KMDH_NADH_I,KMDH_NADH_m,KMDH_NAD_I,KMDH_NAD_II,KMDH_NAD_m,KMDH_OAA_I;
    double KMDH_OAA_II,KMDH_OAA_m,KMDH_eq,KMS_AcCoA,KMS_GOX,KMS_GOX_AcCoA,KMez_AcCoA,KMez_MAL,KMez_cAMP,KNonPTS_I,KNonPTS_S,KPDH_AcCoA_m,KPDH_CoA_m,KPDH_NADH_m,KPDH_NAD_m,KPDH_PYR_m,KPDH_i,KPTS_EIIA,KPTS_GLC,KPck_ADP_i;
    double KPck_ATP_I,KPck_ATP_i,KPck_OAA,KPck_OAA_I,KPck_PEP,KPck_PEP_i,KPdhRPYR,KPfk_ADP_a,KPfk_ADP_b,KPfk_ADP_c,KPfk_AMP_a,KPfk_AMP_b,KPfk_ATP_s,KPfk_F6P_s,KPfk_PEP,KPgi_F6P,KPgi_F6P_6pginh,KPgi_G6P,KPgi_G6P_6pginh,KPgi_eq;
    double KPgl_6PGL_m,KPgl_6PG_m,KPgl_eq,KPgl_h1,KPgl_h2,KPpc_FBP,KPpc_PEP,KPps_PEP,KPps_PYR,KPta_AcCoA_i,KPta_AcP_i,KPta_AcP_m,KPta_CoA_i,KPta_Pi_i,KPta_Pi_m,KPta_eq,KPyk_ADP,KPyk_AMP,KPyk_ATP,KPyk_FBP;
    double KPyk_PEP,KR5PI_eq,KRu5P_eq,KSDH_SUC_m,KSDH_eq,KTal_eq,KTktA_eq,KTktB_eq,KaKGDH_CoA_m,KaKGDH_NADH_I,KaKGDH_NAD_m,KaKGDH_SUC_I,KaKGDH_Z,KaKGDH_aKG_I,KaKGDH_aKG_m,KaceBAK_Cra,KaceBAK_Crp,KaceBAK_DNA,KaceBAK_GOX,KaceBAK_PYR;
    double KaceBAK_PYRprime,Kacs_Crp,KcAMPdegr_cAMP,KfbaA_Cra,KfbaA_Crp,Kfbp_Cra,KfumABC_Crp,KgapA_Cra,KgapA_Crp,Kglk_Cra,KgltA_Crp,KicdA_Cra,Kmdh_Crp,KpckA_Cra,Kpdh_PdhR,KpfkA_Cra,KppsA_Cra,KpykF_Cra,KsdhCDAB_Crp,KsucAB_Crp;
    double LAceK,LFbp,LICDH,LIcl,LMez,LPfk,LPpc,LPps,LPyk,LaceBAK,NAD,NADH,NADP,NADPH,POratio,POratio_prime,PdhRtotal,Pi,SS_Mez_ACE,SS_Mez_GLC;
    double SS_Ppc_ACE,SS_Ppc_GLC,VFba_blf,aceBAK_DNA,kATP,kAceKki_cat,kAceKph_cat,kAcs_cat,kBM_ACE_AcCoA,kBM_ACE_E4P,kBM_ACE_F6P,kBM_ACE_FUM,kBM_ACE_G6P,kBM_ACE_GAP,kBM_ACE_OAA,kBM_ACE_PEP,kBM_ACE_PYR,kBM_ACE_R5P,kBM_ACE_SUC,kBM_ACE_aKG;
    double kBM_GLC_AcCoA,kBM_GLC_E4P,kBM_GLC_F6P,kBM_GLC_FUM,kBM_GLC_G6P,kBM_GLC_GAP,kBM_GLC_OAA,kBM_GLC_PEP,kBM_GLC_PYR,kBM_GLC_R5P,kBM_GLC_SUC,kBM_GLC_aKG,kCS_cat,kFba_cat,kFbp_cat,kFum1_cat,kFum2_cat,kGAPDH_cat,kGlk_cat,kICDH_cat;
    double kIcl_cat,kMDH1_cat,kMDH2_cat,kMS_cat,kMez_cat,kPDH_cat,kPTS1,kPck_cat,kPfk_cat,kPpc_cat,kPps_cat,kPyk_cat,kSDH1_cat,kSDH2_cat,kaKGDH_cat,kaceBAK_cat_IclR,kdegr,kexpr,kmPTS1,nAceK;
    double nCraFBP,nCrpcAMP,nFbp,nICDH,nIcl,nMez,nPdhRPYR,nPfk,nPpc,nPps,nPyk,nacs,nfumABC,ngltA,nsdhCDAB,nsucAB,pH,pH_Eda_m,pH_Edd_m,pK_Eda;
    double pK_Edd,rho,v6PGDH_max,vAck_max,vCya_max,vEda_max,vEdd_max,vG6PDH_max,vNonPTS_max,vPTS4_max,vPgi_max,vPgl_max,vPta_max,vR5PI_max,vRu5P_max,vTal_max,vTktA_max,vTktB_max,vaceBAK_Cra_bound,vaceBAK_Cra_unbound;
    double vaceBAK_Crp_bound,vaceBAK_Crp_unbound,vacs_Crp_bound,vacs_Crp_unbound,vcAMPdegr_max,vfbaA_Cra_bound,vfbaA_Cra_unbound,vfbaA_Crp_bound,vfbaA_Crp_unbound,vfbp_Cra_bound,vfbp_Cra_unbound,vfumABC_Crp_bound,vfumABC_Crp_unbound,vgapA_Cra_bound,vgapA_Cra_unbound,vgapA_Crp_bound,vgapA_Crp_unbound,vglk_Cra_bound,vglk_Cra_unbound,vgltA_Crp_bound;
    double vgltA_Crp_unbound,vicdA_Cra_bound,vicdA_Cra_unbound,vmdh_Crp_bound,vmdh_Crp_unbound,vpckA_Cra_bound,vpckA_Cra_unbound,vpdh_PdhR_bound,vpdh_PdhR_unbound,vpfkA_Cra_bound,vpfkA_Cra_unbound,vppsA_Cra_bound,vppsA_Cra_unbound,vpykF_Cra_bound,vpykF_Cra_unbound,vsdhCDAB_Crp_bound,vsdhCDAB_Crp_unbound,vsucAB_Crp_bound,vsucAB_Crp_unbound,default0;
    double CraFBP,B,A,Cra,CrpcAMP,Crp,H,EIIA,alphaGLC,alphaACE,SS_Ppc,SS_Mez,OP_NADH,OP_FADH2,PdhRPYR,PdhR,QEda_pH,QEdd_pH,vMDH1_max,vFum2_max;
    double vSDH2_max,vSDH1_max,vATP,mu,vMDH2_max,vFum1_max;
    ADP = paramdataPtr->parametervector[0]; /* 0.56 */
    AMP = paramdataPtr->parametervector[1]; /* 0.28 */
    ATP = paramdataPtr->parametervector[2]; /* 9.6 */
    CoA = paramdataPtr->parametervector[3]; /* 1.4 */
    Cratotal = paramdataPtr->parametervector[4]; /* 0.0003 */
    Crptotal = paramdataPtr->parametervector[5]; /* 0.0115 */
    D = paramdataPtr->parametervector[6]; /* 0 */
    EIIAtotal = paramdataPtr->parametervector[7]; /* 0.0769 */
    EXTERNAL = paramdataPtr->parametervector[8]; /* 0 */
    Factor_aceB = paramdataPtr->parametervector[9]; /* 0.509128 */
    Factor_aceK = paramdataPtr->parametervector[10]; /* 0.0293109 */
    IclRtotal = paramdataPtr->parametervector[11]; /* 8.3e-05 */
    K6PGDH_6PG = paramdataPtr->parametervector[12]; /* 0.0999956 */
    K6PGDH_ATP_inh = paramdataPtr->parametervector[13]; /* 3.03764 */
    K6PGDH_NADP = paramdataPtr->parametervector[14]; /* 0.00774325 */
    K6PGDH_NADPH_inh = paramdataPtr->parametervector[15]; /* 0.0180077 */
    KAceK_3PG = paramdataPtr->parametervector[16]; /* 1.0927 */
    KAceK_GOX = paramdataPtr->parametervector[17]; /* 0.415957 */
    KAceK_ICDH = paramdataPtr->parametervector[18]; /* 0.731925 */
    KAceK_ICDHP = paramdataPtr->parametervector[19]; /* 8.19069 */
    KAceK_ICIT = paramdataPtr->parametervector[20]; /* 0.0658771 */
    KAceK_OAA = paramdataPtr->parametervector[21]; /* 0.0345864 */
    KAceK_PEP = paramdataPtr->parametervector[22]; /* 0.232436 */
    KAceK_PYR = paramdataPtr->parametervector[23]; /* 0.0311645 */
    KAceK_aKG = paramdataPtr->parametervector[24]; /* 0.723859 */
    KAck_ACE_m = paramdataPtr->parametervector[25]; /* 6.33002 */
    KAck_ADP_m = paramdataPtr->parametervector[26]; /* 0.296006 */
    KAck_ATP_m = paramdataPtr->parametervector[27]; /* 0.0881512 */
    KAck_AcP_m = paramdataPtr->parametervector[28]; /* 0.0742685 */
    KAck_eq = paramdataPtr->parametervector[29]; /* 222.087 */
    KAcs_ACE = paramdataPtr->parametervector[30]; /* 0.0233726 */
    KCS_AcCoA = paramdataPtr->parametervector[31]; /* 0.172622 */
    KCS_OAA = paramdataPtr->parametervector[32]; /* 0.0127217 */
    KCS_OAA_AcCoA = paramdataPtr->parametervector[33]; /* 0.0177911 */
    KCS_aKG = paramdataPtr->parametervector[34]; /* 0.271462 */
    KCraFBP = paramdataPtr->parametervector[35]; /* 0.0322704 */
    KCrpcAMP = paramdataPtr->parametervector[36]; /* 0.368586 */
    KCya_EIIAP = paramdataPtr->parametervector[37]; /* 0.00242536 */
    KEda_GAP_m = paramdataPtr->parametervector[38]; /* 1.37717 */
    KEda_KDPG_m = paramdataPtr->parametervector[39]; /* 0.350113 */
    KEda_PYR_m = paramdataPtr->parametervector[40]; /* 11.6552 */
    KEda_eq = paramdataPtr->parametervector[41]; /* 0.500041 */
    KEdd_6PG_m = paramdataPtr->parametervector[42]; /* 0.347444 */
    KEdd_KDPG_m = paramdataPtr->parametervector[43]; /* 0.648821 */
    KEdd_eq = paramdataPtr->parametervector[44]; /* 1000.01 */
    KFba_DHAP = paramdataPtr->parametervector[45]; /* 0.0879915 */
    KFba_FBP = paramdataPtr->parametervector[46]; /* 0.13302 */
    KFba_GAP = paramdataPtr->parametervector[47]; /* 0.0879908 */
    KFba_GAP_inh = paramdataPtr->parametervector[48]; /* 0.771367 */
    KFba_eq = paramdataPtr->parametervector[49]; /* 0.33082 */
    KFbp_FBP = paramdataPtr->parametervector[50]; /* 0.00255028 */
    KFbp_PEP = paramdataPtr->parametervector[51]; /* 0.225362 */
    KFum_FUM_m = paramdataPtr->parametervector[52]; /* 0.0554577 */
    KFum_eq = paramdataPtr->parametervector[53]; /* 18.9482 */
    KG6PDH_G6P = paramdataPtr->parametervector[54]; /* 0.0699847 */
    KG6PDH_NADP = paramdataPtr->parametervector[55]; /* 0.00607658 */
    KG6PDH_NADPH_g6pinh = paramdataPtr->parametervector[56]; /* 0.192067 */
    KG6PDH_NADPH_nadpinh = paramdataPtr->parametervector[57]; /* 0.0251647 */
    KGAPDH_GAP = paramdataPtr->parametervector[58]; /* 0.104745 */
    KGAPDH_NAD = paramdataPtr->parametervector[59]; /* 0.449961 */
    KGAPDH_NADH = paramdataPtr->parametervector[60]; /* 0.0200034 */
    KGAPDH_PGP = paramdataPtr->parametervector[61]; /* 0.0995525 */
    KGAPDH_eq = paramdataPtr->parametervector[62]; /* 0.217878 */
    KGlk_ATP_m = paramdataPtr->parametervector[63]; /* 0.799718 */
    KGlk_G6P_i = paramdataPtr->parametervector[64]; /* 14.9942 */
    KGlk_GLC_m = paramdataPtr->parametervector[65]; /* 0.219961 */
    KICDH_ICIT = paramdataPtr->parametervector[66]; /* 9.55731e-05 */
    KICDH_PEP = paramdataPtr->parametervector[67]; /* 0.209919 */
    KIcl_3PG = paramdataPtr->parametervector[68]; /* 0.492465 */
    KIcl_ICIT = paramdataPtr->parametervector[69]; /* 0.0257752 */
    KIcl_PEP = paramdataPtr->parametervector[70]; /* 0.025628 */
    KIcl_aKG = paramdataPtr->parametervector[71]; /* 0.202044 */
    KMDH_MAL_I = paramdataPtr->parametervector[72]; /* 1.32753 */
    KMDH_MAL_m = paramdataPtr->parametervector[73]; /* 1.33028 */
    KMDH_NADH_I = paramdataPtr->parametervector[74]; /* 0.0242784 */
    KMDH_NADH_m = paramdataPtr->parametervector[75]; /* 0.0168477 */
    KMDH_NAD_I = paramdataPtr->parametervector[76]; /* 0.309937 */
    KMDH_NAD_II = paramdataPtr->parametervector[77]; /* 0.330748 */
    KMDH_NAD_m = paramdataPtr->parametervector[78]; /* 0.099998 */
    KMDH_OAA_I = paramdataPtr->parametervector[79]; /* 0.25207 */
    KMDH_OAA_II = paramdataPtr->parametervector[80]; /* 0.396105 */
    KMDH_OAA_m = paramdataPtr->parametervector[81]; /* 0.269873 */
    KMDH_eq = paramdataPtr->parametervector[82]; /* 0.557468 */
    KMS_AcCoA = paramdataPtr->parametervector[83]; /* 0.381655 */
    KMS_GOX = paramdataPtr->parametervector[84]; /* 0.361087 */
    KMS_GOX_AcCoA = paramdataPtr->parametervector[85]; /* 0.311119 */
    KMez_AcCoA = paramdataPtr->parametervector[86]; /* 2.99104 */
    KMez_MAL = paramdataPtr->parametervector[87]; /* 0.00293864 */
    KMez_cAMP = paramdataPtr->parametervector[88]; /* 2.40366 */
    KNonPTS_I = paramdataPtr->parametervector[89]; /* 0.009166 */
    KNonPTS_S = paramdataPtr->parametervector[90]; /* 2.34579 */
    KPDH_AcCoA_m = paramdataPtr->parametervector[91]; /* 0.0080012 */
    KPDH_CoA_m = paramdataPtr->parametervector[92]; /* 0.00698447 */
    KPDH_NADH_m = paramdataPtr->parametervector[93]; /* 0.248943 */
    KPDH_NAD_m = paramdataPtr->parametervector[94]; /* 0.399784 */
    KPDH_PYR_m = paramdataPtr->parametervector[95]; /* 1.20561 */
    KPDH_i = paramdataPtr->parametervector[96]; /* 21.7028 */
    KPTS_EIIA = paramdataPtr->parametervector[97]; /* 0.00660253 */
    KPTS_GLC = paramdataPtr->parametervector[98]; /* 0.0127905 */
    KPck_ADP_i = paramdataPtr->parametervector[99]; /* 0.0388343 */
    KPck_ATP_I = paramdataPtr->parametervector[100]; /* 0.0396311 */
    KPck_ATP_i = paramdataPtr->parametervector[101]; /* 0.0415227 */
    KPck_OAA = paramdataPtr->parametervector[102]; /* 0.669902 */
    KPck_OAA_I = paramdataPtr->parametervector[103]; /* 0.532385 */
    KPck_PEP = paramdataPtr->parametervector[104]; /* 0.0701171 */
    KPck_PEP_i = paramdataPtr->parametervector[105]; /* 0.059769 */
    KPdhRPYR = paramdataPtr->parametervector[106]; /* 0.0856307 */
    KPfk_ADP_a = paramdataPtr->parametervector[107]; /* 238.935 */
    KPfk_ADP_b = paramdataPtr->parametervector[108]; /* 0.250106 */
    KPfk_ADP_c = paramdataPtr->parametervector[109]; /* 0.359964 */
    KPfk_AMP_a = paramdataPtr->parametervector[110]; /* 6.05077 */
    KPfk_AMP_b = paramdataPtr->parametervector[111]; /* 0.0207231 */
    KPfk_ATP_s = paramdataPtr->parametervector[112]; /* 0.160031 */
    KPfk_F6P_s = paramdataPtr->parametervector[113]; /* 0.0231651 */
    KPfk_PEP = paramdataPtr->parametervector[114]; /* 1.93645 */
    KPgi_F6P = paramdataPtr->parametervector[115]; /* 0.199983 */
    KPgi_F6P_6pginh = paramdataPtr->parametervector[116]; /* 0.19997 */
    KPgi_G6P = paramdataPtr->parametervector[117]; /* 2.45964 */
    KPgi_G6P_6pginh = paramdataPtr->parametervector[118]; /* 0.200129 */
    KPgi_eq = paramdataPtr->parametervector[119]; /* 0.717321 */
    KPgl_6PGL_m = paramdataPtr->parametervector[120]; /* 0.0229252 */
    KPgl_6PG_m = paramdataPtr->parametervector[121]; /* 9.96181 */
    KPgl_eq = paramdataPtr->parametervector[122]; /* 42.7993 */
    KPgl_h1 = paramdataPtr->parametervector[123]; /* 0.00225241 */
    KPgl_h2 = paramdataPtr->parametervector[124]; /* 9.71558e-06 */
    KPpc_FBP = paramdataPtr->parametervector[125]; /* 0.126971 */
    KPpc_PEP = paramdataPtr->parametervector[126]; /* 0.0385597 */
    KPps_PEP = paramdataPtr->parametervector[127]; /* 0.000604094 */
    KPps_PYR = paramdataPtr->parametervector[128]; /* 0.0006405 */
    KPta_AcCoA_i = paramdataPtr->parametervector[129]; /* 0.088354 */
    KPta_AcP_i = paramdataPtr->parametervector[130]; /* 0.0678373 */
    KPta_AcP_m = paramdataPtr->parametervector[131]; /* 0.523402 */
    KPta_CoA_i = paramdataPtr->parametervector[132]; /* 0.0244299 */
    KPta_Pi_i = paramdataPtr->parametervector[133]; /* 1.61647 */
    KPta_Pi_m = paramdataPtr->parametervector[134]; /* 0.699177 */
    KPta_eq = paramdataPtr->parametervector[135]; /* 0.0326788 */
    KPyk_ADP = paramdataPtr->parametervector[136]; /* 0.213774 */
    KPyk_AMP = paramdataPtr->parametervector[137]; /* 0.209251 */
    KPyk_ATP = paramdataPtr->parametervector[138]; /* 22.4085 */
    KPyk_FBP = paramdataPtr->parametervector[139]; /* 0.190009 */
    KPyk_PEP = paramdataPtr->parametervector[140]; /* 0.370704 */
    KR5PI_eq = paramdataPtr->parametervector[141]; /* 0.484213 */
    KRu5P_eq = paramdataPtr->parametervector[142]; /* 1.62946 */
    KSDH_SUC_m = paramdataPtr->parametervector[143]; /* 0.220018 */
    KSDH_eq = paramdataPtr->parametervector[144]; /* 8.67636 */
    KTal_eq = paramdataPtr->parametervector[145]; /* 1.17897 */
    KTktA_eq = paramdataPtr->parametervector[146]; /* 1.2001 */
    KTktB_eq = paramdataPtr->parametervector[147]; /* 10 */
    KaKGDH_CoA_m = paramdataPtr->parametervector[148]; /* 0.00391074 */
    KaKGDH_NADH_I = paramdataPtr->parametervector[149]; /* 0.0180044 */
    KaKGDH_NAD_m = paramdataPtr->parametervector[150]; /* 0.0700101 */
    KaKGDH_SUC_I = paramdataPtr->parametervector[151]; /* 0.533941 */
    KaKGDH_Z = paramdataPtr->parametervector[152]; /* 1.41777 */
    KaKGDH_aKG_I = paramdataPtr->parametervector[153]; /* 1.78176 */
    KaKGDH_aKG_m = paramdataPtr->parametervector[154]; /* 0.99931 */
    KaceBAK_Cra = paramdataPtr->parametervector[155]; /* 0.00018881 */
    KaceBAK_Crp = paramdataPtr->parametervector[156]; /* 4.56257 */
    KaceBAK_DNA = paramdataPtr->parametervector[157]; /* 1.13147e-06 */
    KaceBAK_GOX = paramdataPtr->parametervector[158]; /* 0.00566376 */
    KaceBAK_PYR = paramdataPtr->parametervector[159]; /* 0.249242 */
    KaceBAK_PYRprime = paramdataPtr->parametervector[160]; /* 0.00858996 */
    Kacs_Crp = paramdataPtr->parametervector[161]; /* 0.00846341 */
    KcAMPdegr_cAMP = paramdataPtr->parametervector[162]; /* 0.0662446 */
    KfbaA_Cra = paramdataPtr->parametervector[163]; /* 0.00373129 */
    KfbaA_Crp = paramdataPtr->parametervector[164]; /* 0.03915 */
    Kfbp_Cra = paramdataPtr->parametervector[165]; /* 0.000122915 */
    KfumABC_Crp = paramdataPtr->parametervector[166]; /* 0.122193 */
    KgapA_Cra = paramdataPtr->parametervector[167]; /* 0.00494581 */
    KgapA_Crp = paramdataPtr->parametervector[168]; /* 0.0296522 */
    Kglk_Cra = paramdataPtr->parametervector[169]; /* 2.91919e-08 */
    KgltA_Crp = paramdataPtr->parametervector[170]; /* 0.0312059 */
    KicdA_Cra = paramdataPtr->parametervector[171]; /* 5.3399e-05 */
    Kmdh_Crp = paramdataPtr->parametervector[172]; /* 0.255211 */
    KpckA_Cra = paramdataPtr->parametervector[173]; /* 0.000275018 */
    Kpdh_PdhR = paramdataPtr->parametervector[174]; /* 1.24734e-05 */
    KpfkA_Cra = paramdataPtr->parametervector[175]; /* 2.41644e-08 */
    KppsA_Cra = paramdataPtr->parametervector[176]; /* 0.00062288 */
    KpykF_Cra = paramdataPtr->parametervector[177]; /* 0.000117808 */
    KsdhCDAB_Crp = paramdataPtr->parametervector[178]; /* 0.045294 */
    KsucAB_Crp = paramdataPtr->parametervector[179]; /* 0.117486 */
    LAceK = paramdataPtr->parametervector[180]; /* 7.59785e+07 */
    LFbp = paramdataPtr->parametervector[181]; /* 1.11817e+07 */
    LICDH = paramdataPtr->parametervector[182]; /* 126.46 */
    LIcl = paramdataPtr->parametervector[183]; /* 33725.8 */
    LMez = paramdataPtr->parametervector[184]; /* 197779 */
    LPfk = paramdataPtr->parametervector[185]; /* 1.26952e+06 */
    LPpc = paramdataPtr->parametervector[186]; /* 3.8947e+06 */
    LPps = paramdataPtr->parametervector[187]; /* 4.03223e-80 */
    LPyk = paramdataPtr->parametervector[188]; /* 1000.15 */
    LaceBAK = paramdataPtr->parametervector[189]; /* 1243.52 */
    NAD = paramdataPtr->parametervector[190]; /* 2.6 */
    NADH = paramdataPtr->parametervector[191]; /* 0.083 */
    NADP = paramdataPtr->parametervector[192]; /* 0.0021 */
    NADPH = paramdataPtr->parametervector[193]; /* 0.12 */
    POratio = paramdataPtr->parametervector[194]; /* 3.30385 */
    POratio_prime = paramdataPtr->parametervector[195]; /* 1.50583 */
    PdhRtotal = paramdataPtr->parametervector[196]; /* 6.66e-05 */
    Pi = paramdataPtr->parametervector[197]; /* 10 */
    SS_Mez_ACE = paramdataPtr->parametervector[198]; /* 0.0334723 */
    SS_Mez_GLC = paramdataPtr->parametervector[199]; /* 0.0133738 */
    SS_Ppc_ACE = paramdataPtr->parametervector[200]; /* 0.000852337 */
    SS_Ppc_GLC = paramdataPtr->parametervector[201]; /* 0.00323703 */
    VFba_blf = paramdataPtr->parametervector[202]; /* 2.00107 */
    aceBAK_DNA = paramdataPtr->parametervector[203]; /* 5.15e-07 */
    kATP = paramdataPtr->parametervector[204]; /* 1.30699e-05 */
    kAceKki_cat = paramdataPtr->parametervector[205]; /* 9.19676e+15 */
    kAceKph_cat = paramdataPtr->parametervector[206]; /* 8.85104e+12 */
    kAcs_cat = paramdataPtr->parametervector[207]; /* 116474 */
    kBM_ACE_AcCoA = paramdataPtr->parametervector[208]; /* 164.317 */
    kBM_ACE_E4P = paramdataPtr->parametervector[209]; /* 284.114 */
    kBM_ACE_F6P = paramdataPtr->parametervector[210]; /* 321.291 */
    kBM_ACE_FUM = paramdataPtr->parametervector[211]; /* 143.731 */
    kBM_ACE_G6P = paramdataPtr->parametervector[212]; /* 118.486 */
    kBM_ACE_GAP = paramdataPtr->parametervector[213]; /* 734.492 */
    kBM_ACE_OAA = paramdataPtr->parametervector[214]; /* 7011.8 */
    kBM_ACE_PEP = paramdataPtr->parametervector[215]; /* 202.218 */
    kBM_ACE_PYR = paramdataPtr->parametervector[216]; /* 12876.3 */
    kBM_ACE_R5P = paramdataPtr->parametervector[217]; /* 240.315 */
    kBM_ACE_SUC = paramdataPtr->parametervector[218]; /* 150.808 */
    kBM_ACE_aKG = paramdataPtr->parametervector[219]; /* 324.664 */
    kBM_GLC_AcCoA = paramdataPtr->parametervector[220]; /* 2656.7 */
    kBM_GLC_E4P = paramdataPtr->parametervector[221]; /* 603.434 */
    kBM_GLC_F6P = paramdataPtr->parametervector[222]; /* 966.423 */
    kBM_GLC_FUM = paramdataPtr->parametervector[223]; /* 1091.82 */
    kBM_GLC_G6P = paramdataPtr->parametervector[224]; /* 52.0836 */
    kBM_GLC_GAP = paramdataPtr->parametervector[225]; /* 375.616 */
    kBM_GLC_OAA = paramdataPtr->parametervector[226]; /* 19606.2 */
    kBM_GLC_PEP = paramdataPtr->parametervector[227]; /* 708.44 */
    kBM_GLC_PYR = paramdataPtr->parametervector[228]; /* 707.651 */
    kBM_GLC_R5P = paramdataPtr->parametervector[229]; /* 247.339 */
    kBM_GLC_SUC = paramdataPtr->parametervector[230]; /* 2467.94 */
    kBM_GLC_aKG = paramdataPtr->parametervector[231]; /* 4673.85 */
    kCS_cat = paramdataPtr->parametervector[232]; /* 792390 */
    kFba_cat = paramdataPtr->parametervector[233]; /* 1.60675e+06 */
    kFbp_cat = paramdataPtr->parametervector[234]; /* 2.41475e+06 */
    kFum1_cat = paramdataPtr->parametervector[235]; /* 681983 */
    kFum2_cat = paramdataPtr->parametervector[236]; /* 681983 */
    kGAPDH_cat = paramdataPtr->parametervector[237]; /* 8.79789e+07 */
    kGlk_cat = paramdataPtr->parametervector[238]; /* 1.98397e+06 */
    kICDH_cat = paramdataPtr->parametervector[239]; /* 211215 */
    kIcl_cat = paramdataPtr->parametervector[240]; /* 1.09342e+06 */
    kMDH1_cat = paramdataPtr->parametervector[241]; /* 328277 */
    kMDH2_cat = paramdataPtr->parametervector[242]; /* 328277 */
    kMS_cat = paramdataPtr->parametervector[243]; /* 13791.5 */
    kMez_cat = paramdataPtr->parametervector[244]; /* 1.65343e+06 */
    kPDH_cat = paramdataPtr->parametervector[245]; /* 6.59379e+07 */
    kPTS1 = paramdataPtr->parametervector[246]; /* 42555.9 */
    kPck_cat = paramdataPtr->parametervector[247]; /* 2.60105e+06 */
    kPfk_cat = paramdataPtr->parametervector[248]; /* 1.49637e+10 */
    kPpc_cat = paramdataPtr->parametervector[249]; /* 1.0471e+07 */
    kPps_cat = paramdataPtr->parametervector[250]; /* 1070.78 */
    kPyk_cat = paramdataPtr->parametervector[251]; /* 17530.3 */
    kSDH1_cat = paramdataPtr->parametervector[252]; /* 22311.9 */
    kSDH2_cat = paramdataPtr->parametervector[253]; /* 22311.9 */
    kaKGDH_cat = paramdataPtr->parametervector[254]; /* 2.04051e+08 */
    kaceBAK_cat_IclR = paramdataPtr->parametervector[255]; /* 3.36933 */
    kdegr = paramdataPtr->parametervector[256]; /* 0.265324 */
    kexpr = paramdataPtr->parametervector[257]; /* 4.85939 */
    kmPTS1 = paramdataPtr->parametervector[258]; /* 12994.3 */
    nAceK = paramdataPtr->parametervector[259]; /* 2 */
    nCraFBP = paramdataPtr->parametervector[260]; /* 2 */
    nCrpcAMP = paramdataPtr->parametervector[261]; /* 1 */
    nFbp = paramdataPtr->parametervector[262]; /* 4 */
    nICDH = paramdataPtr->parametervector[263]; /* 2 */
    nIcl = paramdataPtr->parametervector[264]; /* 4 */
    nMez = paramdataPtr->parametervector[265]; /* 1.92774 */
    nPdhRPYR = paramdataPtr->parametervector[266]; /* 1 */
    nPfk = paramdataPtr->parametervector[267]; /* 4 */
    nPpc = paramdataPtr->parametervector[268]; /* 3 */
    nPps = paramdataPtr->parametervector[269]; /* 2 */
    nPyk = paramdataPtr->parametervector[270]; /* 4 */
    nacs = paramdataPtr->parametervector[271]; /* 2.31 */
    nfumABC = paramdataPtr->parametervector[272]; /* 0.74 */
    ngltA = paramdataPtr->parametervector[273]; /* 1.07 */
    nsdhCDAB = paramdataPtr->parametervector[274]; /* 0.74 */
    nsucAB = paramdataPtr->parametervector[275]; /* 0.74 */
    pH = paramdataPtr->parametervector[276]; /* 7.5 */
    pH_Eda_m = paramdataPtr->parametervector[277]; /* 7.49978 */
    pH_Edd_m = paramdataPtr->parametervector[278]; /* 7.15931 */
    pK_Eda = paramdataPtr->parametervector[279]; /* 13.1433 */
    pK_Edd = paramdataPtr->parametervector[280]; /* 3.42317 */
    rho = paramdataPtr->parametervector[281]; /* 564 */
    v6PGDH_max = paramdataPtr->parametervector[282]; /* 193585 */
    vAck_max = paramdataPtr->parametervector[283]; /* 533045 */
    vCya_max = paramdataPtr->parametervector[284]; /* 13.2427 */
    vEda_max = paramdataPtr->parametervector[285]; /* 582.202 */
    vEdd_max = paramdataPtr->parametervector[286]; /* 944.737 */
    vG6PDH_max = paramdataPtr->parametervector[287]; /* 53780.5 */
    vNonPTS_max = paramdataPtr->parametervector[288]; /* 4679.47 */
    vPTS4_max = paramdataPtr->parametervector[289]; /* 10119.9 */
    vPgi_max = paramdataPtr->parametervector[290]; /* 3.62077e+06 */
    vPgl_max = paramdataPtr->parametervector[291]; /* 45257.5 */
    vPta_max = paramdataPtr->parametervector[292]; /* 5915.4 */
    vR5PI_max = paramdataPtr->parametervector[293]; /* 42445.8 */
    vRu5P_max = paramdataPtr->parametervector[294]; /* 23336.1 */
    vTal_max = paramdataPtr->parametervector[295]; /* 42057.1 */
    vTktA_max = paramdataPtr->parametervector[296]; /* 7656.02 */
    vTktB_max = paramdataPtr->parametervector[297]; /* 285535 */
    vaceBAK_Cra_bound = paramdataPtr->parametervector[298]; /* 0.0596534 */
    vaceBAK_Cra_unbound = paramdataPtr->parametervector[299]; /* 7.53474e-05 */
    vaceBAK_Crp_bound = paramdataPtr->parametervector[300]; /* 1.59563e-05 */
    vaceBAK_Crp_unbound = paramdataPtr->parametervector[301]; /* 0.00262305 */
    vacs_Crp_bound = paramdataPtr->parametervector[302]; /* 0.000592461 */
    vacs_Crp_unbound = paramdataPtr->parametervector[303]; /* 0 */
    vcAMPdegr_max = paramdataPtr->parametervector[304]; /* 9.91179 */
    vfbaA_Cra_bound = paramdataPtr->parametervector[305]; /* 0 */
    vfbaA_Cra_unbound = paramdataPtr->parametervector[306]; /* 0.0173187 */
    vfbaA_Crp_bound = paramdataPtr->parametervector[307]; /* 0.0150909 */
    vfbaA_Crp_unbound = paramdataPtr->parametervector[308]; /* 0 */
    vfbp_Cra_bound = paramdataPtr->parametervector[309]; /* 0.000701946 */
    vfbp_Cra_unbound = paramdataPtr->parametervector[310]; /* 0 */
    vfumABC_Crp_bound = paramdataPtr->parametervector[311]; /* 0.0540414 */
    vfumABC_Crp_unbound = paramdataPtr->parametervector[312]; /* 0 */
    vgapA_Cra_bound = paramdataPtr->parametervector[313]; /* 0 */
    vgapA_Cra_unbound = paramdataPtr->parametervector[314]; /* 0.0156125 */
    vgapA_Crp_bound = paramdataPtr->parametervector[315]; /* 0.0242728 */
    vgapA_Crp_unbound = paramdataPtr->parametervector[316]; /* 0 */
    vglk_Cra_bound = paramdataPtr->parametervector[317]; /* 0.00139121 */
    vglk_Cra_unbound = paramdataPtr->parametervector[318]; /* 0.184677 */
    vgltA_Crp_bound = paramdataPtr->parametervector[319]; /* 0.0505062 */
    vgltA_Crp_unbound = paramdataPtr->parametervector[320]; /* 0 */
    vicdA_Cra_bound = paramdataPtr->parametervector[321]; /* 0.0159016 */
    vicdA_Cra_unbound = paramdataPtr->parametervector[322]; /* 0.00401129 */
    vmdh_Crp_bound = paramdataPtr->parametervector[323]; /* 0.158286 */
    vmdh_Crp_unbound = paramdataPtr->parametervector[324]; /* 0 */
    vpckA_Cra_bound = paramdataPtr->parametervector[325]; /* 0.0110092 */
    vpckA_Cra_unbound = paramdataPtr->parametervector[326]; /* 0 */
    vpdh_PdhR_bound = paramdataPtr->parametervector[327]; /* 9.45249e-06 */
    vpdh_PdhR_unbound = paramdataPtr->parametervector[328]; /* 0.00156527 */
    vpfkA_Cra_bound = paramdataPtr->parametervector[329]; /* 0.00272437 */
    vpfkA_Cra_unbound = paramdataPtr->parametervector[330]; /* 0.0754715 */
    vppsA_Cra_bound = paramdataPtr->parametervector[331]; /* 0.113449 */
    vppsA_Cra_unbound = paramdataPtr->parametervector[332]; /* 0 */
    vpykF_Cra_bound = paramdataPtr->parametervector[333]; /* 8.98315e-05 */
    vpykF_Cra_unbound = paramdataPtr->parametervector[334]; /* 0.0038128 */
    vsdhCDAB_Crp_bound = paramdataPtr->parametervector[335]; /* 0.357936 */
    vsdhCDAB_Crp_unbound = paramdataPtr->parametervector[336]; /* 0 */
    vsucAB_Crp_bound = paramdataPtr->parametervector[337]; /* 0.0282708 */
    vsucAB_Crp_unbound = paramdataPtr->parametervector[338]; /* 0 */
    default0 = paramdataPtr->parametervector[339]; /* 1 */
    CraFBP = Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP));
    B = 1.0+ADP/KPfk_ADP_a+AMP/KPfk_AMP_a;
    A = 1.0+PEP/KPfk_PEP+ADP/KPfk_ADP_b+AMP/KPfk_AMP_b;
    Cra = Cratotal-CraFBP;
    CrpcAMP = Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP));
    Crp = Crptotal-CrpcAMP;
    H = pow(10.0,-pH)*1000.0;
    EIIA = EIIAtotal-EIIAP;
    alphaGLC = GLCex/(GLCex+KPTS_GLC);
    alphaACE = ACEex/(ACEex+KAcs_ACE)*(1.0-alphaGLC);
    SS_Ppc = alphaGLC*SS_Ppc_GLC+alphaACE*SS_Ppc_ACE;
    SS_Mez = alphaGLC*SS_Mez_GLC+alphaACE*SS_Mez_ACE;
    OP_NADH = (vE_GAPDH+vE_PDH+vE_aKGDH+vE_MDH)*POratio;
    OP_FADH2 = vE_SDH*POratio_prime;
    PdhRPYR = PdhRtotal*pow(PYR,nPdhRPYR)/(pow(PYR,nPdhRPYR)+pow(KPdhRPYR,nPdhRPYR));
    PdhR = PdhRtotal-PdhRPYR;
    QEda_pH = 1.0+2.0*pow(10.0,pH_Eda_m-pK_Eda)/(1.0+pow(10.0,pH-pK_Eda)+pow(10.0,2.0*pH_Eda_m-pH-pK_Eda));
    QEdd_pH = 1.0+2.0*pow(10.0,pH_Edd_m-pK_Edd)/(1.0+pow(10.0,pH-pK_Edd)+pow(10.0,2.0*pH_Edd_m-pH-pK_Edd));
    vMDH1_max = MDH*kMDH1_cat;
    vFum2_max = Fum*kFum2_cat;
    vSDH2_max = SDH*kSDH2_cat;
    vSDH1_max = SDH*kSDH1_cat;
    vATP = OP_NADH+OP_FADH2-vE_Glk-vE_Pfk+vE_GAPDH+vE_Pyk-vE_Pps+vE_Ack-vE_Acs+vE_aKGDH-vE_Pck-vE_AceKki-vE_Cya;
    mu = kATP*vATP;
    vMDH2_max = MDH*kMDH2_cat;
    vFum1_max = Fum*kFum1_cat;
    ACEex = 1.11524e-10;
    AcCoA = 0.422054;
    AcP = 0.0237175;
    AceK = 0.00107739;
    Acs = 0.000406765;
    CS = 0.0200874;
    E4P = 0.0261385;
    EIIAP = 0.00367648;
    F6P = 0.650133;
    FBP = 0.0992524;
    FUM = 0.0372591;
    Fba = 0.0651811;
    Fbp = 0.000460524;
    Fum = 0.0152583;
    G6P = 0.913391;
    GAP = 0.173332;
    GAPDH = 0.0652304;
    GLC = 0.012;
    GLCex = 22.2222;
    GLCfeed = 0.0;
    GOX = 3.77198e-10;
    Glk = 0.00548198;
    ICDH = 0.0154742;
    ICDHP = 0.0128716;
    ICIT = 0.00736962;
    Icl = 0.0367581;
    KDPG = 1.38744;
    MAL = 0.417978;
    MDH = 0.00976235;
    MS = 0.0187144;
    Mez = 0.0133654;
    OAA = 0.00254143;
    PDH = 0.00183994;
    PEP = 1.16312;
    PYR = 0.144336;
    Pck = 0.00360495;
    Pfk = 0.00967106;
    Ppc = 0.00323499;
    Pps = 0.0173145;
    Pyk = 0.0107065;
    R5P = 0.156479;
    RU5P = 0.334105;
    S7P = 0.141549;
    SDH = 0.193564;
    SUC = 0.116135;
    X = 0.011;
    X5P = 0.486039;
    aKG = 0.0010271;
    aKGDH = 0.00819793;
    cAMP = 0.246959;
    sixPG = 0.391757;
    sixPGL = 0.0109159;
    icVector[0] = ACEex;
    icVector[1] = AcCoA;
    icVector[2] = AcP;
    icVector[3] = AceK;
    icVector[4] = Acs;
    icVector[5] = CS;
    icVector[6] = E4P;
    icVector[7] = EIIAP;
    icVector[8] = F6P;
    icVector[9] = FBP;
    icVector[10] = FUM;
    icVector[11] = Fba;
    icVector[12] = Fbp;
    icVector[13] = Fum;
    icVector[14] = G6P;
    icVector[15] = GAP;
    icVector[16] = GAPDH;
    icVector[17] = GLC;
    icVector[18] = GLCex;
    icVector[19] = GLCfeed;
    icVector[20] = GOX;
    icVector[21] = Glk;
    icVector[22] = ICDH;
    icVector[23] = ICDHP;
    icVector[24] = ICIT;
    icVector[25] = Icl;
    icVector[26] = KDPG;
    icVector[27] = MAL;
    icVector[28] = MDH;
    icVector[29] = MS;
    icVector[30] = Mez;
    icVector[31] = OAA;
    icVector[32] = PDH;
    icVector[33] = PEP;
    icVector[34] = PYR;
    icVector[35] = Pck;
    icVector[36] = Pfk;
    icVector[37] = Ppc;
    icVector[38] = Pps;
    icVector[39] = Pyk;
    icVector[40] = R5P;
    icVector[41] = RU5P;
    icVector[42] = S7P;
    icVector[43] = SDH;
    icVector[44] = SUC;
    icVector[45] = X;
    icVector[46] = X5P;
    icVector[47] = aKG;
    icVector[48] = aKGDH;
    icVector[49] = cAMP;
    icVector[50] = sixPG;
    icVector[51] = sixPGL;
}

