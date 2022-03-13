/* This is the main file in this program. */
/* It orchestrates everything else. */

#include "con2prim.h"
#include "beta_equilibrium.c"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define SQR(x) ((x) * (x))

int n_tab_beta;
double Y_e_tab[MAX_NTAB], log_rho0_tab_beta[MAX_NTAB];
int n_nearest_beta;

void prepare_c2p(CCTK_ARGUMENTS){

  printf("\n ~~~ INSIDE prepare_c2p ~~~ \n");
  /* Here we prepare the EOS format needed for con2prim. */
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  c2p.eoskey = EOS_Omni_GetHandle_(c2p_EOS_type);
  c2p.key_temp = c2p_rhoT_key;
  *c2p_eos_eoskey = c2p.eoskey;

  // general EOS settings
  c2p.eoskey_polytrope = 1;
  c2p.eos_prec = 1e-10;
//  c2p.eos_prec = c2p_precEOS;

  // tolerance for convergence in recovery scheme
  c2p.tol_x = c2p_tolerance;
  c2p.tol_x_retry = c2p_tolerance_retry; // reduced tolerance for retry if first attempt failed

  // maximum number of iterations for c2p scheme
  c2p.max_iterations   = c2p_max_iterations;
  c2p.extra_iterations = c2p_extra_iterations;

  // further EOS settings
  double rhomi,rhoma,tmi,tma,yemi,yema,emi,ema,pmi,pma;
  tmi = 0.0;
  tma = 200.0;
  yemi = 0.0;
  yema = 1.0;
  rhomi = 0.0;
  rhoma = 1.0e16;
  emi = 0.0;
  ema = 1.0e50;
  pmi = 0.0;
  pma = 1.0e50;

  if (c2p.eoskey == 4){
    EOS_nuc_eos_get_lims(&rhomi,&rhoma,&tmi,&tma,&yemi,&yema,&emi,&ema,&pmi,&pma);
  }

  c2p.eos_rho_min = rhomi;
  c2p.eos_rho_max = rhoma;
  c2p.eos_temp_min = tmi;
  c2p.eos_temp_max = tma;
  c2p.eos_eps_min = emi; //-1.518793e-03; //emi;
  c2p.eos_eps_max = ema;
  c2p.eos_press_min = pmi;
  c2p.eos_press_max = pma;
  c2p.eos_ye_min = yemi;
  c2p.eos_ye_max = yema;

  // atmosphere settings
  c2p.rho_atmo_tol = c2p_rho_atmo_tolerance; // relative tolerance
  c2p.rho_atmo = (c2p_rho_atmo > c2p.eos_rho_min) ? c2p_rho_atmo : c2p.eos_rho_min;    // in M_Sun = c = G = 1
  c2p.T_atmo   = (c2p_T_atmo > c2p.eos_temp_min)  ? c2p_T_atmo   : c2p.eos_temp_min;  // in MeV, 8.62e-06 MeV = 1e5K
//  c2p.Ye_atmo  = (c2p_ye_atmo > c2p.eos_ye_min)   ? c2p_ye_atmo  : c2p.eos_ye_min;   // electron fraction
  c2p.Ye_atmo  = c2p_ye_atmo;   // electron fraction
  c2p.retain_B_atmo = c2p_retain_B_atmo; // keep the magnetic field in the atmo?


  /* Compute atmospheric pressure */
  CCTK_REAL xrho = c2p.rho_atmo;
  CCTK_REAL xeps, xpress;
  CCTK_REAL xtemp  = c2p.T_atmo;
  CCTK_REAL xye    = c2p.Ye_atmo;
  CCTK_INT  anyerr, keyerr;
  CCTK_INT  havetemp = 1;
  CCTK_INT  npoints  = 1;
  EOS_Omni_press(c2p.eoskey, havetemp, c2p.eos_prec, npoints, &xrho, &xeps, &xtemp, &xye, &xpress, &keyerr, &anyerr);
  c2p.P_atmo = xpress;
  c2p.eps_atmo = xeps;
  if(c2p_use_eps_atmo){c2p.eps_atmo = c2p_eps_atmo;}
//  EOS_Omni_press_from_rhoenthalpy(eoskey, keytemp, precEOS, 1, &xrho, &xeps, &xtemp, &xye, &xprs, &cH, &anyerr, &keyerr);

  /* Sharing with other thorns */
  *c2p_eos_prec     = c2p.eos_prec;
  *c2p_eos_P_atmo   = c2p.P_atmo;
  *c2p_eos_ye_atmo  = c2p.Ye_atmo;
  *c2p_eos_T_atmo   = c2p.T_atmo;
  *c2p_eos_eps_atmo = c2p.eps_atmo;
  *c2p_eos_rho_atmo = c2p.rho_atmo;

  /* Sharing with other thorns - INTS*/
  *c2p_eos_eoskey = c2p.eoskey;
  *c2p_eos_keytemp = c2p.key_temp;
  // number of primitive and conserved quantities
  c2p.numprims = 14; // rho, vel123, eps, B123, ye, temp, press, W, ent, abar
  c2p.numcons = 10; // D, S123, tau, B123, ye_con, temp

  // stricly enforce v^2 < 1
  if (c2p_enforce_v2_smaller_than_1)
    c2p.enforce_v2 = 1;
  else
    c2p.enforce_v2 = 0;

  // set recovery method
  if      (CCTK_EQUALS(c2p_preferred_algorithm, "auto"))        c2p.c2p_method = 5;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "NR_3DE"))      c2p.c2p_method = 1;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "NR_3DP"))      c2p.c2p_method = 2;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "NR_2DR"))      c2p.c2p_method = 3;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "NR_2DW"))      c2p.c2p_method = 4;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "palenzuela"))  c2p.c2p_method = 5;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "NR_2D_Noble")) c2p.c2p_method = 6;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "Newman"))      c2p.c2p_method = 7;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "HARM")) c2p.c2p_method = 11;
  else CCTK_PARAMWARN("The specified solution method does not exist or is not implemented!");

  c2p.use_c2p_method_backup = c2p_use_backup_scheme;
  if      (CCTK_EQUALS(c2p_algorithm_retry, "auto"))        c2p.c2p_method_backup = 5;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "NR_3DE"))      c2p.c2p_method_backup = 1;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "NR_3DP"))      c2p.c2p_method_backup = 2;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "NR_2DR"))      c2p.c2p_method_backup = 3;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "NR_2DW"))      c2p.c2p_method_backup = 4;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "palenzuela"))  c2p.c2p_method_backup = 5;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "NR_2D_Noble")) c2p.c2p_method_backup = 6;
  else if (CCTK_EQUALS(c2p_algorithm_retry, "Newman"))      c2p.c2p_method_backup = 7;
  else if (CCTK_EQUALS(c2p_preferred_algorithm, "HARM")) c2p.c2p_method = 11;
  if ((c2p.c2p_method >= 12) || (c2p.c2p_method_backup >= 12)){
    printf("Invalid c2p method(s). Terminating...\n");
    exit(1);
  }

  /* T and Ye don't matter if we are dealing with an ideal EOS */
  // ideal gas and Hybrid EOS
  c2p.evolve_T  = 0;
  c2p.evolve_Ye = 0;
  if (c2p.eoskey == 4) {
    c2p.evolve_Ye = 1;
//    c2p.evolve_T  = 1;		/*PE: added this. why wouldn't we evolve T for a realistic?*/
  }

  if(Ye_force_cold_beta_equil){
    c2p.force_cold_beta_equil = 1;
    load_beta_equilibrium_file(beta_equil_file, log_rho0_tab_beta, Y_e_tab, &n_tab_beta);;
  }

  printf("CON2PRIM This is the c2p I am working with:\n");
  printf("CON2PRIM c2p_method: %d\n",c2p.c2p_method);
  printf("CON2PRIM c2p_method_backup: %d\n",c2p.c2p_method_backup);
  printf("CON2PRIM use_c2p_method_backup: %d\n",c2p.use_c2p_method_backup);
  printf("CON2PRIM eoskey: %d\n",c2p.eoskey);
  printf("CON2PRIM eoskey_polytrope: %d\n",c2p.eoskey_polytrope);
  printf("CON2PRIM eos_prec: %g\n",c2p.eos_prec);
  printf("CON2PRIM eos_rho_min: %g\n",c2p.eos_rho_min);
  printf("CON2PRIM eos_rho_max: %g\n",c2p.eos_rho_max);
  printf("CON2PRIM eos_temp_min: %g\n",c2p.eos_temp_min);
  printf("CON2PRIM eos_temp_max: %g\n",c2p.eos_temp_max);
  printf("CON2PRIM eos_eps_min: %g\n",c2p.eos_eps_min);
  printf("CON2PRIM eos_eps_max: %g\n",c2p.eos_eps_max);
  printf("CON2PRIM eos_press_min: %g\n",c2p.eos_press_min);
  printf("CON2PRIM eos_press_max: %g\n",c2p.eos_press_max);
  printf("CON2PRIM evolve_T: %d\n",c2p.evolve_T);
  printf("CON2PRIM evolve_Ye: %d\n",c2p.evolve_Ye);
  printf("CON2PRIM tol_x: %g\n",c2p.tol_x);
  printf("CON2PRIM tol_x_retry: %g\n",c2p.tol_x_retry);
  printf("CON2PRIM max_iterations: %d\n",c2p.max_iterations);
  printf("CON2PRIM extra_iterations: %d\n",c2p.extra_iterations);
  printf("CON2PRIM rho_atmo: %g\n",c2p.rho_atmo);
  printf("CON2PRIM rho_atmo_tol: %g\n",c2p.rho_atmo_tol);
  printf("CON2PRIM T_atmo: %g\n",c2p.T_atmo);
  printf("CON2PRIM Ye_atmo: %g\n",c2p.Ye_atmo);
  printf("CON2PRIM P_atmo: %g\n",c2p.P_atmo);
  printf("CON2PRIM eps_atmo: %g\n",c2p.eps_atmo);
  printf("CON2PRIM retain_B_atmo: %d\n",c2p.retain_B_atmo);
  printf("CON2PRIM numprims: %d\n",c2p.numprims);
  printf("CON2PRIM numcons: %d\n",c2p.numcons);
  printf("CON2PRIM enforce_v2: %d\n",c2p.enforce_v2);
  printf("CON2PRIM force_cold_beta_equil: %d\n",c2p.force_cold_beta_equil);
  printf("CON2PRIM n_tab_beta: %d\n \n",n_tab_beta);

}

struct c2p_report extern_report;
CCTK_INT con2primMHD(CCTK_REAL * prim,
                     const CCTK_REAL * g_up,
                     const CCTK_REAL * g_lo,
                     CCTK_REAL * cons,
                     const CCTK_INT  cac_iter,
                     const CCTK_INT  index,
		     const CCTK_REAL xco,
		     const CCTK_REAL yco, 
		     const CCTK_REAL zco,
		     const CCTK_REAL * METRIC, 
		     const CCTK_REAL * METRIC_PHYS, 
		     const CCTK_REAL * METRIC_LAP_PSI4){
  /* This routine converts the conservatives to the primitives,
     it returns 1 if everything was okay.  */
  /* This routine works on a given point of the grid. */

  DECLARE_CCTK_PARAMETERS

  // con2prim process
  struct c2p_report report;
  bool excise = c2p_excise;
  bool grace = c2p_grace;

  if(c2p.force_cold_beta_equil == 1){
    cons[YE] = prim[YE] * cons[D];
  }
 
  con2prim_MHD_(prim, cons, g_up, g_lo, excise, grace, &report, cac_iter, index,xco,yco,zco,METRIC,METRIC_PHYS,METRIC_LAP_PSI4);

  if((c2p.force_cold_beta_equil==1)&&(prim[RHO]>=1e-7)){
    double yeres = interp_c2p(log_rho0_tab_beta, Y_e_tab, n_tab_beta,log10(prim[RHO]), &n_nearest_beta);
    prim[YE] = yeres;
  }
  extern_report = report;
//  extern_report.failed=report.failed;
//  extern_report.adjust_cons=report.adjust_cons;
//  extern_report.count=err_msg.err_msg;
//  extern_report.count=report.count;
//  extern_report.retry=report.retry;
//  extern_report.c2p_keyerr=report.c2p_keyerr;
//  extern_report.nEOScalls=report.nEOScalls;

  /* for (int i = 0; i < 14; i++) */
  /*   printf("%d %g\n", i, prim[i]); */

  /* if (report.failed) { */
  /*   printf("FAILED: con2prim reconstruction failed for the following reason:\n", report.failed); */
  /*   printf("Report: %s\n", report.err_msg); */
  /* } else { */
  /*   printf("SUCCESS: con2prim converged in %i iterations!\n",report.failed, report.count); */
  /* } */

  if (report.failed) return 1; else return 0;

}
