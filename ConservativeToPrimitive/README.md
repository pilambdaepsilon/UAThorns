# ConservativeToPrimitive

This thorn containes routines which drive the conservative-to-primitive inversion during GR(M)HD evolution. The thorn was originally designed to work with IllinoisGRMHD, but in principle it can be interfaced with any GR(M)HD code that uses the Cactus infrastructure. At its core, ConservativeToPrimitive drives the conservative-to-primitive algorithms implemented in the GRMHD_con2prim code developed by Daniel Siegel and Phillip Moesta. Below we explain the main parameters used to control ConservativeToPrimitive:

* gamma_th: The adiabatic constant used in Hybrid EOSs (such that the EOS is described by a cold component and an appended thermal Gamma-law EOS).
* c2p_EOS_type: type of EOS used.
* c2p_excise: forces the solution to be set to atmosphere.
* c2p_grace: Allow for some additional room around the accepted relative tolerance when setting values to atmosphere.
* c2p_rhoT_key: Set to true if the temperature is prescribed.
* c2p_use_eps_atmo: Use an atmpheric value for the energy density. If false, it recovers the atmosphere based on the atmospheric values of other variables.
* c2p_preferred_algorithm: The primary C2P algorithm used.
* c2p_max_iterations: Maximum number of iterations in C2P algorithm before terminating attempt.
* c2p_extra_iterations: Allow for some extra itreations beyond the maximum, only if algorithm fails to converge on a solution before the maximum number of iterations.
* c2p_tolerance: The relative error threshold for accepting a solution as successful
* c2p_use_backup_scheme: If primary C2P algorithm fails, use a backup scheme.
* c2p_tolerance_retry: control the relative error threshold for the backup scheme.
* c2p_algorithm_retry: Name of the backup C2P scheme.
* c2p_enforce_v2_smaller_than_1: Allow only causal solutions from the recovery.
* c2p_rho_atmo: Atmoshperic value of the rest mass density, used when resetting to atmoshphere after failure. It is best if this is consistent with the value set in the GRMHD code.
* c2p_eps_atmo: Same as above, but for the energy density. Only used if c2p_use_eps_atmo is True.
* c2p_rho_atmo_tolerance: Relative error threshold in setting a point to atmosphere. We reset to atmosphere if RHO <= c2p_rho_atmo * (1 + c2p_rho_atmo_tolerance) (with a little wiggle room if c2p_grace is set to true)
* c2p_T_atmo: atmospheric value for the temperature in MeV
* c2p_ye_atmo: atmospheric value for the electron fraction
* c2p_retain_B_atmo: Leave the value of the magnetic field untouched in the atmosphere.
* Ye_force_cold_beta_equil: Force beta equilibrium solution for the electron fraction. After rho is recovered, assume T=0 and get the corresponding value of Ye in beta equilibrium.
* beta_equil_file: Directory/file containing the beta-equilibrium Ye(rho) function in the following format: 
   \[number of lines in file\]
   \[rest mass density in cgs\] \[Ye\]
   \[...\]                      \[...\]

## Examples
Please see the parameter files provided in the examples directory for examples of how to use ConservativeToPrimitive with IllinoisGRMHD in the case of BNS systems with either a Gamma=2 polytropic EOS or a tabulated finite temperature EOS (LS220)

## LICENSE
ConservativeToPrimitive: a thorn which drives the use of routines from the GRMHD_con2prim code (see below) and allows for the recovery of primitives from conservatives in GRMHD evolutions.

Copyright (C) 2021 by Pedro Espino <pespino@berkeley.edu> and Gabriele Bozzola <gabrielebozzola@email.arizona.edu>

GRMHD_con2prim is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License. You should have received a copy of the license along with this work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.

# GRMHD_con2prim
GRMHD_con2prim is a general-purpose framework that comes with various methods for the recovery of primitive variables from conservative variables in general-relativistic magnetohydrodynamics (GRMHD) with dynamical spacetimes.

GRMHD_con2prim is available at the Zenodo repository with DOI 10.5281/zenodo.1213306.
A separate updated development version is available at https://bitbucket.org/dsiegel/grmhd_con2prim.

The associated detailed methods paper is:

Siegel, D., Moesta, P., Desai, D., Wu, S. (2018), "Recovery schemes for primitive variables in general-relativistic magnetohydrodynamics", ApJS accepted, https://arxiv.org/abs/1712.07538.

The DOI of the Zenodo repository (and the url of the development repository if appropriate) as well as the paper mentioned above should be cited in any academic work that uses GRMHD_con2prim or is based on GRMHD_con2prim, i.e., the code published in this repository.

While this code is free to use (see license below), we ask everyone who uses it to report bugs and to share further development based on this code, for the benefit of everyone else.

## LICENSE

GRMHD_con2prim: a general-purpose framework with various methods for the recovery of primitive variables from conservative variables in general-relativistic magnetohydrodynamics with dynamical spacetimes.

Copyright (C) 2018 by Daniel Siegel <dsiegel@astro.columbia.edu> and Philipp
Moesta <pmoesta@berkeley.edu>

GRMHD_con2prim is licensed under a
Creative Commons Attribution-ShareAlike 4.0 International License.

You should have received a copy of the license along with this
work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.

## DOCUMENTATION

GRMHD_con2prim can be used in GRMHD evolution codes (e.g., Siegel & Metzger (2017), PRL 119, 231102, https://arxiv.org/abs/1705.05473) or in a standalone  fashion to implement new schemes and compare them to existing ones, as in the aforementioned methods paper. Please see USAGE file for more details on how to use this code.

## ACKNOWLEDGEMENTS

Dhruv Desai and Samantha Wu have contributed to this work. Together with other equations of state (EOS), GRMHD_con2prim comes with a modified version of the Helmholtz EOS, originally developed by Timmes & Arnett 1999, ApJS, 125, 277 and Timmes & Swesty 2000, ApJS 126, 501, see http://cococubed.asu.edu/code_pages/eos.shtml, including modifications by Rodrigo Fernandez used in Fernandez & Metzger 2013, MNRAS 435, 502.
