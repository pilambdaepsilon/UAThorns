Authors: Gabriele Bozzola and Pedro Espino

# ConservativeToPrimitive
This thorn drives the conservatives-to-primitives conversion using the algorithms implemented in GRMHD_con2prim,
which was originally written by Daniel Siegel and Phillip Moesta. This thorn contains the aforementioned implementations of C2P
algorithms, with modifications required to interface it with the EOS driver used within the EinsteinToolkit (EOS_Omni) and 
IllinoisGRMHD. In principle, the thorn may be used with any GR(M)HD but one must be careful in using consistent definitions
of the ingoing conservatives and outgoing primitives.

# LICENSE

ConservativeToPrimitive: a driver of the several conservative-to-primitive routines and algorithms implemented in the GRMHD_con2prim code.
Copyright (C) 2021 by Gabriele Bozzola <gabrielebozzola@email.arizona.edu> and Pedro Espino <pespino@berkeley.edu>
ConservativeToPrimitive is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.
You should have received a copy of the license along with this work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.

#=============================================================================================================================
#=============================================================================================================================

# GRMHD_con2prim
GRMHD_con2prim is a general-purpose framework that comes with various methods for the recovery of primitive variables from 
conservative variables in general-relativistic magnetohydrodynamics (GRMHD) with dynamical spacetimes. GRMHD_con2prim is available 
at the Zenodo repository with DOI 10.5281/zenodo.1213306.
A separate updated development version is available at https://bitbucket.org/dsiegel/grmhd_con2prim. The associated detailed methods paper is:

Siegel, D., Moesta, P., Desai, D., Wu, S. (2018), "Recovery schemes for primitive variables in general-relativistic magnetohydrodynamics", ApJS accepted, https://arxiv.org/abs/1712.07538.

The DOI of the Zenodo repository (and the url of the development repository if appropriate) as well as the paper mentioned above should be cited 
in any academic work that uses GRMHD_con2prim or is based on GRMHD_con2prim, i.e., the code published in this repository.

While this code is free to use (see license below), we ask everyone who uses it to report bugs and to share further development based on this code, for the benefit of everyone else.

# LICENSE

GRMHD_con2prim: a general-purpose framework with various methods for the recovery of primitive variables from conservative variables in general-relativistic magnetohydrodynamics with dynamical spacetimes.
Copyright (C) 2018 by Daniel Siegel <dsiegel@astro.columbia.edu> and Philipp
Moesta <pmoesta@berkeley.edu>
GRMHD_con2prim is licensed under a
Creative Commons Attribution-ShareAlike 4.0 International License.
You should have received a copy of the license along with this
work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.

# DOCUMENTATION

GRMHD_con2prim can be used in GRMHD evolution codes (e.g., Siegel & Metzger (2017), PRL 119, 231102, https://arxiv.org/abs/1705.05473) or in a standalone  fashion to implement new schemes and compare them to existing ones, as in the aforementioned methods paper. Please see USAGE file for more details on how to use this code.

# ACKNOWLEDGEMENTS

Dhruv Desai and Samantha Wu have contributed to this work. Together with other equations of state (EOS), GRMHD_con2prim comes with a modified version of the Helmholtz EOS, originally developed by Timmes & Arnett 1999, ApJS, 125, 277 and Timmes & Swesty 2000, ApJS 126, 501, see http://cococubed.asu.edu/code_pages/eos.shtml, including modifications by Rodrigo Fernandez used in Fernandez & Metzger 2013, MNRAS 435, 502.
