#include <math.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loopcontrol.h"

void RhoBFinder_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *RhoBFinder_rho_b_max = 99999.;
}

/* Reduce grid functions */
void RhoBFinder_Reduction(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (RhoBFinder_compute_only_first && cctk_iteration != 0)
    return;

  CCTK_INT handle_max = CCTK_ReductionHandle("maximum");
//  CCTK_INT handle_max = CCTK_ReductionArrayHandle("maximum");
  if (handle_max < 0)
    CCTK_WARN(0, "Unable to get reduction handle 'max'.");
  if (RhoBFinder_comp_rho_b_max_every != 0 &&
      cctk_iteration % RhoBFinder_comp_rho_b_max_every == 0)
  {
    if (CCTK_Reduce(cctkGH, -1, handle_max, 1, CCTK_VARIABLE_REAL,
                    RhoBFinder_rho_b_max, 1,
                    CCTK_VarIndex("IllinoisGRMHD::rho_b")))
      CCTK_WARN(0, "Error while reducing IllinoisGRMHD::rho_b");
  }
}

/* Look for the location of the rho_b maximum */
static CCTK_REAL local_rho_b_max_loc[4]; /* these two are shared between the next three routines */
static int have_warned_about_multiple_maxima = 0;
void RhoBFinder_LocationSearch_Setup(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (RhoBFinder_compute_only_first && cctk_iteration != 0)
    return;

  if(RhoBFinder_comp_rho_b_max_every == 0 ||
     cctk_iteration % RhoBFinder_comp_rho_b_max_every != 0)
    return;

  local_rho_b_max_loc[0] = 0.0;
  local_rho_b_max_loc[1] = 0.0;
  local_rho_b_max_loc[2] = 0.0;
  local_rho_b_max_loc[3] = 0.0; // Information if it was found at all
}

void RhoBFinder_LocationSearch_Search(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(RhoBFinder_comp_rho_b_max_every == 0 ||
     cctk_iteration % RhoBFinder_comp_rho_b_max_every != 0)
    return;

  /* Look for the location of the global maximum.
   * This algorithm will have problems when that occurs at more than
   * one location. */
  #pragma omp parallel
  {
    LC_LOOP3(location_reduction, i,j,k,
             cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2],
             cctk_lsh[0] - cctk_nghostzones[0], cctk_lsh[1] - cctk_nghostzones[1],
             cctk_lsh[2] - cctk_nghostzones[2], cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
      if (rho_b[i3D] == *RhoBFinder_rho_b_max &&
          (!RhoBFinder_rho_b_max_loc_only_positive_x || x[i3D] >= 0.0))
      {
        #pragma omp critical
        {
          // have to already warn here when we still have the actual coordinates around
          if (verbosity_level >= 1 && local_rho_b_max_loc[3] >= 1.)
          {
            if(!have_warned_about_multiple_maxima && local_rho_b_max_loc[3] == 1.) // *second* time we get here
              CCTK_WARN(1, "Found more than one identical maximum on single processor.");
            if (verbosity_level >= 2)
            {
              if (round(local_rho_b_max_loc[3]) == 1.) { // once we detect the second maximum, output the first as well
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "my candidate: (%g,%g,%g) with value %g.",
                           local_rho_b_max_loc[0],local_rho_b_max_loc[1],local_rho_b_max_loc[2],
                           *RhoBFinder_rho_b_max);
              }
              CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "my candidate: (%g,%g,%g) with value %g.",
                         x[i3D], y[i3D], z[i3D],
                         *RhoBFinder_rho_b_max);
            }
          }
            local_rho_b_max_loc[0] += x[i3D];
            local_rho_b_max_loc[1] += y[i3D];
            local_rho_b_max_loc[2] += z[i3D];
            local_rho_b_max_loc[3] += 1.;
          }
        }
      }
    LC_ENDLOOP3(location_reduction);
    }
  }
void RhoBFinder_LocationSearch_Combine(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int set_maximum_location;

  CCTK_REAL global_rho_b_max_loc[4];

  if (RhoBFinder_compute_only_first && cctk_iteration != 0)
    return;


  if(RhoBFinder_comp_rho_b_max_every == 0 ||
     cctk_iteration % RhoBFinder_comp_rho_b_max_every != 0)
    return;

  /* At this point at least one MPI process should have the location set.
   * The others should have a value of 0.0, so that the sum reduction
   * gives the value on each MPI process. This however is not true sometimes
   * and this is the point where the algorithm can fail. */
//  CCTK_INT handle_sum = CCTK_LocalArrayReductionHandle("sum");
  CCTK_INT handle_sum = CCTK_ReductionArrayHandle("sum");
  if (handle_sum < 0)
     CCTK_WARN(0, "Unable to get reduction handle 'sum'.");
  if (CCTK_ReduceLocArrayToArray1D(cctkGH, -1, handle_sum, &local_rho_b_max_loc,
                                   global_rho_b_max_loc, 4, CCTK_VARIABLE_REAL))
    CCTK_WARN(0, "Error while reducing local_rho_b_max_loc");
  if (round(global_rho_b_max_loc[3]) > 0.) {
    if (round(global_rho_b_max_loc[3]) == 1.)
    {
      set_maximum_location = 1;
      CCTK_VInfo(CCTK_THORNSTRING, "Found maximum at (%f, %f, %f) with value %4.3e",
                global_rho_b_max_loc[0], global_rho_b_max_loc[1], global_rho_b_max_loc[2],
                *RhoBFinder_rho_b_max);
    } else {
      if (RhoBFinder_average_multiple_maxima_locations)
      {
        global_rho_b_max_loc[0] /= global_rho_b_max_loc[3];
        global_rho_b_max_loc[1] /= global_rho_b_max_loc[3];
        global_rho_b_max_loc[2] /= global_rho_b_max_loc[3];

        if (verbosity_level >= 1 && !have_warned_about_multiple_maxima)
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Found more than one (%d) identical maximum, using average location (%g,%g,%g).",
                     (int)lrint(global_rho_b_max_loc[3]),
                     global_rho_b_max_loc[0], global_rho_b_max_loc[1],
                     global_rho_b_max_loc[2]);
          have_warned_about_multiple_maxima = 1;
        }

        set_maximum_location = 1;
      } else {
        if (verbosity_level >= 1)
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Found more than one (%d) identical maximum, not setting anything.",
                     (int)lrint(global_rho_b_max_loc[3]));
          if (verbosity_level >= 2 && local_rho_b_max_loc[3] == 1.0) // in this case we did not already warn about this locally
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "my candidate: (%g,%g,%g) with value %g.",
                       local_rho_b_max_loc[0],local_rho_b_max_loc[1],local_rho_b_max_loc[2],
                       *RhoBFinder_rho_b_max);
          }
        }
        set_maximum_location = 0;
      }
    }

  } else {
    if (verbosity_level >= 1) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not find grid point with maximum density %g on grid. Using old location  (%g,%g,%g).",
                 *RhoBFinder_rho_b_max,
                 RhoBFinder_rho_b_max_loc[0], RhoBFinder_rho_b_max_loc[1],
                 RhoBFinder_rho_b_max_loc[2]);
    }
    set_maximum_location = 0;
  }


  if (set_maximum_location)
  {
      RhoBFinder_rho_b_max_loc[0] = global_rho_b_max_loc[0];
      RhoBFinder_rho_b_max_loc[1] = global_rho_b_max_loc[1];
      RhoBFinder_rho_b_max_loc[2] = global_rho_b_max_loc[2];
  }
/*  CCTK_VInfo(CCTK_THORNSTRING, "New location: %g,%g,%g",
    RhoBFinder_rho_b_max_loc[0],
    RhoBFinder_rho_b_max_loc[1],
    RhoBFinder_rho_b_max_loc[2]); */
}
