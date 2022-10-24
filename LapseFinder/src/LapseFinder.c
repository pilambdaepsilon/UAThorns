#include <math.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loopcontrol.h"

void LapseFinder_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *LapseFinder_alp_min = 99999.;
}

/* Reduce grid functions */
void LapseFinder_Reduction(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT handle_min = CCTK_ReductionHandle("minimum");
//  CCTK_INT handle_min = CCTK_ReductionArrayHandle("minimum");
  if (handle_min < 0)
    CCTK_WARN(0, "Unable to get reduction handle 'min'.");
  if (LapseFinder_comp_alp_min_every != 0 &&
      cctk_iteration % LapseFinder_comp_alp_min_every == 0)
  {
    if (CCTK_Reduce(cctkGH, -1, handle_min, 1, CCTK_VARIABLE_REAL,
                    LapseFinder_alp_min, 1,
                    CCTK_VarIndex("ADMBase::alp")))
      CCTK_WARN(0, "Error while reducing ADMBase::alp");
  }
}

/* Look for the location of the alpha minimum */
static CCTK_REAL local_alp_min_loc[4]; /* these two are shared between the next three routines */
static int have_warned_about_multiple_minima = 0;
void LapseFinder_LocationSearch_Setup(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(LapseFinder_comp_alp_min_every == 0 ||
     cctk_iteration % LapseFinder_comp_alp_min_every != 0)
    return;

  local_alp_min_loc[0] = 0.0;
  local_alp_min_loc[1] = 0.0;
  local_alp_min_loc[2] = 0.0;
  local_alp_min_loc[3] = 0.0; // Information if it was found at all
}

void LapseFinder_LocationSearch_Search(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(LapseFinder_comp_alp_min_every == 0 ||
     cctk_iteration % LapseFinder_comp_alp_min_every != 0)
    return;

  /* Look for the location of the global minimum.
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
      if ( (alp[i3D] <= *LapseFinder_alp_min*(1 + LapseFinder_tolerance)) &&
           (alp[i3D] >= *LapseFinder_alp_min*(1 - LapseFinder_tolerance)) &&
          (!LapseFinder_alp_min_loc_only_positive_x || x[i3D] >= 0.0))
      {
        #pragma omp critical
        {
          // have to already warn here when we still have the actual coordinates around
          if (verbosity_level >= 1 && local_alp_min_loc[3] >= 1.)
          {
            if(!have_warned_about_multiple_minima && local_alp_min_loc[3] == 1.) // *second* time we get here
              CCTK_WARN(1, "Found more than one identical minimum on single processor.");
            if (verbosity_level >= 2)
            {
              if (round(local_alp_min_loc[3]) == 1.) { // once we detect the second minimum, output the first as well
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "my candidate: (%g,%g,%g) with value %g.",
                           local_alp_min_loc[0],local_alp_min_loc[1],local_alp_min_loc[2],
                           *LapseFinder_alp_min);
              }
              CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "my candidate: (%g,%g,%g) with value %g.",
                         x[i3D], y[i3D], z[i3D],
                         *LapseFinder_alp_min);
            }
          }
            local_alp_min_loc[0] += x[i3D];
            local_alp_min_loc[1] += y[i3D];
            local_alp_min_loc[2] += z[i3D];
            local_alp_min_loc[3] += 1.;
          }
        }
      }
    LC_ENDLOOP3(location_reduction);
    }
  }
void LapseFinder_LocationSearch_Combine(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int set_minimum_location;

  CCTK_REAL global_alp_min_loc[4];

  if(LapseFinder_comp_alp_min_every == 0 ||
     cctk_iteration % LapseFinder_comp_alp_min_every != 0)
    return;

  /* At this point at least one MPI process should have the location set.
   * The others should have a value of 0.0, so that the sum reduction
   * gives the value on each MPI process. This however is not true sometimes
   * and this is the point where the algorithm can fail. */
  //CCTK_INT handle_sum = CCTK_LocalArrayReductionHandle("sum");
  CCTK_INT handle_sum = CCTK_ReductionArrayHandle("sum");
  if (handle_sum < 0)
     CCTK_WARN(0, "Unable to get reduction handle 'sum'.");
  if (CCTK_ReduceLocArrayToArray1D(cctkGH, -1, handle_sum, &local_alp_min_loc,
                                   global_alp_min_loc, 4, CCTK_VARIABLE_REAL))
    CCTK_WARN(0, "Error while reducing local_alp_min_loc");
  if (round(global_alp_min_loc[3]) > 0.) {
    if (round(global_alp_min_loc[3]) == 1.)
    {
      set_minimum_location = 1;
    } else {
      if (LapseFinder_average_multiple_minima_locations)
      {
        global_alp_min_loc[0] /= global_alp_min_loc[3];
        global_alp_min_loc[1] /= global_alp_min_loc[3];
        global_alp_min_loc[2] /= global_alp_min_loc[3];

        if (verbosity_level >= 1 && !have_warned_about_multiple_minima)
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Found more than one (%d) identical minimum, using average location (%g,%g,%g).",
                     (int)lrint(global_alp_min_loc[3]),
                     global_alp_min_loc[0], global_alp_min_loc[1],
                     global_alp_min_loc[2]);
          have_warned_about_multiple_minima = 1;
        }

        set_minimum_location = 1;
      } else {
        if (verbosity_level >= 1)
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Found more than one (%d) identical minimum, not setting anything.",
                     (int)lrint(global_alp_min_loc[3]));
          if (verbosity_level >= 2 && local_alp_min_loc[3] == 1.0) // in this case we did not already warn about this locally
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "my candidate: (%g,%g,%g) with value %g.",
                       local_alp_min_loc[0],local_alp_min_loc[1],local_alp_min_loc[2],
                       *LapseFinder_alp_min);
          }
        }
        set_minimum_location = 0;
      }
    }

  } else {
    if (verbosity_level >= 1) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not find grid point with minimum density %g on grid. Using old location  (%g,%g,%g).",
                 *LapseFinder_alp_min,
                 LapseFinder_alp_min_loc[0], LapseFinder_alp_min_loc[1],
                 LapseFinder_alp_min_loc[2]);
    }
    set_minimum_location = 0;
  }


  if (set_minimum_location)
  {
      LapseFinder_alp_min_loc[0] = global_alp_min_loc[0];
      LapseFinder_alp_min_loc[1] = global_alp_min_loc[1];
      LapseFinder_alp_min_loc[2] = global_alp_min_loc[2];
  }
/*  CCTK_VInfo(CCTK_THORNSTRING, "New location: %g,%g,%g",
    LapseFinder_alp_min_loc[0],
    LapseFinder_alp_min_loc[1],
    LapseFinder_alp_min_loc[2]); */
}
