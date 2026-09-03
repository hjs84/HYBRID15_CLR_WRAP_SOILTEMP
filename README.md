# HYBRID15_CLR_WRAP_SOILTEMP

To compile and run TEST_DRIVER.F90 from command line on CSD3:

nice -19 ./compile_HYBRID15_CLR.sh

## File Architecture
* src/HYBRID15_CLR.F90: Main driver module and F2PY entry point (run_model).
* src/INIT.F90: Model initialization subroutine handling site parameters and initial carbon/water pools.
* src/PARS_MOD.F90: Physical constants, array dimensions, and model parameter declarations.
* src/VARS_MOD.F90: Global variable declarations and allocated storage.
* driver.txt: Site-specific initialization data, default pool values.

The run_model subroutine accepts a 21-element real(8) array for parameters (listed in order of highest sensitivity to lowest):

Index | Parameter            | Description                             | Units (I think)
------|----------------------|-----------------------------------------|----------------
1     | wfps_threshold       | Water-filled pore space threshold       | %
2     | T_ref                | Reference temperature for respiration   | deg C
3     | Vcmax_top            | Max carboxylation rate at canopy top    | mol m-2 s-1
4     | theta_sat            | Saturated volumetric soil water content | m3 m-3
5     | q10                  | Temperature sensitivity coefficient     | -
6     | dz(2)                | Soil layer 2 thickness                  | mm
7     | Topt_J               | Temperature optimum for Jmax            | deg C
8     | moisture_dry_width   | Moisture stress function width parameter| -
9     | pool_initial(7,1)    | Initial slow SOM pool (Layer 1)         | gC m-2
10    | pool_initial(7,2)    | Initial slow SOM pool (Layer 2)         | gC m-2
11    | omega_J              | Temperature response curve width        | deg C
12    | Kx                   | Xylem hydraulic conductance             | mol m-2 s-1 MPa-1
13    | lwp_crit             | Critical leaf water potential           | MPa
14    | b_perc               | Percolation exponent                    | -
15    | LAI                  | Leaf area index                         | m2 m-2
16    | pool_initial(3,1)    | Initial active SOM pool (Layer 1)       | gC m-2
17    | pool_initial(3,2)    | Initial active SOM pool (Layer 2)       | gC m-2
18    | asw                  | Available soil water parameter          | -
19    | perc_max             | Maximum percolation rate                | mm d-1
20    | dz(1)                | Soil layer 1 thickness                  | mm
21    | KPAR                 | PAR extinction coefficient              | -

## Repository Variations

* Origin model: https://github.com/hjs84/HYBRID15_CLR.git forked from https://github.com/adfriend45/HYBRID15_CLR.git.
* This Branch: Uses soil temperature model (based on 1D diffusion equation).
* Soil temperature modelled as daily average air temperature is available at:
  https://github.com/hjs84/HYBRID15_CLR_WRAP.git
