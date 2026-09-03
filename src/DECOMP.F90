!======================================================================!
subroutine decomp
!----------------------------------------------------------------------!
! This routine is based on the 'EightPoolCenturyMod' code of Manas as
! on 2026-08-31.
!----------------------------------------------------------------------!
use PARS_MOD
use VARS_MOD
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
!T_soil = TC_day
!----------------------------------------------------------------------!
denom = saturation_to_field_capacity * swc_field_capacity
!----------------------------------------------------------------------!
co2 = zero
!----------------------------------------------------------------------!
! swc is volumetric soil water I think. guess calc like this
! really need soil depth. need to make consistent with hydro
!----------------------------------------------------------------------!
do kl = 1, nlayers
  !--------------------------------------------------------------------!
  ! Temperature and sm modifiers from Manas code (EightPoolCenturyMod.F90)
  !--------------------------------------------------------------------!
  !tmod = q10 ** ((T_soil (kl) - T_ref) / 10.0)
  tmod = q10 ** ((T_soil_daily_mean (kl) - T_ref) / 10.0)
  tmod = min (one, tmod)
  tmod = max (zero, tmod)
  !--------------------------------------------------------------------!
  ! Volumetric soil water content (mm/mm)
  !--------------------------------------------------------------------!
  theta (kl) = sm_day (kl) / dz (kl)
  !--------------------------------------------------------------------!
  !wfps = 100.0 * swc (kl) / denom
  ! Eqn. 54 of Friend et a (1997).
  !wfps = 100.0 * theta (kl) / (saturation_to_field_capacity * swc_field_capacity)
  wfps = 100.0 * theta (kl) / theta_sat
  wmod = exp (((wfps - wfps_threshold) ** 2) / (-moisture_dry_width))
  wmod = min (one, wmod)
  wmod = max (zero, wmod)
  !--------------------------------------------------------------------!
  ! Combined temperature and water decay modifier.
  !--------------------------------------------------------------------!
  amod = tmod * wmod
  !--------------------------------------------------------------------!
  ! Active SOM turns over more slowly in fine-textured soils.  The term
  ! 1.0 - 0.75 * (silt + clay) reduces active SOM decay when silt + clay
  ! is high.
  !----------------------------------------------------------------------!
  texture_modifier = one - active_texture_coefficient * &
    (silt_fraction (kl) + clay_fraction (kl))
  texture_modifier = max(zero, texture_modifier)
  !--------------------------------------------------------------------!
  ! Construct effective daily decay rates.
  !--------------------------------------------------------------------!
  k_decay = k_decay_max * amod
  !--------------------------------------------------------------------!
  ! Structural shoot litter: lignin slows decay  using exp (-3 Lignin)
  !--------------------------------------------------------------------!
  k_decay (ip_surface_structural) = k_decay (ip_surface_structural) * &
    exp (-structural_lignin_decay_coefficient * shoot_lignin_frac)
  !--------------------------------------------------------------------!
  ! Structural root litter: root lignin generally produces slower decay.
  !--------------------------------------------------------------------!
  k_decay (ip_soil_structural) = k_decay (ip_soil_structural) * &
    exp (-structural_lignin_decay_coefficient * root_lignin_frac)
  !--------------------------------------------------------------------!
  ! Active SOM: texture stabilization reduces turnover.
  !--------------------------------------------------------------------!
  k_decay (ip_active_som) = k_decay (ip_active_som) * texture_modifier
  !--------------------------------------------------------------------!
  ! Initial state in this layer.
  !--------------------------------------------------------------------!
  c_start (:) = c_state (:,kl)
  !--------------------------------------------------------------------!
  ! Calculate amount decayed from each pool.
  ! Decay = rate * pool_size * time. The min() prevents removing more
  ! carbon than exists in a pool during one step.
  !----------------------------------------------------------------------!
  do ip = 1, n_pools
    decay (ip) = min (k_decay (ip) * c_start (ip) * dt_years, &
                      c_start (ip))
  end do
  !--------------------------------------------------------------------!
  ! Calculate daily plant carbon input and split it into shoots and roots.
  !----------------------------------------------------------------------!
  root_input  = total_input (kl) * root_fraction
  shoot_input = total_input (kl) * (one - root_fraction)
  !------------------------------------------------------------------!
  ! Split shoot and root litter into structural/metabolic pools.
  !--------------------------------------------------------------------!
  ! Fm is the metabolic fraction.  1-Fm is structural fraction.
  ! Shoot and root can have different L:N values.
  !--------------------------------------------------------------------!
  fm_shoot = metabolic_fraction_intercept - &
    metabolic_fraction_lignin_n_slope * shoot_lignin_to_n
  fm_root  = metabolic_fraction_intercept - &
    metabolic_fraction_lignin_n_slope * root_lignin_to_n
  !--------------------------------------------------------------------!
  ! Aboveground plant input goes to surface litter pools.
  !--------------------------------------------------------------------!
  input_vec (ip_surface_structural) = shoot_input * (one - fm_shoot)
  input_vec (ip_surface_metabolic)  = shoot_input * fm_shoot
  !--------------------------------------------------------------------!
  ! Belowground plant input goes to soil/root litter pools.
  !--------------------------------------------------------------------!
  input_vec (ip_soil_structural) = root_input * (one - fm_root)
  input_vec (ip_soil_metabolic)  = root_input * fm_root
  !--------------------------------------------------------------------!
  ! ?
  !--------------------------------------------------------------------!
  transfer_vec = zero
  leached_c = zero
  !--------------------------------------------------------------------!
  ! Pool 1: surface structural litter.
  !--------------------------------------------------------------------!
  co2 = co2 + pool1_surface_structural_co2_fraction * &
        decay (ip_surface_structural)
  transfer_vec (ip_surface_microbe) = &
    transfer_vec (ip_surface_microbe) + &
    pool1_surface_structural_transfer_fraction * &
    decay (ip_surface_structural) * (one - shoot_lignin_frac)
  transfer_vec (ip_slow_som) = transfer_vec (ip_slow_som) + &
    pool1_surface_structural_transfer_fraction * &
    decay (ip_surface_structural) * shoot_lignin_frac
  !--------------------------------------------------------------------!
  ! Pool 2: soil/root structural litter.
  !--------------------------------------------------------------------!
  co2 = co2 + pool2_soil_structural_co2_fraction * &
    decay (ip_soil_structural)
  transfer_vec (ip_active_som) = transfer_vec (ip_active_som) + &
    pool2_soil_structural_transfer_fraction * &
    decay (ip_soil_structural) * (one - root_lignin_frac)
  transfer_vec (ip_slow_som) = transfer_vec (ip_slow_som) + &
    pool2_soil_structural_transfer_fraction * &
    decay (ip_soil_structural) * root_lignin_frac
  !--------------------------------------------------------------------!
  ! Pool 3: active SOM.
  !--------------------------------------------------------------------!
  silt_plus_clay = silt_fraction (kl) + clay_fraction (kl)
  ! CO2 fraction from active SOM decay.  More fine texture lowers ft.
  ft = max (zero, min (active_som_co2_intercept - &
    active_som_co2_silt_clay_slope * silt_plus_clay, one))
  ! Leached C fraction.  More leaching water and sand increase this term.
  cal = max (zero, min ((leaching_cm_day / &
    active_som_leach_water_scale) * (active_som_leach_base + &
    active_som_leach_sand_multiplier * sand_fraction (kl)), one))
  ! Passive SOM formation fraction.  Clay increases passive stabilization.
  cap = max (zero, min (active_som_passive_base_fraction + &
    active_som_passive_clay_multiplier * clay_fraction (kl), one))
  ! Remaining active decay carbon goes to slow SOM.
  cas = max (zero, one - ft - cal - cap)
  ! Safety normalization.  The fractions should normally sum to <= 1.
  ! This guard prevents mass-balance errors if unusual parameter values or
  ! very high leaching make the sum exceed 1.
  total = ft + cal + cas + cap
  if (total > one) then
    ft  = ft  / total
    cal = cal / total
    cas = cas / total
    cap = cap / total
  end if
  co2 = co2 + ft * decay (ip_active_som)
  leached_c = leached_c + cal * decay (ip_active_som)
  transfer_vec (ip_slow_som) = transfer_vec (ip_slow_som) + &
    cas * decay (ip_active_som)
  transfer_vec (ip_passive_som) = transfer_vec (ip_passive_som) + &
    cap * decay (ip_active_som)
  !--------------------------------------------------------------------!
  ! Pool 4 - surface microbial C.
  !----------------------------------------------------------------------!
  co2 = co2 + surface_microbe_co2_fraction * decay (ip_surface_microbe)
  transfer_vec (ip_slow_som) = transfer_vec (ip_slow_som) + &
    surface_microbe_slow_fraction * decay (ip_surface_microbe)
  !--------------------------------------------------------------------!
  ! Pool 5 - surface metabolic litter.
  !--------------------------------------------------------------------!
  co2 = co2 + surface_metabolic_co2_fraction * &
    decay (ip_surface_metabolic)
  transfer_vec (ip_surface_microbe) = &
    transfer_vec (ip_surface_microbe) + &
    surface_metabolic_microbe_fraction * decay (ip_surface_metabolic)
  !--------------------------------------------------------------------!
  ! Pool 6 - soil/root metabolic litter.
  !----------------------------------------------------------------------!
  co2 = co2 + soil_metabolic_co2_fraction * decay (ip_soil_metabolic)
  transfer_vec (ip_active_som) = transfer_vec (ip_active_som) + &
    soil_metabolic_active_fraction * decay(ip_soil_metabolic)
  !--------------------------------------------------------------------!
  ! Pool 7 - slow SOM.
  !----------------------------------------------------------------------!
  csp = max (zero, min (slow_som_passive_base_fraction + &
    slow_som_passive_clay_multiplier * clay_fraction (kl), one))
  co2 = co2 + slow_som_co2_fraction * decay (ip_slow_som)
  transfer_vec (ip_passive_som) = transfer_vec (ip_passive_som) + &
    csp * decay (ip_slow_som)
  transfer_vec (ip_active_som) = transfer_vec (ip_active_som) + &
    (one - slow_som_co2_fraction - csp) * decay (ip_slow_som)
  !--------------------------------------------------------------------!
  ! Pool 8 - passive SOM.
  !--------------------------------------------------------------------!
  co2 = co2 + passive_som_co2_fraction * decay (ip_passive_som)
  transfer_vec (ip_active_som) = transfer_vec (ip_active_som) + &
    passive_som_active_fraction * decay (ip_passive_som)
  !--------------------------------------------------------------------!
  ! Final state.
  !--------------------------------------------------------------------!
  c_end = max (zero, c_start - decay + transfer_vec + input_vec)
  c_state (:,kl) = c_end
  !--------------------------------------------------------------------!
end do ! kl
!--------------------------------------------------------------------!
SOM = sum (c_state)
!--------------------------------------------------------------------!
! Heterotrophic respiration flux (kg[C]/m2/s)
!----------------------------------------------------------------------!
Rhet = co2 / day_s
!----------------------------------------------------------------------!
end subroutine decomp
!======================================================================!
