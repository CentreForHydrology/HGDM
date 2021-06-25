  ! Converted from R -> Fortran by Kevin Shook, June 2020
  ! 
  ! This function applies net precip to water, and calculates runoff depth
  ! Because the areas of application (of direct precip., evap, and runoff)
  ! depend on the storage
  
  ! The function arguments are mostly dimensionless fractions or 
  ! areal-average depths, so that that they are more easily shared among
  ! basins. However, they need to be dimensionalised within the function.
  
  ! Arguments:
  ! current_depth - state variable. Current depth of water (m) ) in 
  ! depressional storage
  
  ! delta_depth - flux. Change of depth of water to be applied (m) 
  ! directly to the small pond water area (i.e. direct precip or evap.), 
  ! either positive or negative
  
  ! runoff_depth - flux. Depth of upland runoff running into 
  ! small ponds (m)

  ! max_depth - state parameter. Maximum depth of water (m) that can be
  ! stored in the small ponds
  
  ! max_water_area_frac - state parameter. Max possible areal fraction of
  ! the small ponds basin/GRU that can be covered with water.
  
  ! application_iterations - subroutine execution parameter. Number of 
  ! iterations used for applying fluxes 
  
  ! area_mult, area_power - scaling parameters. Multiplier and exponent
  ! relating small pond water area to water depth.
  
  ! vol_frac - state variable. Fraction of maximum storage volume filled
  ! by water. Updated by subroutine.
  
  ! area_frac - state variable. Fraction of maximum water area covered
  ! by water. Updated by subroutine.
  
  ! contrib_frac - state variable. Fracation of the GRU/basin affected
  ! by small ponds, which is contributing runoff. Updated by subroutine.
  
  ! depth - state variable. Current depth of water (m) on the GRU/basin.
  ! Updated by subroutine.
  
  ! outflow_depth - flux. Depth of outflow (m) from small ponds. Output 
  ! by subroutine.

subroutine small_depressions_delta_water(current_depth, delta_depth, runoff_depth, contrib_frac,
                            max_depth, max_water_area_frac, application_iterations,
                            area_mult, area_power, vol_frac, area_frac,
                            contrib_frac, depth, outflow_depth) 

! arguments
  real :: current_depth, delta_depth, contrib_frac, max_depth
  real ::  max_water_area_frac, area_mult, area_power
  integer :: application_iterations

! internal variables
  real :: current_vol_frac, iteration_delta, iteration_runoff, iteration_contrib_frac
  real :: iteration_area_frac, total_outflow_depth, water_area_frac, upland_area_frac
  real :: delta_depth_applied, current_depth

  current_vol_frac = current_depth / max_depth  ! current fractional storage volume
  iteration_delta = delta_depth / application_iterations
  iteration_runoff = runoff_depth / application_iterations
  iteration_contrib_frac = contrib_frac
  iteration_area_frac = water_frac_area(current_vol_frac, area_mult, area_power)
  total_outflow_depth = 0
  
  ! check to see if worth doing
  if ((iteration_delta /= 0) .or. (iteration_runoff > 0)) then
    ! apply fluxes iteratively, changing the areas of application
    do i = 1, application_iterations
      water_area_frac = max_water_area_frac * iteration_area_frac   ! basin/HRU areal fraction of water
      upland_area_frac = 1 - water_area_frac
      delta_depth_applied = (iteration_delta * water_area_frac)  +
        iteration_runoff * upland_area_frac * (1.0 - iteration_contrib_frac)

      current_depth = current_depth + delta_depth_applied
      
      ! check to see if max depth of storage is exceeded
      if (current_depth > max_depth) then
      ! max depth is exceeded, so there will be outflows
        excess_depth = max(current_depth - max_depth, 0)
        current_depth = max_depth
        iteration_contrib_frac = 1
        total_outflow_depth = total_outflow_depth +
          (excess_depth * water_area_frac) +
          (iteration_runoff * upland_area_frac * iteration_contrib_frac)

      else 
      ! max depth is not exceeded, so no outflows
        excess_depth = 0
        total_outflow_depth = total_outflow_depth +
          (iteration_runoff * upland_area_frac * iteration_contrib_frac)
        ! update contributing fraction for next iteration
        iteration_contrib_frac = PHM_contrib_frac(iteration_contrib_frac, current_depth,
                                                   delta_depth_applied, max_depth)
      endif

      ! update state variables
      current_vol_frac = current_depth / max_depth
      iteration_area_frac = PHM_water_frac_area(current_vol_frac, area_mult, area_power)
    enddo

    vol_frac = current_depth  / max_depth
    area_frac = iteration_area_frac
    contrib_frac = iteration_contrib_frac
    depth = current_depth
    outflow_depth = total_outflow_depth

  else 
    ! nothing to do - return initial values
    vol_frac = current_depth / max_depth
    area_frac = iteration_area_frac
    contrib_frac = iteration_contrib_frac
    depth = current_depth
    outflow_depth = 0
  endif

end subroutine
