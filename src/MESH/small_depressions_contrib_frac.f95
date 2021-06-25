  ! A function for calculating the contributing fraction of 
  ! parameterized small depressions
  ! Copyright (C) 2020  Kevin Shook

  !  This program is free software; you can redistribute it and/or modify
  !  it under the terms of the GNU General Public License as published by
  !  the Free Software Foundation; either version 3 of the License, or
  !  (at your option) any later version.
  !
  !  This program is distributed in the hope that it will be useful,
  !  but WITHOUT ANY WARRANTY; without even the implied warranty of
  !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  !  GNU General Public License for more details.
  !
  !  You should have received a copy of the GNU General Public License
  !  along with this program; if not, write to the Free Software
  !  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
  !  MA 02110-1301, USA. 

  ! Converted from R -> Fortran by Kevin Shook, June 2020

  ! This function uses a simple linear form of the
  ! hysteretic relationship between the volumetric water fraction
  ! and the contributing fraction of a basin/grid/GRU.
  
  ! In this function, the contributing fraction refers to the
  ! areal fraction of a basin/grid/GRU which can contribute
  ! overland flow to its outlet, whether or not flow is actually
  ! occurring.
  
  ! The contributing fraction is estimated inside the space defined
  ! by the ranges of the volume fraction (x) and 
  ! the contributing fraction (y), i.e. between (0,0) and (1,1).
  
  ! If there is evaporation, then the contributing fraction is zero.
  
  ! If there is an addition of water, then the contributing fraction
  ! is estimated from the current value, on a trajectory to (1,1).
  
  ! Note that because the contributing fraction of the small depressions
  ! varies continuously with the application of water, this function
  ! needs to be called repeatedly, as water is added to the depressions.
  ! This is shown in the subroutine small_ponds_delta_water.
  
  ! Arguments:
  ! current_contribfrac - state variable. Current version of the
  ! contributing fraction of small depressions in the basin/grid/GRU
  
  ! current_depth - state variable. Depth of water (m) in the fraction
  ! of the basin/grid/GRU occupied by small depressions.
  
  ! delta_depth - flux. Depth of water (m) applied to the fraction of the
  ! basin/grid/GRU occupied by small depressions.
  
  ! max-depth - limiting value. Maximum possible depth of water (m) in 
  ! the fraction of the basin/grid/GRU occupied by small depressions.
  
  

function small_depressions_contrib_frac(current_contrib_frac, &
                                  current_depth, delta_depth, max_depth)
  real :: small_depressions_contrib_frac
  
  ! arguments
  real, intent(in) :: current_contrib_frac, current_depth, delta_depth, max_depth
  
  ! internal variables
  real :: new_contribfrac, volfrac, volfrac_delta, slope

  if (delta_depth == 0.0 ) then
    new_contribfrac = current_contrib_frac
  else 
    if (delta_depth < 0.0) then
      ! evaporation
      new_contribfrac = 0.0
    else 
      ! direct precip and/or runoff
      volfrac = current_depth / max_depth  ! current fractional water vol
      volfrac_delta = delta_depth / max_depth ! change in water vol frac
      
      ! calculate the slope of the line connecting the current point
      ! to (1, 1)
      slope = (1.0 - current_contrib_frac) / (1.0 - volfrac)
      slope = max(slope, 0.0)  ! slope can't be < 0 for addition of water
      
      ! using the slope, calculate the new contributing fraction
      new_contribfrac = current_contrib_frac + (volfrac_delta * slope)
    endif
  endif
  
  small_depressions_contrib_frac = min(max(new_contribfrac, 0.0), 1.0)

end function small_depressions_contrib_frac
