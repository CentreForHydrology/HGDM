  ! A function for calculating the water area of parameterized small ponds
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

  ! The function uses a simple exponential relationship to estimate the
  ! fractional water area from the fractional water volume.

  ! The fractional water area is the fraction of the maximum possible
  ! water area occupied by small depressions in a basin/grid/GRU which 
  ! is occupied by water.
  
  ! Arguments:
  ! vol_frac - state variable. The fraction of the depressional
  ! storage in a basin/grid/GRU which is occupied by water.
  
  ! mult - calculation parameter. The multiplier in the relationship 
  ! between the volume and area fractions. Optional.
  
    
  ! power - calculation parameter. The exponent in the relationship 
  ! between the volume and area fractions. Optional.
  
  ! If mult or power are not specified, the function will use values
  ! estimated from Smith Creek Research Basin. These will NOT be
  ! appropriate for your location.
  

function  small_ponds_water_frac_area (vol_frac, mult, power)
  real :: small_ponds_water_frac_area
  ! arguments
  real, intent(in) :: vol_frac
  real, optional :: mult, power
  
  ! internal variables
  real :: water_frac
  
  ! check for parameter values - if not specified use SCRB values
  if (.not. present(mult)) then
    mult = 1.0058485
  endif
  
  if (.not. present(power)) then
    power = 1.5198
  endif
    
  ! use simple exponential equation to calculate water areal fraction
  if (vol_frac <= 0.0) then
    water_frac = 0.0
  else
    water_frac = mult * vol_frac ** power
  endif
  
  small_ponds_water_frac_area = max(min(water_frac, 1.0), 0.0)
  
end function small_ponds_water_frac_area

