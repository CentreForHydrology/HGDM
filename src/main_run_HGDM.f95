! A fortran95 program for the Hysteretic and Gatekeeping Depressions Model (HGDM) model
! By Mohamed Ahmed
! This file reads the hypothetical inputs and calls the HGDM algorithm
program main


    use HGDM

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


    ! arguments
    real :: current_depth, delta_depth, contrib_frac, max_depth
    real ::  max_water_area_frac, area_mult, area_power
    !integer :: application_iterations

    ! internal variables
    !for main only
    character(len=100) :: trash
    real :: in_data(20,2)
    integer :: i

    !outputs
    real :: depth, outflow_depth


    ! define initial conditions
    current_depth = 0.2 !m

    !define model parameters
    max_depth = 1.25 !m
    max_water_area_frac = 0.35 !fraction
    area_mult = 1.0058485 !for SCRB
    area_power = 1.5198 !for SCRB


    !initialize variables
    contrib_frac = 0.0
    vol_frac = 0.0
    area_frac = 0.0
    depth = 0.0
    outflow_depth = 0.0


    !read fluxes (inputs) from file
    !!read input file

    open(1,file='input_depth.txt', status='old')
    read(1,*)trash
    do i=1,20
        read(1,*) in_data(i,:)
    end do
    close(1)

    !write(*,*) in_data


    !open file for output
    open(1,file='output.txt', status='unknown')
    !loop each time step
    do i=1,20
        if(i .eq. 20)then
            k=1
        end if
        !input fluxes
        delta_depth=in_data(i,1) ! in m
        runoff_depth=in_data(i,2) !in m
        call small_depressions_delta_water(current_depth, delta_depth, runoff_depth, contrib_frac, &
                                max_depth, max_water_area_frac, & !application_iterations,
                                area_mult, area_power, vol_frac, area_frac, & !vol and area frac state variables output
                                depth, outflow_depth) !output


        if(i .eq. 1) write(1,*)'i ','depth ','vol_frac ','area_frac ','contrib_frac ','outflow_depth '

        write(1,*)i,depth,vol_frac,area_frac,contrib_frac,outflow_depth

    end do

    close(1)

stop
end program main
