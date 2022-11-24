module dyna_param_setup
contains
    !
    !===============================================================
    !  ROUTINES FOR SETTING UP THE PARAMETER VALUES FOR THE DIFFERENT
    !  HRUs
    !  Sarah Shannon (2021): Monte Carlo friendly sampling of 
    !  snow and glacier params.  
    !===============================================================
    !
    subroutine param_setup (chv, &
         chvdt,                  &
         dt,                     &
         lnto,                   &
         num_par_types,          &
         mcpar,                  &
         smax,                   &
         srinit,                 &
         srmax,                  &
         szm,                    &
         t0dt,                   &
         td,                     &
         min_ddf_snow,           &
         max_ddf_snow,           &
         tlapse,                 &
         plapse,                 &
         rain2snow_temp,         &
         melt_temp,              &
         tlag_snow,              &
         min_ddf_ice,            &
         max_ddf_ice,            &
         tlag_ice,               &
         snow2ice,               &
         snow_ice_sublim,        &
         rcf,                    &
         scf)
         

        use dyna_common_types

        implicit none

        ! Argument Declares

        double precision, dimension(:,:) :: mcpar

        double precision, allocatable, dimension(:) :: chv
        double precision, allocatable, dimension(:) :: chvdt
        double precision, allocatable, dimension(:) :: lnto
        double precision, allocatable, dimension(:) :: smax
        double precision, allocatable, dimension(:) :: srinit
        double precision, allocatable, dimension(:) :: srmax
        double precision, allocatable, dimension(:) :: szm
        double precision, allocatable, dimension(:) :: t0dt
        double precision, allocatable, dimension(:) :: td
        double precision, allocatable, dimension(:) :: min_ddf_snow   !--- m w e oC-1 day-1 
        double precision, allocatable, dimension(:) :: max_ddf_snow   !--- m w e oC-1 day-1 
        double precision, allocatable, dimension(:) :: tlapse         !--- degC per km 
        double precision, allocatable, dimension(:) :: plapse         !--- % per 100 m  
        double precision, allocatable, dimension(:) :: rain2snow_temp !--- degC
        double precision, allocatable, dimension(:) :: melt_temp      !--- degC [set to zero]
        double precision, allocatable, dimension(:) :: tlag_snow      !--- dimensionless
        double precision, allocatable, dimension(:) :: min_ddf_ice    !--- m w e oC-1 day-1  
        double precision, allocatable, dimension(:) :: max_ddf_ice    !--- m w e oC-1 day-1   
        double precision, allocatable, dimension(:) :: tlag_ice       !--- dimensionless
        double precision, allocatable, dimension(:) :: snow2ice       !--- year [1/(nyear*365) nyear=100]
        double precision, allocatable, dimension(:) :: snow_ice_sublim!--- dimensionless
        double precision, allocatable, dimension(:) :: rcf            !--- dimensionless
        double precision, allocatable, dimension(:) :: scf            !--- dimensionless
        double precision :: dt
        integer :: num_par_types

        ! Local Declares
        integer :: i

        ! End declares

        allocate(chv(num_par_types))
        allocate(chvdt(num_par_types))
        allocate(lnto(num_par_types))
        allocate(smax(num_par_types))
        allocate(srinit(num_par_types))
        allocate(srmax(num_par_types))
        allocate(szm(num_par_types))
        allocate(t0dt(num_par_types))
        allocate(td(num_par_types))
        allocate(min_ddf_snow(num_par_types))
        allocate(max_ddf_snow(num_par_types))
        allocate(tlapse(num_par_types))
        allocate(plapse(num_par_types))
        allocate(rain2snow_temp(num_par_types))
        allocate(melt_temp(num_par_types))
        allocate(tlag_snow(num_par_types))
        allocate(min_ddf_ice(num_par_types))
        allocate(max_ddf_ice(num_par_types))
        allocate(tlag_ice(num_par_types))
        allocate(snow2ice(num_par_types))
        allocate(snow_ice_sublim(num_par_types))
        allocate(rcf(num_par_types))
        allocate(scf(num_par_types))

        do i = 1, num_par_types

           !--- default params
           szm(i)               = mcpar(i,1)
           lnto(i)              = mcpar(i,2)
           srmax(i)             = mcpar(i,3)
           srinit(i)            = mcpar(i,4)
           chv(i)               = mcpar(i,5)
           td(i)                = mcpar(i,6)
           smax(i)              = mcpar(i,7)
           chvdt(i)             = CHV (i) * dt
           t0dt(i)              = LnTo (i) + Log (dt)
           
           !--- snow params. Use scale factor to find min ddf of snow    
           max_ddf_snow(i)      = mcpar(i,8)                   
           min_ddf_snow(i)      = max_ddf_snow(i) * mcpar(i,9) 
           tlapse(i)            = mcpar(i,10)                  
           plapse(i)            = mcpar(i,11)
           rain2snow_temp(i)    = mcpar(i,12)
           melt_temp(i)         = mcpar(i,13)
           tlag_snow(i)         = mcpar(i,14)
           
           !--- glacier params. Use scale factor to find min and max ddf of ice and tlag_ice
           tlag_ice(i)          = tlag_snow(i) * mcpar(i,15)
           min_ddf_ice(i)       = min_ddf_snow(i) * mcpar(i,16)
           max_ddf_ice(i)       = max_ddf_snow(i) * mcpar(i,16)
           snow2ice(i)          = mcpar(i,17) 
           snow_ice_sublim(i)   = mcpar(i,18) 
           rcf(i)               = mcpar(i,19)
           scf(i)               = mcpar(i,20)
             
!--- To do error check and force crash if any of these are true
!mcpar(i,16) < 1, mcpar(i,9) > 1 and mcpar(i,15) > 1

                
        end do

        end subroutine

    end module dyna_param_setup
