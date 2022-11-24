module dyna_main_loop
    implicit none
contains
    !
    !===========================================================
    !  MainLoop subroutine to set parameters, 
    !  initialise model and run 
    !===========================================================
    !--- Sarah Shannon 2021: Edit to input observed discharge. 
    !--- Use this to calculate NSE.  
    !--- Move to next simulation if NSE < threshold
    subroutine mainloop (all_pm_names,  &
        nac,                            &
        nstep,                          &
        num_rivers,                     &
        acc,                            &
        dt,                             &
        dtt,                            &
        dyna_hru,                       &
        mcpar,                          &
        ntt,                            &
        node_to_flow_mapping,           &
        num_par_types,                  &
        pe_step,                        &
        qobs_riv_step_start,            &
        qobs_riv_step,                  &
        sim_id,                         &
        out_dir_path,                   &        
        r_gau_step,                     &
        rivers,                         &
        route_riv,                      &
        route_tdh,                      &
        sum_ac_riv,                     &
        wt,                             &
        temp_step,                      &
        elev_climate,                   &
        year,                           &
        month,                          &    
        day,                            &
        julian_day)

        use dyna_common_types
        use dyna_param_setup
        use dyna_initialise_run
        use dyna_init_satzone
        use dyna_init_glacier
        use dyna_topmod
        use dyna_write_output
        use dta_route_processing

        implicit none

        ! Argument Declares

        integer :: nac
        integer :: nstep
        integer :: num_rivers
        doubleprecision :: dt
        doubleprecision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
       
        ! Parameters
        character(len=20), dimension(21) :: all_pm_names !names of all possible parameters
        doubleprecision, allocatable, dimension(:,:) :: mcpar
        integer :: num_par_types
        integer :: ntt
        doubleprecision :: wt
        doubleprecision :: acc
        doubleprecision, allocatable, dimension(:) :: srmax
        doubleprecision, allocatable, dimension(:) :: chv
        doubleprecision, allocatable, dimension(:) :: chvdt
        doubleprecision, allocatable, dimension(:) :: lnto
        doubleprecision, allocatable, dimension(:) :: smax
        doubleprecision, allocatable, dimension(:) :: srinit
        doubleprecision, allocatable, dimension(:) :: szm
        doubleprecision, allocatable, dimension(:) :: t0dt
        doubleprecision, allocatable, dimension(:) :: td
        doubleprecision, allocatable, dimension(:) :: min_ddf_snow
        doubleprecision, allocatable, dimension(:) :: max_ddf_snow
        doubleprecision, allocatable, dimension(:) :: tlapse
        doubleprecision, allocatable, dimension(:) :: plapse
        doubleprecision, allocatable, dimension(:) :: rain2snow_temp
        doubleprecision, allocatable, dimension(:) :: melt_temp
        doubleprecision, allocatable, dimension(:) :: tlag_snow
        doubleprecision, allocatable, dimension(:) :: min_ddf_ice
        doubleprecision, allocatable, dimension(:) :: max_ddf_ice
        doubleprecision, allocatable, dimension(:) :: tlag_ice
        doubleprecision, allocatable, dimension(:) :: snow2ice
        doubleprecision, allocatable, dimension(:) :: snow_ice_sublim
        doubleprecision, allocatable, dimension(:) :: rcf
        doubleprecision, allocatable, dimension(:) :: scf
       
    
        ! Input and flow data
        double precision, allocatable, dimension(:,:) :: pe_step
        double precision, allocatable, dimension(:,:) :: r_gau_step
        double precision, allocatable, dimension(:,:) :: temp_step
        double precision, allocatable, dimension(:)   :: elev_climate 
        
        double precision :: qobs_riv_step_start(num_rivers)
        double precision :: qobs_riv_step(num_rivers,nstep)
        integer :: sim_id
        character(900) :: out_dir_path
        integer, dimension(:), allocatable :: julian_day
        integer, dimension(:), allocatable :: year
        integer, dimension(:), allocatable :: month
        integer, dimension(:), allocatable :: day
        
        ! River Data
        doubleprecision :: rivers(nac, num_rivers)
        double precision, dimension(:,:), allocatable :: q
        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh
        integer :: node_to_flow_mapping(:)
        doubleprecision :: sum_ac_riv (num_rivers)

        ! End declares

        !
        !=========================================================
        !  PUT THE PARAM VALUES INTO THE VARIABLE NAMES FOR THE MODEL TYPE
        !=========================================================
        
        call param_setup (chv, &
            chvdt,             &
            dt,                &
            lnto,              &
            num_par_types,     &
            mcpar,             &
            smax,              &
            srinit,            &
            srmax,             &
            szm,               &
            t0dt,              &
            td,                &
            min_ddf_snow,      &
            max_ddf_snow,      &
            tlapse,            &
            plapse,            &
            rain2snow_temp,    &
            melt_temp,         &
            tlag_snow,         &
            min_ddf_ice,       &
            max_ddf_ice,       &
            tlag_ice,          &
            snow2ice,          &
            snow_ice_sublim,   &
            rcf,               &
            scf)
        
        !=========================================================
        !  INITIALISE VARIABLES, ROUTING AND STORES
        !=========================================================
        
        call initialise_run (dyna_hru, &
            dyna_riv,                  &
            nac,                       &
            num_rivers)

        call init_satzone(nac,         &
            nstep,                     &
            num_rivers,                &
            chvdt,                     &
            dt,                        &
            dyna_hru,                  &
            mcpar,                     &
            q,                         &
            qobs_riv_step_start,       &
            rivers,                    &
            route_riv,                 &
            route_tdh,                 &
            smax,                      &
            srinit,                    &
            srmax,                     &
            sum_ac_riv,                &
            szm,                       &
            t0dt)
                
        !=========================================================
        ! Initialise glacier depth and temperature
        !=========================================================    
        call init_glacier(dyna_hru,   &
                          nac)
        
        !=========================================================
        !  CALL THE MAIN TOPMODEL STRUCTURE
        !=========================================================
        !
        call topmod (nac,       &
            nstep,              &
            num_rivers,         &
            acc,                &
            dt,                 &
            dtt,                &
            dyna_hru,           &
            dyna_riv,           &
            node_to_flow_mapping,&
            ntt,                &
            pe_step,            &
            q,                  &
            r_gau_step,         &
            rivers,             &
            route_riv,          &
            route_tdh,          &
            smax,               &
            srmax,              &
            szm,                &
            t0dt,               &
            td,                 &
            wt,                 &
            qobs_riv_step,      &
            sim_id,             &
            out_dir_path,       & 
            temp_step,          &
            elev_climate,       &
            min_ddf_snow,       &
            max_ddf_snow,       &
            tlapse,             &
            plapse,             &
            rain2snow_temp,     &
            melt_temp,          &
            tlag_snow,          &
            year,               &
            month,              &
            day,                &
            julian_day,         &
            min_ddf_ice,        &
            max_ddf_ice,        &
            tlag_ice,           &
            snow2ice,           &
            snow_ice_sublim,    &
            rcf,                &
            scf)               
          
        !
        !=========================================================
        !  WRITE OUTPUTS
        !=========================================================

        call write_output(dyna_hru, &
            mcpar,                  &
            nac,                    &
            all_pm_names,           &
            num_par_types,          &
            num_rivers,             &
            q)                      
            

        deallocate(mcpar)

        return
                                                                        
    end subroutine mainloop

end module dyna_main_loop
