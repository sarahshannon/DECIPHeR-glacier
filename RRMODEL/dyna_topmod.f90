module dyna_topmod
contains
    !
    !=================================================================
    !  DynaTOP 1-2 (more modular structure for development)
    !=================================================================
    !  Keith Beven, Lausanne, 1997: Jim Freer Lancaster 12/1/00
    !  Updated by Gemma Coxon and Toby Dunne - update to Fortran 2003, made more modular
    !  Works with new routing routines and routes flows through a regional river netowrk
    !  Updated: Sarah Shannon to include degree day snow and glacier melting
    ! (2020)
!=================================================================
    subroutine topmod (nac,     &
        nstep,                  &
        num_rivers,             &
        acc,                    &
        dt,                     &
        dtt,                    &
        dyna_hru,               &
        dyna_riv,               &
        node_to_flow_mapping,   &
        ntt,                    &
        pe_step,                &
        q,                      &
        r_gau_step,             &
        rivers,                 &
        route_riv,              &
        route_tdh,              &
        smax,                   &
        srmax,                  &
        szm,                    &
        t0dt,                   &
        td,                     &
        wt,                     &
        qobs_riv_step,          &
        sim_id,                 &
        out_dir_path,           & 
        temp_step,              &
        elev_climate,           &
        min_ddf_snow,           &
        max_ddf_snow,           &
        tlapse,                 &
        plapse,                 &
        rain2snow_temp,         &
        melt_temp,              &
        tlag_snow,              &
        year,                   &
        month,                  &
        day,                    &
        julian_day,             &
        min_ddf_ice,            &
        max_ddf_ice,            &
        tlag_ice,               &
        snow2ice,               &
        snow_ice_sublim,        &
        rcf,                    &
        scf)                   

        use dyna_common_types
        use dyna_climate
        use dyna_snowmelt
        use dyna_glacier
        use dyna_glacier_flux
        use dyna_dynamic_dist
        use dyna_root_zone1
        use dyna_unsat_zone1
        use dyna_results_ac
        use dyna_river
        use dyna_results_ts
        use dta_route_processing
        use dyna_satzone_analyticalsolver
        use dyna_nse
        use dyna_file_open    
          
        implicit none

        ! Argument Declares
        integer :: nac
        integer :: nstep
        integer :: num_rivers
        doubleprecision :: acc
        doubleprecision :: dt
        doubleprecision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        integer :: ntt
        integer :: node_to_flow_mapping(:)

        double precision, allocatable, dimension(:,:) :: pe_step
        double precision, allocatable, dimension(:,:) :: r_gau_step
        double precision, allocatable, dimension(:,:) :: temp_step
        double precision, allocatable, dimension(:)   :: elev_climate
        integer, dimension(:), allocatable            :: julian_day
        integer, dimension(:), allocatable            :: year
        integer, dimension(:), allocatable            :: month
        integer, dimension(:), allocatable            :: day
        
        double precision :: q(num_rivers, nstep)
        doubleprecision :: rivers(nac, num_rivers)

        ! Parameters
        doubleprecision, dimension(:) :: smax
        doubleprecision, dimension(:) :: srmax
        doubleprecision, dimension(:) :: szm
        doubleprecision, dimension(:) :: t0dt
        doubleprecision, dimension(:) :: td
        doubleprecision, dimension(:) :: min_ddf_snow
        doubleprecision, dimension(:) :: max_ddf_snow
        doubleprecision, dimension(:) :: tlapse
        doubleprecision, dimension(:) :: plapse
        doubleprecision, dimension(:) :: rain2snow_temp
        doubleprecision, dimension(:) :: melt_temp
        doubleprecision, dimension(:) :: tlag_snow
        doubleprecision, dimension(:) :: min_ddf_ice
        doubleprecision, dimension(:) :: max_ddf_ice
        doubleprecision, dimension(:) :: tlag_ice
        doubleprecision, dimension(:) :: snow2ice
        doubleprecision, dimension(:) :: snow_ice_sublim
        doubleprecision, dimension(:) :: rcf
        doubleprecision, dimension(:) :: scf
        doubleprecision, parameter    :: nspin = 9
        doubleprecision, parameter    :: nse_threshold = 0.3
        
         
        doubleprecision  :: wt
        double precision :: qobs_riv_step(num_rivers, nstep)
        doubleprecision  :: nse(num_rivers)
        integer :: sim_id
        character(900) :: out_dir_path
        
        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh

        ! Local Declares
        integer :: ia
        integer :: it
        doubleprecision :: sat_ac

        integer :: i
        
        character(len=1024) :: line_format
                
        ! End declares

        dtt = dt / dble (ntt)

        !==========================================================
        !  START LOOP ON TIME STEPS (MAIN MODEL LOOP)
        !=========================================================
        
        do it = 1, nstep
           
            !  INITIALISE QB FOR EVERY RIVER REACH FOR THIS TIMESTEP
            !
            dyna_riv%qof = 0
            dyna_riv%qb  = 0
            dyna_hru%qin = 0

            
            !  CALL THE RAIN ROUTINE TO SPLIT UP THE RAIN INPUTS
            !  AND DISTRIBUTE THE EVAPORATION DATA, and temp if 
            !  required
            call climate (nac,  &
                dyna_hru,       &
                it,             &
                pe_step,        &
                r_gau_step,     &
                temp_step,      &
                elev_climate)

            !================================================================
            !
            !  FIRST PART OF THE DYNAMIC TOPMODEL FORMULATION
            !
            !  Distribute current QBF values to other elements
            !  Old time step QBF used to form new time step QIN values
            !
            !  CALL DYNAMIC_DIST() TO DO THESE CALCULATIONS
            !===========================================================
            !
            call dynamic_dist (nac, &
                num_rivers, &
                dtt, &
                dyna_hru, &
                dyna_riv, &
                rivers)

            !=================================================================
            !  START LOOP ON Fractional Area Cells
            !=================================================================

            sat_ac = 0.D0
              
            do ia = 1, nac
              
               !=============================================================
               !  Degree day snow melt model. This adds snow melt to the 
               !  precip field
               !=============================================================
            
               if (dyna_hru(ia)%ims_snow.eq.1) then
                  call snowmelt (nac,  &
                       dyna_hru,       &
                       it,             &
                       ia,             &
                       dt,             &
                       min_ddf_snow,   &
                       max_ddf_snow,   &
                       tlapse,         &
                       plapse,         &
                       rain2snow_temp, &
                       melt_temp,      &
                       tlag_snow,      &
                       snow_ice_sublim,&
                       rcf,            &
                       scf,            &
                       nstep,          &
                       julian_day)

                  !=============================================================
                  ! Degree day glacier melt. Glacier melt is added to 
                  ! the precip field 
                  !=============================================================
                 
                  if (dyna_hru(ia)%ims_glacier.eq.1) then
                     call glacier (nac,    &
                          dyna_hru,        &
                          it,              &
                          ia,              &
                          dt,              &
                          min_ddf_ice,     &
                          max_ddf_ice,     &
                          snow2ice,        &
                          melt_temp,       &
                          tlag_ice,        &
                          snow_ice_sublim, & 
                          nstep,           &
                          julian_day)      
                    end if  
                  
                 
                 
               
               end if !--- end if snow switch is on
       
                !--- write out hru fields
              ! if (month(it).eq.5.and.day(it).eq.1.and.year(it).ge.2007) then
              !      write(43,*) ia, dyna_hru(ia)%elev, dyna_hru(ia)%temp, dyna_hru(ia)%p, &
              !      dyna_hru(ia)%ep, dyna_hru(ia)%snowdepth, dyna_hru(ia)%glacierdepth, &
              !      dyna_hru(ia)%init_glacierdepth*0.917 
              !  end if  
                
               
               dyna_hru(ia)%quz  = 0.D0
               dyna_hru(ia)%uz   = 0.D0
               dyna_hru(ia)%ex   = 0.D0
               dyna_hru(ia)%exs  = 0.D0
               dyna_hru(ia)%exus = 0.D0
               
                
                !=============================================================
                !  CALL THE APPROPRIATE ROOT ZONE FORMULATION (1 is standard)
                !  THIS CALCULATES THE MAIN P OR EVAP INPUTS/OUTPUTS TO SRZ
                !=============================================================
               if (dyna_hru(ia)%ims_rz.eq.1) then

                  call root_zone1 (nac, &
                       dyna_hru, &
                       ia, &
                       it, &
                       srmax)
                  
               end if
                !
                !=============================================================
                !  CALL THE APPROPRIATE SUZ ZONE FOMULATION (1 is standard)
                !  THIS CALCULATES QUZ AND CHANGES IN THE SUZ STORE
                !=============================================================

                if (dyna_hru(ia)%ims_uz.eq.1) then

                    call unsat_zone1 (nac, &
                        dt, &
                        dtt, &
                        dyna_hru, &
                        ia, &
                        it, &
                        td)

                end if

                !=============================================================
                !  CALL THE KINEMATIC SOLUTION
                !  Implicit time stepping. Iterative solution of kinematic equation
                !  for downslope flows at this time step
                !=============================================================
                

                    call satzone_analyticalsolver (nac, &
                        dtt, &
                        dyna_hru, &
                        ia, &
                        it, &
                        ntt, &
                        smax, &
                        szm)

                
                !=============================================================
                !  CALL THE FRACTIONAL AREA RESULTS SUBROUTINE
                !=============================================================
                !
                call results_ac (nac, &
                    num_rivers, &
                    dyna_hru, &
                    dyna_riv, &
                    ia)

                !==========================================================
                !  END OF Fractional area LOOP
                !==========================================================
              
                end do !---end hru loop

                !=============================================================
                !  CALL THE RIVER ROUTING SUBROUTINE
                !=============================================================

                call river (dyna_riv,      &
                     nstep,                &
                     num_rivers,           &
                     it,                   &
                     node_to_flow_mapping, &
                     q,                    &
                     route_riv,            &
                     route_tdh)

                !=============================================================
                !  Calc NSE. If NSE < threshold then exit and do the 
                !  next mc_sim 
                !=============================================================
                
                if (it ==  5476) then                
                    call calc_nse(it,    &
                        nstep,           &
                        num_rivers,      &
                        q,               &
                        qobs_riv_step,   &    
                        nse)
                                
                      !--- Only open output files if the NSE is good.
                     
                      if (any(nse > nse_threshold)) then
                            print *, 'open output file NSE is', sim_id, nint(nse * 1000.0) * 1E-3
                            call file_open (sim_id,out_dir_path)
                      elseif (any(nse < nse_threshold)) then
                            print *, 'bad NSE go to next sim', sim_id, nint(nse * 1000.0) * 1E-3
                            exit
                     end if 
                    
                 end if
                
                !=============================================================
                !  CALL THE TIMESTEP RESULTS SUBROUTINE - CHECKS WATER BALANCE
                !=============================================================

                call results_ts (dtt, &
                     nac,             &
                     num_rivers,      &
                     dyna_hru,        &
                     dyna_riv)

                !=============================================================
                ! Sarah: Update glacier thickness outside of hru loop  
                ! Removed glacier thickness adjustment following   
                ! discussions about hrus being non-contiguous. 
		! This can possibly be picked up in future versions                 
                !============================================================= 
                if (any(dyna_hru(:)%ims_glacier.eq.1)) then
                    if (it >=  3651) then 
                    call glacier_flux (nac,     &
                          dyna_hru,         &
                          it,               &
                          dt,               &
                          nstep,            &
                          year,             &
                          month,            &
                          day,              &
                          julian_day)
                    end if          
                end  if
                
    
                !==========================================================
                !    END OF TIME STEP LOOP
                !==========================================================     
               
                end do !--- end time loop
                
                return

            end subroutine

        end module dyna_topmod
