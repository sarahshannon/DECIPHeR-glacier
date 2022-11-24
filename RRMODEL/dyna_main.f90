!                                                                       
!  =================================================================    
!  DECIPHeR VERSION 1
!
!  Bristol University 2018 (Gemma Coxon, Jim Freer, Rosie Lane and Toby Dunne)
!  Based on fortran 77 version of dynamic TOPMODEL produced in
!  Lancaster University 12/01/00 (Keith Beven & Jim Freer)              
!  Migrated to std=F2003 Toby Dunne 25/01/2015
!  Modified extensively by Gemma Coxon 2016-2018
!  Modifed by Sarah Shannon 2019 for snow and glacier melt
!  =================================================================
!                                                                       
program main

    use dyna_common_types
    use dyna_project
    use dyna_random
    use dyna_tread_dyna
    use dyna_inputs
    use dyna_mc_setup
    use dyna_modelstruct_setup
    use dyna_genpar
    use dyna_main_loop
    use dyna_file_open
    use dyna_read_routingdata

    use dta_utility
    use dta_route_processing
    use dta_riv_tree_node

    implicit none

    !  Local variable declares
    integer :: seed_1, seed_2
    integer :: i, h
    integer :: nac
    integer :: nstep
    integer :: num_rivers

    ! allocated from 'Input'
    double precision, dimension(:,:), allocatable :: pe_step
    double precision, dimension(:,:), allocatable :: temp_step
    double precision, dimension(:),   allocatable :: qobs_riv_step_start
    double precision, dimension(:,:),   allocatable :: qobs_riv_step
    double precision, dimension(:,:), allocatable :: r_gau_step
    integer,          dimension(:),   allocatable :: year
    integer,          dimension(:),   allocatable :: month
    integer,          dimension(:),   allocatable :: day
    integer,          dimension(:),   allocatable :: julian_day
    double precision, dimension(:),   allocatable :: elev_climate

    ! allocated after 'Input'
    type(dyna_hru_type), dimension(:), allocatable :: dyna_hru
    type(dyna_riv_type), dimension(:), allocatable :: dyna_riv

    character(900),  dimension(:), allocatable :: temp_file
    character(900),  dimension(:), allocatable :: elev_file
    character(900),  dimension(:), allocatable :: initglacier_file
    
    
    ! MULTI POINT RIVER ROUTING
    ! Toby Dunne April 2016 + GC June 2016
    ! these type are defined in modules as part of the dta files
    type(route_river_info_type) :: route_riv
    type(route_time_delay_hist_type) :: route_tdh
    integer :: cat_route_vmode
    integer :: routing_mode
    integer, dimension(:), allocatable :: node_to_flow_mapping
    character(900) :: out_dir_path

    ! Monte Carlo Simulations
    integer :: numsim, i_mc, num_par_types
    integer, allocatable, dimension(:) :: num_mcpar
    integer, parameter  :: nparam =  23 
    character(len=20), dimension(nparam) :: all_pm_names !names of all possible parameters - must be added to if new params are developed!
    doubleprecision, allocatable, dimension(:,:) :: mcpar_ll
    doubleprecision, allocatable, dimension(:,:) :: mcpar_ul
    doubleprecision, allocatable, dimension(:,:) :: mcpar
    
    !--- lhs params
    !doubleprecision, allocatable, dimension(:,:) :: lhs_data
    doubleprecision, dimension(20)               :: dummy

    ! Kinematic solution parameters and time step
    doubleprecision :: dt
    doubleprecision :: acc
    doubleprecision :: wt
    integer :: ntt
    doubleprecision :: dtt

    double precision, dimension(:,:), allocatable :: rivers !(nac, n_riv)
    double precision, dimension(:), allocatable :: sum_ac_riv  !

    character(1024) :: arg
    character(1024) :: auto_start_file
    character(1024) :: lhs_params_filechar
    integer         :: lhs_params_filenum, sim_id
    
    auto_start_file = ''
    i = 0

    print *, '--- Starting DECIPHeR rainfall-runoff modelling ---'
    
    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-auto')) then
            CALL get_command_argument(i+1, auto_start_file)
        endif
        i = i + 1
    end do

    !--- Sarah hacking in here to batch run on bp
    !--- pass LHS params file via the command line arg
    call get_command_argument(1, lhs_params_filechar)
    
    read(lhs_params_filechar,*)lhs_params_filenum
    print *, 'lhs_params_filenum',lhs_params_filenum
    
    !  Call project to read in project files
    call project (auto_start_file, &
         out_dir_path,             &
         temp_file,                &
         elev_file,                &
         initglacier_file)
    
    write(999,*) ''
    write(999,*) 'Reading in HRU files'

    ! Call tread_dyna to read in the HRU file
    call tread_dyna (nac, &
        num_rivers, &
        dyna_hru, &
        dyna_riv, &
        rivers, &
        sum_ac_riv)

    write(999,*) ''
    write(999,*) 'Reading in Model Structure files'

    ! Read in model structure info - this has been moved to before inputs subroutine 
    ! we can error check that temp and elev are provided for the snow model
    ! before going any further
    call modelstruct_setup (dyna_hru, &
         nac,                         &
         temp_file,                   &
         elev_file,                   &
         initglacier_file)

    write(999,*) ''
    write(999,*) 'Reading in Input files'

    ! Call inputs to read in the input data
    call inputs (nstep,         &
        num_rivers,             &
        dt,                     &
        pe_step,                &
        qobs_riv_step_start,    &
        qobs_riv_step,          &
        r_gau_step,             &
        temp_step,              &
        year,                   &
        month,                  &    
        day,                    &
        julian_day,             &
        elev_climate)
        

    write(999,*) ''
    write(999,*) 'Reading in Parameter files'

        ! Read in parameter data
    call mc_setup (ACC, &
        dyna_hru, &
        mcpar_ll, &
        mcpar_ul, &
        NTT, &
        num_mcpar, &
        num_par_types, &
        numsim, &
        seed_1, &
        seed_2, &
        WT, &
        all_pm_names)

    !  Initialise randowm number generator
    !
    call ranin (seed_1, seed_2)

    write(999,*) ''
    write(999,*) 'Reading in Model Structure files'

    
    write(999,*) ''
    write(999,*) 'Reading in Routing Data'

    ! call read_routingdata to read in the routing data
    call read_routingdata (cat_route_vmode, &
        dyna_hru, &
        nac, &
        node_to_flow_mapping, &
        num_rivers, &
        route_riv, &
        route_tdh, &
        routing_mode)

    write(999,*) ''
    write(999,*) 'Looping through simulations...'

    !--- Loop through for each simulation
    
    open (unit=99, file=lhs_params_filechar, status='old', action='read')
    
    do i_mc = 1, numsim
 
        ! Get a new set of random parameters
        call genpar (mcpar, &
            mcpar_ll, &
            mcpar_ul, &
            num_mcpar, &
            num_par_types)
        
        !--- override the mc params with lhs params
        read(99, *) sim_id, dummy
        
        do h = 1, 1453
            mcpar(h,:) = dummy
        end do
        
        print *, 'sim_id', sim_id
        
        !--- name output file with sim_id for calibrating on bluecrystal
        !call file_open (i_mc,out_dir_path)
        !call file_open (sim_id,out_dir_path)
	
        !--- end lhs params
                                                        
        call mainloop (all_pm_names,&
            nac,                    &
            nstep,                  &
            num_rivers,             &
            acc,                    &
            dt,                     &
            dtt,                    &
            dyna_hru,               &
            mcpar,                  &
            ntt,                    &
            node_to_flow_mapping,   &
            num_par_types,          &
            pe_step,                &
            qobs_riv_step_start,    &
            qobs_riv_step,          &
            sim_id,                 &
            out_dir_path,           &
            r_gau_step,             &
            rivers,                 &
            route_riv,              &
            route_tdh,              &
            sum_ac_riv,             &
            wt,                     &
            temp_step,              &
            elev_climate,           &
            year,                   &
            month,                  &    
            day,                    &
            julian_day)

        close (40)
        close (41)
        
    end do
    close(99)
    close(999)
    !
    !============================================================
    !  END MAIN MC LOOP
    !============================================================
    !
    !stop

end program main
