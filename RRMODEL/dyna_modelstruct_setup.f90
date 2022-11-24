module dyna_modelstruct_setup
contains
    !
    !===============================================================
    !  ROUTINE FOR SETTING UP MODEL STRUCTURE PER HRU
    !===============================================================
    !
    subroutine modelstruct_setup (dyna_hru, &
         nac,                               &
         temp_file,                         &
         elev_file,                         &
         initglacier_file)

        use dyna_common_types

        implicit none

        ! Argument Declares

        type(dyna_hru_type), dimension(:), allocatable :: dyna_hru
        integer :: nac
        character(900),  dimension(:), allocatable :: temp_file
        character(900),  dimension(:), allocatable :: elev_file
        character(900),  dimension(:), allocatable :: initglacier_file
        
        ! Local Declares
        character(len=20), dimension(:), allocatable :: ms_names
        integer, dimension(:,:), allocatable :: ms_values
        integer :: num_ms_types, num_ms_names
        integer :: ms_rz_found, ms_uz_found, ms_snow_found, ms_glacier_found 
        integer :: i, j
        
        ! End declares

        ! Skip header lines
        read(11,*)
        read(11,*)

        ! Read in the number of model structure types and check this is the same as read in the from the HRU meta files
        read(11,*) num_ms_types, num_ms_names

        !  Make sure that the number of different parameter types from HRU Meta File matches number in parameter file
        if (maxval(dyna_hru%ims).ne.num_ms_types) then
            print *, 'Number of model structure types must match the number given from the HRU analysis'
            print *, 'Num_ms_types from model structure file = ', num_ms_types
            print *, 'Number of model structure types from HRU Meta File = ', maxval(dyna_hru%ims)
            STOP
        end if

        !allocate param_names list and param_values array, now we know number of params
        allocate(ms_names(num_ms_names+1))
        allocate(ms_values(num_ms_types, num_ms_names))

        !read in list of parameter names
        read(11, *) ms_names

        ! Read in whole parameter values table
        do i = 1, num_ms_types
            read(11, *) j, ms_values(i,:)
        end do

       
        !loop through model structure list, checking that each model component name corresponds to an existing component
        !if components exist, add them dyna_Hru structure
        ! NEW MODEL COMPONENTS MUST BE ADDED TO THE LIST AT THE BEGINNING OF THIS MODULE

        ! First check that all

        ! Check for Root Zone - First set to default option (rootzone = 1)

        dyna_hru%ims_rz      = 1
        dyna_hru%ims_uz      = 1
        dyna_hru%ims_snow    = 0 
        dyna_hru%ims_glacier = 0
        
        ms_rz_found=0

        do i = 2, num_ms_names+1

            select case (ms_names(i))

                case ('rootzone')

                    do j = 1, nac
                        dyna_hru(j)%ims_rz = ms_values(dyna_hru(j)%ims, i-1)
                    end do
                    ms_rz_found = 1

                case ('unsatzone')

                    do j = 1, nac
                        dyna_hru(j)%ims_uz = ms_values(dyna_hru(j)%ims, i-1)
                    end do
                    ms_uz_found = 1

                 case ('snow')
                    
                    do j = 1, nac
                       dyna_hru(j)%ims_snow = ms_values(dyna_hru(j)%ims, i-1)
                    end do
                    ms_snow_found = 1

                 case ('glacier')
                    
                    do j = 1, nac
                       dyna_hru(j)%ims_glacier = ms_values(dyna_hru(j)%ims, i-1)
                    end do
                    ms_glacier_found = 1
                    
                 case default

                print *, 'Model Structure Component ', trim(ms_names(i)), ' is not recognised'
                print *, 'The following model structure components can be set in the model structure file:'
                print *, 'rootzone  unsatzone  snow'

            end select

         end do

         if (any(dyna_hru(1:nac)%ims_glacier == 1 ).and.any(dyna_hru(1:nac)%ims_snow == 0 ) ) then
            print *, 'Glacier model needs to have snow model switched'  // &
                 ' on, stopping simulation'
            stop
         end if
            
         !--- check temp and elev files are provided for the snow model
         if (any(dyna_hru(1:nac)%ims_snow == 1 ) ) then
            if (any(temp_file(:)=="" )) then
               print *, 'Snow model switch is on but temperature file is'  // &
                    ' not provided, stopping simulation'
               stop
            elseif (any(elev_file(:)=="" )) then
               print *, 'Snow model switch is on but climate elevation'  // &
                    ' file is not provided, stopping simulation'
               stop
            end if
         end if
         
         !--- THIS DOES NOT WORK YET. COME BACK TO THIS 
         !--- check temp, elev files and initial glacier thickness 
         !--- are provided for the glacier model
         if (any(dyna_hru(1:nac)%ims_glacier == 1 ) ) then
            if (any(initglacier_file(:)=="" )) then
               print *, 'Glacier model switch is on but initial glacier'  // &
                    ' depth is not provided, stopping simulation'
               stop
            elseif (any(elev_file(:)=="" )) then
               print *, 'Glacier model switch is on but climate elevation'  // &
                    ' file is not provided, stopping simulation'
               stop
            elseif (any(temp_file(:)=="" )) then
               print *, 'Glacier model switch is on but temperature file is'  // &
                    ' not provided, stopping simulation'
               stop
            end if
         end if
         
         
        if (ms_rz_found.ne.1) then !i.e. user specified a parameter name that does not exist - print a warning and stop program
            print *,
            print *, 'WARNING: rootzone has not been specified in model structure file'
            print *, 'All HRUs will have the default rootzone model component'
        end if

        if (ms_uz_found.ne.1) then !i.e. user specified a parameter name that does not exist - print a warning and stop program
            print *,
            print *, 'WARNING: unsatzone has not been specified in model structure file'
            print *, 'All HRUs will have the default unsatzone model component'
        end if

        if (ms_snow_found.ne.1) then 
            print *,
            print *, 'WARNING: snow switch has not been specified in model structure file'
            print *, 'All HRUs will have the default of no snow'
         else if (ms_snow_found.eq.1) then
            write(999,*) ''
            write(999,*) 'Snow melt model is switched on'
            write(999,*) ''
         end if

         if (ms_glacier_found.ne.1) then 
            print *,
            print *, 'WARNING: glacier has not been specified in model structure file'
            print *, 'All HRUs will have the default of no glaciers'
         else if (ms_glacier_found.eq.1) then
            write(999,*) ''
            write(999,*) 'Glacier model is switched on'
            write(999,*) ''
        end if
         
        write(999,*) 'Read in model structure file succesfully'

        close (11)

    end subroutine

end module dyna_modelstruct_setup
