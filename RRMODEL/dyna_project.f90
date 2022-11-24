module dyna_project

contains
    !
    !===============================================================
    !  ROUTINES FOR READING IN PROJECTS BASIC DATA
    !  THIS SETS UP ALL THE DIRECTORIES AND FILENAMES TO OPEN
    !  IT ALSO SETS SOME IMPORTANT FLAGS SO THAT THE PROGRAM KNOWS
    !  WHAT TO RUN
    !  Modified by Sarah Shannon to open air temp and forcing elev
    !  if files are provided 
    !===============================================================

    subroutine project (auto_start_file, &
         out_dir_path,                   &
         temp_file,                      &
         elev_file,                      &
         initglacier_file)

      use dyna_common_types

      implicit none

      ! Argument Declares
      character(1024), intent(in) :: auto_start_file
      
      ! Local Declares
      integer :: allocatestatus, deallocatestatus
      integer :: i
      integer :: internal_file_count
      logical :: file_exists
      ! End declares

      parameter (internal_file_count = 50)
      
      ! temp variable to store filename strings
      character(900) :: char_val
      character(900) :: file_open (internal_file_count)
      character(900) :: mstruct_file,mc_par_file, out_dir_path

      integer :: id_project
      integer :: proj_auto
      integer :: id_cat_file
      integer :: id_flow_file

      ! Allocatable arrays - project files
      integer :: num_proj
      character(900), dimension(:), allocatable :: proj_name
      character(900), dimension(:), allocatable :: proj_dir
      character(900), dimension(:), allocatable :: proj_file       !filemanager file (old project file)
      character(900), dimension(:), allocatable :: proj_settings   !settings file for specified project
      character(900), dimension(:), allocatable :: proj_outdir
      character(900)                            :: comment

      ! Allocatable arrays - hru files
      integer :: num_cat_input
      character(900), dimension(:), allocatable :: cat_comment
      character(900), dimension(:), allocatable :: cat_file
      character(900), dimension(:), allocatable :: cat_meta_file
      character(900), dimension(:), allocatable :: cat_route_flow_conn
      character(900), dimension(:), allocatable :: cat_route_riv_data
      character(900), dimension(:), allocatable :: cat_route_flow_point
      character(900), dimension(:), allocatable :: initglacier_file
      
      ! Allocatable arrays - input files
      integer :: num_flow_input
      character(900),  dimension(:), allocatable :: flow_file_comment
      character(900),  dimension(:), allocatable :: flow_file
      character(900),  dimension(:), allocatable :: evap_file
      character(900),  dimension(:), allocatable :: precip_file
      character(900),  dimension(:), allocatable :: temp_file
      character(900),  dimension(:), allocatable :: elev_file
      doubleprecision, dimension(:), allocatable :: flow_err

      integer :: ioerr
      character(1024) :: tmp_char
      
      character(1024) :: proj_id, fname
      
      do i = 1, internal_file_count
         ! initialise to empty string
         file_open (i) = ''
      end do

      !===================================================================
      !  PART 1: READING IN THE LIST OF PROJECTS
      !
      !  READ IN THE PROJECT FILE + NUMBER OF PROJECTS
      !===================================================================

      !char_val = 'project.dat' this was the default 
      !--- read in project file using command line arg. 
      call get_command_argument(2, proj_id)
      
      fname = 'project_' // trim(proj_id) // '.dat' 
      print *, fname
      
      !write(999,*) 'Open:', trim(char_val)
      !open (10, file = char_val, status = 'old')
      open (10, file = fname, status = 'old')
      
      ! Read in the number of projects available
      read(10,*)
      read (10, * ) num_proj

      ! Check to see whether this is an automated read.
      if(len_trim(auto_start_file) > 0) then
         proj_auto = 1
      else
         proj_auto=0
         ! auto_start_file = 'PROJECTS/auto.dat'
      endif

      ! Allocate storage for arrays
      allocate(proj_name (num_proj), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(proj_dir (num_proj), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(proj_file (num_proj), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(proj_settings (num_proj), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(proj_outdir (num_proj), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"

      ! For each project - read in project info
      do i = 1, num_proj

         ! Read in project name, directory, file, where to put output
         read (10,*) !empty line before each different project
         read (10,525) comment, proj_name(i)
         read (10, 525) comment,proj_dir (i)
         read (10, 525) comment,proj_file (i)
         read (10, 525) comment,proj_outdir (i)
         print *, i, trim(proj_name(i)), trim(proj_dir(i))

      end do
                                                                        
      !  Read in the .auto data if required for batch runs

      if (proj_auto.eq.1) then
         print *, 'Open:', trim(auto_start_file)
         open (10, file = auto_start_file, status = 'old')
         read (10, * ) id_project
         read (10, * ) id_cat_file
         read (10, * ) id_flow_file
      else
         PRINT *, 'Input the project you wish to use for this run?'
         !read (*,*) id_project
         id_project = 1 !sarah 
      endif

      !tmp_char = trim(proj_dir(id_project))//'/OUTPUT/DYMOND.log'
	
      !--- remove the line which forces the user to have a dir called OUTPUT
      tmp_char = trim(proj_dir(id_project))//'logfile.dat'
      open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

      if(ioerr/=0) then
         print*,'error opening output file: ', trim(tmp_char)
         print*,'ensure the directory exists and correct write permissions are set'
         stop
      endif

      write(999,*) '--- DYMOND ---'
      write(999,*) ''
      write(999,*) 'Reading in project files'

      close (10)    !
      !===================================================================
      !  PART 2: READ IN THE CHOSEN PROJECT FILE
      !===================================================================
      !
      char_val = trim(proj_dir(id_project))//trim(proj_file(id_project))//'.dat'
      write(999,*) 'Open 1:', char_val
      open (10, file = char_val, status = 'old')

      read (10, *) !skip header lines
      read (10, *)
      read (10, 525) comment,mc_par_file
      read (10, 525) comment,mstruct_file
      read (10, 530) comment,num_cat_input
      read (10, 530) comment,num_flow_input
      read(10,*)      !comment lines
      read(10,*)

      !  OPEN UP DIFFERENT DYNATOP/TOPMOD ATB DISTRIB FILES
      allocate(cat_comment (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(cat_file (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(cat_meta_file (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(cat_route_flow_conn (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(cat_route_riv_data (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(cat_route_flow_point (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(initglacier_file (num_cat_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      print *, 'The following HRU files are available:'

      do i = 1, num_cat_input

         read (10,*) !blank line above each catchment input
         read (10, 525) comment,cat_comment (i)
         read (10, 525) comment,cat_file (i)
         read (10, 525) comment,cat_meta_file (i)
         read (10, 525) comment,cat_route_flow_conn(i)
         read (10, 525) comment,cat_route_riv_data(i)
         read (10, 525) comment,cat_route_flow_point(i)
         read (10, 525) comment,initglacier_file(i)
         print *, i, trim(cat_comment(i))

      end do

      if (proj_auto.eq.0) then
         PRINT *, 'Input the HRU file you wish to use for this run?'
         !  read ( *, * ) id_cat_file
         id_cat_file=1!sarah
      endif

      !  THEN READ IN ALL THE INPUT METADATA
      print *,
      PRINT *, 'The following input files are available:'

      allocate(flow_file_comment (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(flow_file (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(evap_file (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(precip_file (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(temp_file (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(elev_file (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      allocate(flow_err (num_flow_input), stat = AllocateStatus)
      if (allocatestatus /= 0) stop "*** Not enough memory ***"
      
      ! Read in all the input files
      read(10,*) !ignore these header lines
      read(10,*)
      
      do i = 1, num_flow_input

         read(10,*) !blank header line
         read (10, 525) comment, flow_file_comment (i)
         read (10, 525) comment, flow_file (i)
         read (10, 525) comment, precip_file (i)
         read (10, 525) comment, evap_file (i)
         read (10, 525) comment, temp_file (i)
         read (10, 525) comment, elev_file (i)
         print *, i, trim(flow_file_comment (i))
           
      end do

      if (proj_auto.eq.0) then
         PRINT *, 'Input the INPUT file you wish to use for this run?'
         ! read ( *, * ) id_flow_file
         id_flow_file = 1!sarah
      end if

      !===================================================================
      !  PART 3: OPENING THE CHOSEN PROJECT FILES & READING IN THE START.DAT
      !  INFORMATION
      !
      !  THEN USING THE PROJECT DATA SET UP THE NAMES OF THE FILES NEEDING
      !  TO BE OPENED
      !
      !===================================================================
      !
      !  OPEN UP THE MC PARAMETER FILE
      !
      file_open (2) = trim( proj_dir (id_project) ) // &
           mc_par_file
      
      write(999,*) 'Open 2:', file_open (2)
      open (12, file = file_open (2) , status = 'old')
      !
      !  OPEN UP THE MODEL STRUCTURE FILE
      !
      file_open (3) = trim( proj_dir (id_project) ) // &
           mstruct_file
      
      write(999,*) 'Open 3:', file_open (3)
      open (11, file = file_open (3) , status = 'old')
      !
      !  OPEN UP THE CAT ATB FILE
      !
      file_open (4) = trim( proj_dir (id_project) ) // &
           cat_file (id_cat_file)
      
      write(999,*) 'Open 4:', file_open (4)
      open (13, file = file_open (4) , status = 'old')
      
      !  OPEN UP THE HRU Meta FILE
      !
      file_open (5) = trim( proj_dir (id_project) ) // &
           cat_meta_file (id_cat_file)
      
      write(999,*) 'Open 5:', file_open (5)
      open (17, file = file_open (5) , status = 'old')
      !
      !
      !  JEF 30/5/16 Open up the other file names for new routing
      !  NOTE chosen 500 for open but all below needs a major
      !  re-write as will now write output per model simulation
      !  and not enough separation in unit numbers for the number
      !  of potential rivers!
      !  GC 13/06/16 Added additional network routing file

      !  Toby new network routing file - flow connectivity
      file_open (32) = trim( proj_dir (id_project) ) // &
           cat_route_flow_conn (id_cat_file)
      
      write(999,*) 'Open 6:', file_open (32)
      open (501, file = file_open (32) , status = 'old')
      
      !  Toby new network routing file - river data
      file_open (33) = trim( proj_dir (id_project) ) // &
           cat_route_riv_data (id_cat_file)
      
      write(999,*) 'Open 7:', file_open (33)
      open (502, file = file_open (33) , status = 'old')
      
      !  Toby new network routing file - flow point
      file_open (34) = trim( proj_dir (id_project) ) // &
           cat_route_flow_point (id_cat_file)

      write(999,*) 'Open 8::', file_open (34)
      open (503, file = file_open (34) , status = 'old')

    !--- Sarah initial glacier depth file     
    file_open (35) = trim( proj_dir(id_project))// &
           initglacier_file(id_cat_file)
    
    !--- do error checks for read
    if  (file_open (35) /= "") then
        inquire(file = file_open (35), exist = file_exists)
        if (file_exists) then
            write(999,*) 'Open 9::', file_open (35)
            open (504, file = file_open (35) , status = 'old') 
        else
            print *, file_open (35), ' does not exist' , &
            ' stopping simulation' 
            stop        
        end if 
    end if 
        
      !
      !  OPEN UP THE FLOW FILES(if required - which it should!!)
      !
      file_open (6) = trim( proj_dir (id_project) ) // &
           flow_file (id_flow_file)
      file_open (7) = trim( proj_dir (id_project) ) // &
           evap_file (id_flow_file)
      file_open (8) = trim( proj_dir (id_project) ) // &
           precip_file (id_flow_file)
      file_open (9) = trim( proj_dir (id_project) ) // &
           temp_file (id_flow_file)
      file_open (10) = trim( proj_dir (id_project) ) // &
           elev_file (id_flow_file)
      
         
      write(999,*) 'Open 9:', file_open (6)
      write(999,*) 'Open 10:', file_open (7)
      write(999,*) 'Open 11:', file_open (8)
      write(999,*) 'Open 12:', file_open (9)
      write(999,*) 'Open 13:', file_open (10)
      
            
      open (14, file = file_open (6) , status = 'old')
      open (15, file = file_open (7) , status = 'old')
      open (16, file = file_open (8) , status = 'old')
      
      !--- open air temp file if provided
      if  (temp_file (id_flow_file) /= "") then
         inquire(file = file_open (9), exist = file_exists)
         if (file_exists) then
            open (18, file = file_open (9) , status = 'old')
         else
            print *, 'Temperature file is specified ' // &
                 'but does not exist'
            stop
         end if
      end if

      !--- open climate forcing elev if provided
      if  (elev_file (id_flow_file) /= "") then
         inquire(file = file_open (10), exist = file_exists)
         if (file_exists) then
            open (19, file = file_open (10) , status = 'old')
         else
            print *, 'Climate elevation file is specified ' // &
                 'but does not exist'
            stop
         end if
      end if
   
   
      out_dir_path = trim(proj_dir (id_project))// &
           trim(proj_outdir (id_project))//'_'

      deallocate(cat_comment, stat = DeAllocateStatus)
      if (deallocatestatus /= 0) stop "*** Trouble deallocating ***"
      deallocate(cat_file, stat = DeAllocateStatus)
      if (deallocatestatus /= 0) stop "*** Trouble deallocating ***"

525   format(a13,1X,a700)
530   format(a13,1X,i4)

      return
    end subroutine project

end module dyna_project
