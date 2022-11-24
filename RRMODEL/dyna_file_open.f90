module dyna_file_open
contains
    !
    !===============================================================
    !  ROUTINES FOR OPENING THE OUTPUT FILES NEEDED PER SIMULATION RUN
    !===============================================================
    !
    subroutine file_open (i_mc, &
        out_dir_path)

        implicit none

        ! Argument Declares
        integer :: i_mc
        character(900) :: out_dir_path

        ! Local Declares
        character(900) :: filename
        character(20) :: mc_id_val

        ! End declares

        ! Now open all the files needed
        ! First get the padded i_mc values into characters for
        ! all the filenames to open
        write (mc_id_val, "('mc_id_',I8.8)") i_mc

        ! HYDROL FILES
        ! 1) Then get the filename for the .res file
        filename = trim(out_dir_path) // &
            trim(mc_id_val) // &
            '.res'
        ! Then open this .res file
        open(unit=40,file=filename,status='unknown')

        ! 2) Then get the filename for the .flow file
        filename = trim(out_dir_path) // &
            trim(mc_id_val) // &
            '.flow'
        ! Then open this .flow file
        open(unit=41,file=filename,status='unknown')

        ! 3) filename for glacier area and vol file
        filename = trim(out_dir_path) // &
            trim(mc_id_val) // &
            '_glacier.dat'
        open(unit=42,file=filename,status='unknown')

        ! 4) filename for glacier depth
        
        !--- binary output is less than half size of ascii
        !filename = trim(out_dir_path) // &
        !    trim(mc_id_val) // &
        !    '_hru_snowdepth.bin'
        !open(unit=43,file=filename,form='unformatted',access='stream',status='replace')

        !filename = trim(out_dir_path) // &
        !    trim(mc_id_val) // &
        !    '_hru_glacier.bin'
        !open(unit=44,file=filename,form='unformatted',access='stream',status='replace')
         !5) filename for hru precip
        !filename = trim(out_dir_path) // &
        !    trim(mc_id_val) // &
        !    '_hru_fields.dat'
        !open(unit=43,file=filename,status='replace')

        return
    end subroutine file_open

end module dyna_file_open
