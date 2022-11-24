module dyna_climate
contains
!
!===============================================================
!  ROUTINE TO DISTRBITUE RAIN, PET and AIR TEMP TO EACH HRU
!===============================================================

  subroutine climate (nac, &
       dyna_hru, &
       it, &
       pe_step, &
       r_gau_step,&
       temp_step,&
       elev_climate)

    use dyna_common_types
    implicit none

    integer :: nac
    type(dyna_hru_type) :: dyna_hru(nac)
    integer :: it
    double precision, dimension(:,:) :: pe_step
    double precision, dimension(:,:) :: r_gau_step
    double precision, dimension(:,:) :: temp_step
    double precision, dimension(:)   :: elev_climate

    ! Local Declares
    integer :: ia
    
    do ia = 1, nac

       dyna_hru(ia)%p  = r_gau_step (dyna_hru(ia)%ippt, it) !0.001!0.0007508700!
       dyna_hru(ia)%ep = pe_step (dyna_hru(ia)%ipet, it)

       if (size(temp_step) > 1) then
          dyna_hru(ia)%temp  = temp_step (dyna_hru(ia)%itemp, it)
          
       end if
       
       !--- elev of forcing data has no time dimension.
       !--- We are assuming the climate elev grid is the same as the 
       !--- temperature grid
       if (it == 1.and. size(elev_climate) > 0 ) then
          dyna_hru(ia)%elev_climate = elev_climate(dyna_hru(ia)%itemp)
       end if

    end do

   
  end subroutine climate
        
 end module dyna_climate
