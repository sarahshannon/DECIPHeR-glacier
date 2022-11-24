module dyna_init_glacier

contains
    
!==========================================================
!--- Initialise glacier temperature. Convert glacier 
!--- depth from meters to meters of water equivalent 
!--- Sarah Shannon: Jan 2021 
!==========================================================
    
    subroutine init_glacier(dyna_hru, &
                            nac)    
       
    use dyna_common_types

    implicit none
    
    integer :: nac
    integer :: ia
    integer :: h
    type(dyna_hru_type) :: dyna_hru(nac) 
    real, parameter :: dice = 0.917    !density of ice, fraction of water
    
    !--- set const for all HRUs
    dyna_hru%glaciertemp  = -5.0
   
    do ia = 1, nac
         !--- convert m ice to m water equival
         dyna_hru(ia)%glacierdepth = dyna_hru(ia)%init_glacierdepth*dice        
    end do  
     
    
     
  end subroutine init_glacier

end module dyna_init_glacier
    