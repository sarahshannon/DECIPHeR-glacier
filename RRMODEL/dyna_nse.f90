module dyna_nse

contains
    
!==========================================================
!--- Calculate NSE. Use this to abort to next simulation  
!--- if the NSE is less than threshold. Hard coded time 
!--- steps for which to calculate NSE over. Make 
!--- this more flexible later, this is hack for now 
!--- Sarah Shannon 2021
!==========================================================
    
    subroutine calc_nse(it,         &
                   nstep,           &
                   num_rivers,      &
                   q,               &
                   qobs_riv_step,   &    
                   nse)
    
 
    implicit none
    integer :: it
    integer :: nstep
    integer :: num_rivers
    double precision :: qobs_riv_step(num_rivers, nstep)
    double precision :: q(num_rivers, nstep)
    doubleprecision  :: nse(num_rivers) 
    integer, parameter  :: it_start   = 3650  !--- after 10 year spinup
    integer, parameter  :: it_end     = 5475  !--- year 15 spinup 
    integer :: n, r
    
    double precision, dimension(:), allocatable :: qobs
    double precision, dimension(:), allocatable :: qsim
    double precision, dimension(:), allocatable :: top
    double precision, dimension(:), allocatable :: bottom
    
    !--- calc nse over time steps it_start --> it_end
    !--- choose which catchment. here I use upper Naryn
    n = it_end - it_start + 1
    
    do r = 1, num_rivers    
        
        qobs = qobs_riv_step(r,it_start:it_end)
        
        qsim = q(r,it_start:it_end)
    
        top  = qobs - qsim
    
        bottom = qobs -  (sum(qobs) / n) 
    
        nse(r) = 1 - (sum(top*top)/sum(bottom*bottom))
    end do 
     
  end subroutine calc_nse

end module dyna_nse
    