module dyna_glacier_flux
contains
!==========================================================
!--- Output fraction of catchment with 
!--- glacier depth > threshold ice depth (meters of ice)
!==========================================================
  subroutine glacier_flux (nac, &
       dyna_hru,        &
       it,              &
       dt,              &
       nstep,           &
       year,            &
       month,           &
       day,             &
       julian_day)
         
    use dyna_common_types
    
    implicit none
    integer              :: nac
    type(dyna_hru_type)  :: dyna_hru(nac)
    integer              :: it
    integer              :: nstep 
    integer              :: year(nstep)
    integer              :: month(nstep)
    integer              :: day(nstep)
    integer              :: julian_day(nstep) 
    integer              :: ia, i, iw
    double precision     :: dt
    
    
    double precision     :: sumglacier_area_mmm
    double precision     :: sumglacier_area_mm
    double precision     :: sumglacier_area_cm
    double precision     :: sumglacier_area_m 
    double precision     :: wtranstot 
    double precision, parameter  :: cellsize =  74.0 ! hardcoding this for now
    integer, allocatable, dimension (:) :: ghru
    real, parameter      :: dice = 0.917  !density of ice as fraction of water
      
    !if (it == 1) then 
        sumglacier_area_mmm = 0.0
        sumglacier_area_mm  = 0.0
        sumglacier_area_cm  = 0.0
        sumglacier_area_m   = 0.0
    !end if 
    
    if (any(dyna_hru(1:nac)%ims_glacier == 1 )) then
        allocate(ghru(count(dyna_hru%ims_glacier.eq.1)))
        i = 1
         do ia = 1,nac
            if (dyna_hru(ia)%ims_glacier.eq.1) then
                ghru(i) = ia
                i = i + 1
            end if
        end do
    end if

!--- Fractional area of catchment with glacier depth > ice depth threshold               
    do ia = 1,size(ghru)   
        if  ( dyna_hru(ghru(ia))%glacierdepth/dice  > 0.000001) then      
            sumglacier_area_mmm = sumglacier_area_mmm + dyna_hru(ghru(ia))%ac
        end if  
        if ( dyna_hru(ghru(ia))%glacierdepth/dice  > 0.001) then     
            sumglacier_area_mm = sumglacier_area_mm + dyna_hru(ghru(ia))%ac
        end if
        if ( dyna_hru(ghru(ia))%glacierdepth/dice > 0.01) then     
            sumglacier_area_cm = sumglacier_area_cm + dyna_hru(ghru(ia))%ac
        end if  

    end do !--- end loop through ghrus
   
    !--- write out annual glaciated area fraction [Naryn conversion to km2 * 57892.7] 
    write(42, *) it, sumglacier_area_mmm, sumglacier_area_mm, sumglacier_area_cm
    
    end subroutine glacier_flux

end module dyna_glacier_flux
