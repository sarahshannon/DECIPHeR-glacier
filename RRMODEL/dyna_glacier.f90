module dyna_glacier
contains

!===============================================================
!  Calculate glacier depth: sarah.shannon@bristol.ac.uk 2021
!  Accumulate by converting a fraction of the snowpack to 
!  glacicer ice. Melt using degree day model. Subliminate using   
!  PET. Glacier depth has units m.w.equiv 
!===============================================================

  subroutine glacier (nac, &
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
        
    use dyna_common_types
    implicit none

    integer              :: nac
    type(dyna_hru_type)  :: dyna_hru(nac)
    integer              :: it
    integer              :: nstep
    double precision     :: dt
    integer              :: julian_day(nstep)
         
    !--- local ice variables ----------------
    integer              :: ia
    double precision     :: pot_icemelt
    double precision     :: glacier_accum
    double precision     :: glacier_melt
    double precision     :: glacier_sublim
    double precision     :: ddf_ice_day
    double precision     :: beta                       ! m day-1
    double precision, dimension(:) :: min_ddf_ice      ! m w e oC-1 day-1
    double precision, dimension(:) :: max_ddf_ice      ! m w e oC-1 day-1
    double precision, dimension(:) :: melt_temp        ! oC
    double precision, dimension(:) :: snow2ice         ! basal turnover snow to ice (dimensionless) 
    double precision, dimension(:) :: tlag_ice         ! temp lag 0 - 1  (1 is no lag) 
    double precision, dimension(:) :: snow_ice_sublim  ! fraction to limit pet as sublim
    
    !--- Initialise local variables every timestep
    glacier_melt   = 0.0
    pot_icemelt    = 0.0
    glacier_accum  = 0.0
    glacier_sublim = 0.0
    ddf_ice_day    = 0.0
    beta           = 0.0 
      
    !--- previous timestep depth
    dyna_hru(ia)%glacierdepth_old = dyna_hru(ia)%glacierdepth 
    
    !--- get seasonally varying degree day factor
    ddf_ice_day = (( min_ddf_ice(dyna_hru(ia)%ipar) + max_ddf_ice(dyna_hru(ia)%ipar) ) / 2. &
         + sin((julian_day(it) - 81) / 58.09) * &
         (max_ddf_ice(dyna_hru(ia)%ipar) - min_ddf_ice(dyna_hru(ia)%ipar) )/ 2.)/24*dt
   
    !--- apply a temperature lag to ice. Adjusted hru temp is already calculated
    !--- by snowmelt subroutine
    dyna_hru(ia)%glaciertemp = dyna_hru(ia)%glaciertemp * (1. - tlag_ice(dyna_hru(ia)%ipar) ) &
         + dyna_hru(ia)%temp * tlag_ice(dyna_hru(ia)%ipar)
    
    !--- if there is no snowpack above then the glacier hru can sublimate and melt
    if (dyna_hru(ia)%glacierdepth > 0.0 .and. dyna_hru(ia)%snowdepth == 0) then
       
       !--- reduce the sublim [0 1]. Allowing sublim at any temperature
         dyna_hru(ia)%ep = dyna_hru(ia)%ep * snow_ice_sublim(dyna_hru(ia)%ipar)
             
       !--- set sublim to pet. Note allowing this for all snow pack temperatures
         glacier_sublim = dyna_hru(ia)%ep 
         
         !--- no pet left it has all sublimated
         dyna_hru(ia)%ep = 0. 
             
       !--- potential melt
       if (dyna_hru(ia)%glaciertemp > melt_temp(dyna_hru(ia)%ipar)) then
           pot_icemelt = dyna_hru(ia)%glaciertemp * ddf_ice_day
       
       !---ice melt can't exceed amount of ice present
       glacier_melt = min(dyna_hru(ia)%glacierdepth, pot_icemelt)
       if (glacier_melt < 0.) glacier_melt = 0.
          
    end if
    
 end if
    
 !--- the accum ice is calculated from basal turnover coefficient. 
 !--- Beta is constant with time 
 if (dyna_hru(ia)%snowdepth > 0.0 ) then 
    glacier_accum = snow2ice(dyna_hru(ia)%ipar) * dyna_hru(ia)%snowdepth 
    dyna_hru(ia)%snowdepth = dyna_hru(ia)%snowdepth - glacier_accum 
 end if
 
 !--- keep record of glacier depth = accum - melting - sublim
 dyna_hru(ia)%glacierdepth  = dyna_hru(ia)%glacierdepth + glacier_accum - glacier_melt - glacier_sublim
    
   
 !--- make sure glacier depth never goes negative
 if (dyna_hru(ia)%glacierdepth < 0.)  dyna_hru(ia)%glacierdepth = 0.
 
 !--- change in glacier depth between timesteps
 dyna_hru(ia)%delta_glacierdepth =  dyna_hru(ia)%glacierdepth - dyna_hru(ia)%glacierdepth_old
 
 !--- diagnostic var.
 dyna_hru(ia)%smb = dyna_hru(ia)%delta_snowdepth + dyna_hru(ia)%delta_glacierdepth
 
 !--- Add glacier melt component to precip 
 dyna_hru(ia)%p = dyna_hru(ia)%p  + glacier_melt
   
end subroutine glacier

end module dyna_glacier


             
