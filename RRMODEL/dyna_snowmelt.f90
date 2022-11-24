module dyna_snowmelt
contains

!===============================================================
!  Calculate snowpack depth: sarah.shannon@bristol.ac.uk 2021
!  Lapse rate adjust air temperature and precipitation
!  Scale rain and snow fields 
!  Accumulate snowpack by converting rain to snow. Melt snow 
!  using a degree day model. Subliminate using PET
!  snowdepth has units m.w.eq
!===============================================================

      subroutine snowmelt (nac, &
           dyna_hru,        &
           it,              &
           ia,              &
           dt,              &
           min_ddf_snow,    &
           max_ddf_snow,    & 
           tlapse,          &
           plapse,          &
           rain2snow_temp,  &
           melt_temp,       &
           tlag_snow,       &
           snow_ice_sublim, &
           rcf,             &
           scf,             &
           nstep,           &
           julian_day)
        
          use dyna_common_types
          implicit none

          integer :: nac
          type(dyna_hru_type)  :: dyna_hru(nac)
          integer :: it
          doubleprecision :: dt
          integer :: nstep
          integer :: julian_day(nstep)
          integer :: ia
          
          !--- local snow variables ----------------
          double precision     :: accum
          double precision     :: pot_snowmelt
          double precision     :: snow_melt
          double precision     :: sublim
          double precision     :: ddf_snow_day
          double precision     :: delta_elev
          double precision, dimension(:) :: min_ddf_snow ! m w e oC-1 day-1
          double precision, dimension(:) :: max_ddf_snow ! m w e oC-1 day-1
          double precision, dimension(:) :: tlapse       ! oC per km
          double precision, dimension(:) :: plapse       ! % per 100 m (e.g. 0.5 is 50% per 100m)
          double precision, dimension(:) :: rain2snow_temp !oC
          double precision, dimension(:) :: melt_temp    !  oC
          double precision, dimension(:) :: tlag_snow    ! temp lag (0 - 1 ) 1 is no lag i.e.snow temp = air temp]) 
          double precision, dimension(:) :: snow_ice_sublim ! fraction to limit pet as sublim (1 is all sublim equal pet)
          double precision, dimension(:) :: scf          ! snowfall correct factor (range ~ 1-2) 
          double precision, dimension(:) :: rcf          ! rainfall correct factor (range ~ 1-2)            
          
          !--- Initialise local variables every timestep
          pot_snowmelt = 0.0 
          accum        = 0.0   
          snow_melt    = 0.0
          ddf_snow_day = 0.0
          sublim       = 0.0
          delta_elev   = 0.0
          
          !--- keep track of previous time step snowdepth
          dyna_hru(ia)%snowdepth_old = dyna_hru(ia)%snowdepth
                    
          !--- get seasonally varying degree day factor
          ddf_snow_day = (( min_ddf_snow(dyna_hru(ia)%ipar) + max_ddf_snow(dyna_hru(ia)%ipar) ) / 2. &
               + sin((julian_day(it) - 81) / 58.09) * &
               (max_ddf_snow(dyna_hru(ia)%ipar) - min_ddf_snow(dyna_hru(ia)%ipar) )/ 2.) / 24*dt
          
          !--- hru temp adjusted for elevation 
          dyna_hru(ia)%temp = dyna_hru(ia)%temp + (dyna_hru(ia)%elev &
               - dyna_hru(ia)%elev_climate) * tlapse(dyna_hru(ia)%ipar)/1000.0
                   
          !-- apply a temp lag to snowpack 
          dyna_hru(ia)%snowtemp = dyna_hru(ia)%snowtemp * (1. - tlag_snow(dyna_hru(ia)%ipar) ) &
               + dyna_hru(ia)%temp * tlag_snow(dyna_hru(ia)%ipar)
       
          !--- difference between climate elevation and hru elevation 
          delta_elev = dyna_hru(ia)%elev - dyna_hru(ia)%elev_climate
           
          if (dyna_hru(ia)%p > 0.0) then
            !--- snowing, no rain  
            if  (dyna_hru(ia)%temp <= rain2snow_temp(dyna_hru(ia)%ipar)) then
                accum = scf(dyna_hru(ia)%ipar) *dyna_hru(ia)%p &
                + scf(dyna_hru(ia)%ipar) *dyna_hru(ia)%p*delta_elev * (plapse(dyna_hru(ia)%ipar)/100) 
            !--- adjustments for elev below gridbox mean may result in negative accum/precip 
            !--- for large precip lapse rates. Set equal to zero.              
            if (accum < 0.) accum = 0.
                dyna_hru(ia)%p = 0.0
             
            else
            !---- raining, no snow  
                dyna_hru(ia)%p = rcf(dyna_hru(ia)%ipar) *dyna_hru(ia)%p &
                + rcf(dyna_hru(ia)%ipar)*dyna_hru(ia)%p*delta_elev * (plapse(dyna_hru(ia)%ipar)/100)    
              if (dyna_hru(ia)%p < 0.0) dyna_hru(ia)%p = 0.0
                accum = 0.0   
            end if 
          end if
          
          !--- Only melt and sublimate if there is a snowpack
          if (dyna_hru(ia)%snowdepth > 0.0 ) then
          
             !--- reduce the sublim [0 1]
             dyna_hru(ia)%ep = dyna_hru(ia)%ep * snow_ice_sublim(dyna_hru(ia)%ipar)
             
             !--- set sublim to pet. Note allowing this for all snow pack temperatures
             sublim = dyna_hru(ia)%ep 
             
             !--- no pet left it has all sublimated
             dyna_hru(ia)%ep = 0. 
             
             !--- Calc potential snow melt using degree day factor.
             !--- Melting can only happen if accum is zero i.e its not snowing
             !--- Only works for melt temp = 0
             if (dyna_hru(ia)%snowtemp > melt_temp(dyna_hru(ia)%ipar) .and. accum == 0.) then
                pot_snowmelt = dyna_hru(ia)%snowtemp * ddf_snow_day
             
                !---snow melt can't exceed amount of snow present 
                snow_melt = min(dyna_hru(ia)%snowdepth, pot_snowmelt)
                if (snow_melt < 0.) snow_melt = 0.
                    
             end if
          end if

          !--- keep record of snow depth every time step: depth = accum - melting - sublim
          dyna_hru(ia)%snowdepth  = dyna_hru(ia)%snowdepth + accum - snow_melt - sublim
                       
          !--- make sure snow depth never goes negative. It should not but just to be sure
          if (dyna_hru(ia)%snowdepth < 0.)  dyna_hru(ia)%snowdepth = 0.
                     
          !--- Add snow melt water equivalent to precip 
          dyna_hru(ia)%p = dyna_hru(ia)%p  + snow_melt
          
          !--- change in snowdepth between timesteps
          dyna_hru(ia)%delta_snowdepth = dyna_hru(ia)%snowdepth - dyna_hru(ia)%snowdepth_old
          
          !--- Output 2000 - 2007 HRU snowdepth for MODIS validation ----
    !if (it>=21548) then 
    !    write(43) sngl(dyna_hru(ia)%snowdepth)
    !end if 
              
        end subroutine snowmelt

end module dyna_snowmelt


             
