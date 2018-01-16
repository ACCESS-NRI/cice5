!============================================================================
!
module cpl_arrays_setup
!
! UM-AusCOM coupling involves following fields at the air-ice-sea interface:
!
! A> atm (UM) ==> ice (CICE) [* all at T cell center *]
!                                       
! (1) heatflux + solar radiation (total heatflux)	um_thflx
! (2) penetrating solar flux				um_pswflx
! (3) runoff						um_runoff
! (4) WME (?, for use in KT scheme in nemo)		um_wme
! (5) rainfall						um_rain
! (6) snowfall						um_snow
! (7) evaporation					um_evap
! (8) latent heat flux/evaporation			um_lhflx
! ( 9 - 13) top ice melting				um_tmlt(,,1:5)   
! (14 - 18) bottom ice melting 				um_bmlt(,,1:5)
! (19) windstress 'zonal'				um_taux
! (20) windstress 'meridional'				um_tauy
!     *** the above 20 fields are for the Hadgem3 coupling purpose  ***
!     *** some may not be necessary (such as WME and evap or lhflx) ***
!     *** for ACCESS we may need a few more fields: (temporary...)  ***
! (21) solar radiation (net down)			um_swflx
! (22) long wave radiation (net down)			um_lwflx
! (23) sensible heat flux				um_shflx
! (24) surface pressure					um_press
! (25) co2                                              um_co2
! (26) wind speed                                       um_wnd
! --------- add new fields for ACCESS2 --------
! (27) north (greenland) ice amount                     um_icenth
! (28) south (antarctic) ice amount                     um_icesth
! (29 - 33) ice surface/skin temperature                um_tsfice(,,1:5)
! (34 - 38) ice surface evaporation (sublimation)       um_iceevp(,,1:5)
!
! B> ocn (MOM4) ==> ice (CICE) [* at T or U cell center *]
!                          
! (1) sea surface temperature  (K)              	ocn_sst 
! (2) sea surface salinity   (psu)             	 	ocn_sss
! (3) zonal water speed      (m/s)              	ocn_ssu
! (4) meridional water speed (m/s)              	ocn_ssv
! (5) sea surface gradient (zonal)     (m/m)   	 	ocn_sslx
! (6) sea surface gradient (meridional)(m/m)    	ocn_ssly
! (7) potential ice frm/mlt heatflux (W/m^2)    	ocn_pfmice
! (8) co2 ()                                            ocn_co2
! (9) co2 flux ()                                       ocn_co2fx
!
! C> ice (CICE) ==> atm (UM) [* all from T to T, U or V cell center *]
!                          
! ( 1) ocean surface temperature  (K)       	        ia_sst
! ( 2 - 6 ) ice concentration (fraction)		ia_aicen(,,1:5)
! ( 7 - 11) snow thickness	(m ?)			ia_snown(,,1:5)
! (12 - 16) ice thickness	(m ?)			ia_thikn(,,1:5)
! (17) ice/ocn velocity 'zonal'				ia_uvel
! (18) ice/ocn velocity 'meridional'			ia_vvel
! (19) co2                                              ia_co2
! (20) co2 flux                                         ia_co2fx
! --------- add new fields for ACCESS2 --------
! (21) ocean surface freezing temperature               ia_sstfz
! (22 - 26 ) first order ice concentration              ia_foifr(,,1:5)
! (27 - 31 ) ice top layer temperature                  ia_itopt(,,1:5)
! (32 - 36 ) ice top layer effective conductivity       ia_itopk(,,1:5)
! (37 - 41 ) ice melt pond concentration                ia_pndfn(,,1:5)
! (42 - 46 ) ice melt pond thickness                    ia_pndtn(,,1:5)
!
! D> ice (CICE) ==> ocn (MOM4) [* at T or U cell center *]
!             
! (1) air/ice-ocean stress, x-direction (kg/m s^2)    	io_strsu
! (2) air/ice-ocean stress, y-direction (kg/m s^2)    	io_strsv
! (3) rainfall to ocean (kg/m^2/s)        	        io_rain
! (4) snowfall to ocean (kg/m^2/s)                      io_snow
! (5) salt flux to ocean (kg/m^2/s)               	io_stflx
! (6) 'net' heat flux to ocean (W/m^2)            	io_htflx
!     *(note word 'net' is misleading!) it is actually ice
!     *'melt' heatflux into ocean. (ref: ice_coupling.F, 
!     *it says:
!     *'buffs(n,index_i2c_Fioi_melth) = fhnet(i,j) 
!     *                          ! hf from melting'
! (7) shortwave penetrating to ocean (W/m^2)      	io_swflx
!     *** Also, we pass the following 'atmospheric fluxes': ***
! (8) latent heat flux/evaporation                	io_qflux
! (9) sensible heat flux                          	io_shflx
!(10) long wave radiation                         	io_lwflx
!(11) runoff (kg/m^2/s)                                	io_runof
!(12) pressure                                    	io_press
!(13) ice concentration (fraction)                      io_aice
!
! Seperate ice melting/formation associated water fluxes from the rainfall field:
!
!(14) ice melt waterflux                                io_melt
!(15) ice form waterflux                                io_form
!(16) co2                                               io_co2
!(17) wind speed                                        io_wnd
!
! 2 more added for "iceberg melt" (induced from land ice change):
!
!(18) iceberg melt waterflux                            io_licefw
!(19) iceberg melt heatflux                             io_liceht 
!
! Therefore, currently we have 
! 
! *for ACCESS1.x, 31 in, 33 out => thus jpfldout=33, jpfldin=31 in cpl_parameters.
! for ACCESS-CM2, 47 in, 63 out => thus jpfldout=63, jpfldin=47 in cpl_parameters. 
! now (20171024)  47 in, 65 out                  65          47
!----------------------------------------------------------------------------------
! This module will be largely modified/'simplifed' after ACCESS works !
!============================================================================

!cice stuff
use ice_kinds_mod

implicit none

! Fields in:
!===========
real(kind=dbl_kind), dimension(:,:,:), allocatable :: &   !from atm (UM)
    um_thflx, um_pswflx, um_runoff, um_wme, um_snow, um_rain, &
    um_evap,  um_lhflx,  um_taux,   um_tauy, &
    um_swflx, um_lwflx,  um_shflx,  um_press,um_co2, um_wnd, &
    um_icenth, um_icesth, &
    !!20171024 added for calculation of land ice increment
    lice_nth, lice_sth, msk_nth, msk_sth, amsk_nth, amsk_sth

real(kind=dbl_kind), dimension(:,:,:,:), allocatable :: &   
    um_tmlt, um_bmlt, um_tsfice, um_iceevp

! CORE runoff remapped onto the AusCOM grid (optional)
real(kind=dbl_kind), dimension(:,:,:), allocatable :: & 
    core_runoff

real(kind=dbl_kind), dimension(:,:,:), allocatable :: &   !from ocn (MOM4)
    ocn_sst, ocn_sss, ocn_ssu, ocn_ssv, ocn_sslx, ocn_ssly, ocn_pfmice, &
    ocn_co2, ocn_co2fx  

real(kind=dbl_kind), dimension(:,:), allocatable :: gwork
    !global domain work array, 4 coupling data passing and global data output. 
real(kind=dbl_kind), dimension(:,:,:), allocatable :: vwork  
    !local domain work array.

! Fields out:
!============
real(kind=dbl_kind),dimension(:,:,:), allocatable :: &     !to atm (timeaveraged)
    ia_sst, ia_uvel, ia_vvel, ia_co2, ia_co2fx, ia_sstfz
real(kind=dbl_kind), dimension(:,:,:,:), allocatable :: &
    ia_aicen, ia_snown, ia_thikn, &
    ia_foifr, ia_itopt, ia_itopk, ia_pndfn, ia_pndtn

real(kind=dbl_kind),dimension(:,:,:), allocatable :: &     !to ocn (time averaged)
    io_strsu, io_strsv, io_rain,  io_snow,  io_stflx, io_htflx, io_swflx, &
    io_qflux, io_shflx, io_lwflx, io_runof, io_press, io_aice, &
    io_melt, io_form, io_co2, io_wnd, &
    io_licefw, io_liceht         !2 more added 20171024

! Temporary arrays
!==================

! 1. ice fields averaged over IA cpl interval:
real(kind=dbl_kind),dimension(:,:,:), allocatable :: &
    maiu, muvel, mvvel
real(kind=dbl_kind), dimension(:,:,:,:), allocatable :: &
    maicen, msnown, mthikn, &
    mfoifr, mitopt, mitopk, mpndfn, mpndtn
real(kind=dbl_kind), dimension(:,:,:,:), allocatable :: &       !BX: just in case......
    maicen_saved

! 2. ice fields averaged over IO cpl interval:
real(kind=dbl_kind),dimension(:,:,:), allocatable :: &     
    maice, mstrocnxT, mstrocnyT, mfresh, mfsalt, mfhocn, mfswthru, msicemass
!BX---extra one: (maybe better save the IA interval one......anyway:)
!real(kind=dbl_kind),dimension(:,:,:), allocatable :: &
!    maice_saved

! 3. ocn fields averaged over IA cpl interval:
real(kind=dbl_kind),dimension(:,:,:), allocatable :: &
    msst, mssu, mssv, mco2, mco2fx, msstfz


! other stuff 
!============
real(kind=dbl_kind),dimension(:,:,:), allocatable :: & 
    sicemass   !ice mass

real(kind=dbl_kind),dimension(:,:,:,:), allocatable :: &
    icebergfw   !land ice discharge into ocean as monthly iceberg melt waterflux ==>io_licefw

!===========================================================================
end module cpl_arrays_setup

