MODULE cpl_forcing_handler
!
! It contains subroutines handling coupling fields. 
!
use ice_blocks
use ice_forcing
use ice_read_write
use ice_domain_size
use ice_domain,    only : distrb_info, nblocks 
use ice_flux            !forcing data definition (Tair, Qa, uocn, etc.)
                        !Tn_top, keffn_top ...(for multilayer configuration)   
!use ice_state,     only : aice, aicen, trcr, trcrn, nt_hpnd   !ice concentration and tracers
use ice_state,     only : aice, aicen, trcr, trcrn, nt_hpnd, nt_Tsfc   !ice concentration and tracers
use ice_state,     only: uvel, vvel, vsnon, vicen
use ice_gather_scatter
use ice_constants
use ice_grid,      only : tmask, to_ugrid
use ice_communicate, only : my_task, master_task
!use ice_ocean,     only : cprho
use ice_exit,      only : abort_ice
use ice_shortwave, only : apeffn

use cpl_parameters
use cpl_netcdf_setup
use cpl_arrays_setup

implicit none

real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
  aiiu       ! ice fraction on u-grid 

contains

!=================================================
subroutine get_core_runoff(fname, vname, nrec)
! read in the remapped core runoff data (S.Marsland) which will be used to replace
! the ncep2 runoff sent from matm via coupler 

implicit none

character*(*), intent(in) :: fname, vname
integer(kind=int_kind), intent(in) :: nrec
logical :: dbug
integer(kind=int_kind) :: ncid

dbug = .true.

if ( file_exist(fname) ) then
  call ice_open_nc(fname, ncid)
  call ice_read_global_nc(ncid, nrec, vname, gwork, dbug)
  call scatter_global(core_runoff, gwork, master_task, distrb_info, &
                      field_loc_center, field_type_scalar)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task == 0) write(il_out,*) '(get_core_runoff) file doesnt exist: ', fname
  stop 'CICE stopped: core runoff (remapped) file not found.'
endif

return
end subroutine get_core_runoff

!=================================================
subroutine get_time0_sstsss(fname, nmonth)

! This routine is to be used only once at the beginning at an exp.

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nmonth
logical :: dbug
integer(kind=int_kind) :: ncid

dbug = .true.
!dbug = .false.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) opening ncfile: ',fname
  endif
  call ice_open_nc(fname, ncid)
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) reading in initial SST...'
  endif
  call ice_read_nc(ncid, nmonth, 'TEMP', sst, dbug)
  call gather_global(gwork, sst, master_task, distrb_info)
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) reading in initial SSS...' 
  endif
  call ice_read_nc(ncid, nmonth, 'SALT', sss, dbug)
  call gather_global(gwork, sss, master_task, distrb_info)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) file doesnt exist: ', fname
  endif
  call abort_ice('CICE stopped--initial SST and SSS ncfile not found.')  
endif

return
end subroutine get_time0_sstsss

!=================================================
! temporary use ...
subroutine read_access_a2i_data(fname,nrec,istep) 

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nrec,istep
logical :: dbug
integer(kind=int_kind) :: ncid

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(read_access_a2i_data) opening ncfile: ',fname
  endif
  call ice_open_nc(fname, ncid)
  if (my_task==0) then
    write(il_out,*) '(read_access_a2i_data) reading a2i forcing data...'
  endif
  call ice_read_nc(ncid, nrec, 'thflx_i', um_thflx, dbug)
  call ice_read_nc(ncid, nrec, 'pswflx_i', um_pswflx, dbug)
  call ice_read_nc(ncid, nrec, 'runoff_i', um_runoff, dbug)
  call ice_read_nc(ncid, nrec, 'wme_i', um_wme, dbug)
  call ice_read_nc(ncid, nrec, 'rain_i', um_rain, dbug)
  call ice_read_nc(ncid, nrec, 'snow_i', um_snow, dbug)
  call ice_read_nc(ncid, nrec, 'evap_i', um_evap, dbug)
  call ice_read_nc(ncid, nrec, 'lhflx_i', um_lhflx, dbug)
  call ice_read_nc(ncid, nrec, 'tmlt01_i', um_tmlt(:,:,1,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt02_i', um_tmlt(:,:,2,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt03_i', um_tmlt(:,:,3,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt04_i', um_tmlt(:,:,4,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt05_i', um_tmlt(:,:,5,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt01_i', um_bmlt(:,:,1,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt02_i', um_bmlt(:,:,2,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt03_i', um_bmlt(:,:,3,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt04_i', um_bmlt(:,:,4,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt05_i', um_bmlt(:,:,5,:), dbug)
  call ice_read_nc(ncid, nrec, 'taux_i', um_taux, dbug)
  call ice_read_nc(ncid, nrec, 'tauy_i', um_tauy, dbug)
  call ice_read_nc(ncid, nrec, 'swflx_i', um_swflx,  dbug)
  call ice_read_nc(ncid, nrec, 'lwflx_i', um_lwflx,  dbug)
  call ice_read_nc(ncid, nrec, 'shflx_i', um_shflx,  dbug)
  call ice_read_nc(ncid, nrec, 'press_i', um_press,  dbug)
  call ice_read_nc(ncid, nrec, 'co2_ai', um_co2,  dbug)
  call ice_read_nc(ncid, nrec, 'wnd_ai', um_wnd,  dbug)
!
  call ice_read_nc(ncid, nrec, 'icenth_i', um_icenth, dbug)
  call ice_read_nc(ncid, nrec, 'icesth_i', um_icesth, dbug)
  call ice_read_nc(ncid, nrec, 'tsfice01', um_tsfice(:,:,1,:), dbug)
  call ice_read_nc(ncid, nrec, 'tsfice02', um_tsfice(:,:,2,:), dbug)
  call ice_read_nc(ncid, nrec, 'tsfice03', um_tsfice(:,:,3,:), dbug)
  call ice_read_nc(ncid, nrec, 'tsfice04', um_tsfice(:,:,4,:), dbug)
  call ice_read_nc(ncid, nrec, 'tsfice05', um_tsfice(:,:,5,:), dbug)
  call ice_read_nc(ncid, nrec, 'iceevp01', um_iceevp(:,:,1,:), dbug)
  call ice_read_nc(ncid, nrec, 'iceevp02', um_iceevp(:,:,2,:), dbug)
  call ice_read_nc(ncid, nrec, 'iceevp03', um_iceevp(:,:,3,:), dbug)
  call ice_read_nc(ncid, nrec, 'iceevp04', um_iceevp(:,:,4,:), dbug)
  call ice_read_nc(ncid, nrec, 'iceevp05', um_iceevp(:,:,5,:), dbug)

  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task==0) then
    write(il_out,*) '(ed_access_a2i_data file doesnt exist: ', fname
  endif
  call abort_ice('CICE stopped--ACCESS fields_a2i ncfile not found.')
endif

call check_a2i_fields(istep)

end subroutine read_access_a2i_data

!=================================================
subroutine atm_icefluxes_back2GBM
!convert the a2i fluxes into GBM units for those that are scaled up in UM 
!by "/maicen" before being sent to cice [needed for GSI8 TTI approach].
 
implicit none

um_tmlt(:,:,:,:) = um_tmlt(:,:,:,:) * maicen_saved(:,:,:,:)
um_bmlt(:,:,:,:) = um_bmlt(:,:,:,:) * maicen_saved(:,:,:,:)
um_iceevp(:,:,:,:) = um_iceevp(:,:,:,:) * maicen_saved(:,:,:,:)

end subroutine atm_icefluxes_back2GBM

!=================================================
subroutine read_restart_i2a(fname, sec) !'i2a.nc', 0)

! read ice to atm coupling fields from restart file, and send to atm module

implicit none
character*(*), intent(in) :: fname
integer :: sec

integer(kind=int_kind) :: ncid
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(read_restart_i2a) reading in i2a fields......'
  endif
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, 1, 'icecon01',   ia_aicen(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon02',   ia_aicen(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon03',   ia_aicen(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon04',   ia_aicen(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon05',   ia_aicen(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk01',   ia_snown(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk02',   ia_snown(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk03',   ia_snown(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk04',   ia_snown(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk05',   ia_snown(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk01',   ia_thikn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk02',   ia_thikn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk03',   ia_thikn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk04',   ia_thikn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk05',   ia_thikn(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'isst_ia',    ia_sst,   dbug)
  call ice_read_nc(ncid, 1, 'uvel_ia',    ia_uvel,  dbug)
  call ice_read_nc(ncid, 1, 'vvel_ia',    ia_vvel,  dbug)
  call ice_read_nc(ncid, 1, 'co2_i2',     ia_co2,   dbug)
  call ice_read_nc(ncid, 1, 'co2fx_i2',   ia_co2fx, dbug)
  call ice_read_nc(ncid, 1, 'sstfz_ia',   ia_sstfz, dbug)
  call ice_read_nc(ncid, 1, 'foifr01',    ia_foifr(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'foifr02',    ia_foifr(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'foifr03',    ia_foifr(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'foifr04',    ia_foifr(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'foifr05',    ia_foifr(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopt01',    ia_itopt(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopt02',    ia_itopt(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopt03',    ia_itopt(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopt04',    ia_itopt(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopt05',    ia_itopt(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopk01',    ia_itopk(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopk02',    ia_itopk(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopk03',    ia_itopk(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopk04',    ia_itopk(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'itopk05',    ia_itopk(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndfn01',    ia_pndfn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndfn02',    ia_pndfn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndfn03',    ia_pndfn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndfn04',    ia_pndfn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndfn05',    ia_pndfn(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndtn01',    ia_pndtn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndtn02',    ia_pndtn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndtn03',    ia_pndtn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndtn04',    ia_pndtn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'pndtn05',    ia_pndtn(:,:,5,:),   dbug)

  if (my_task == master_task) then
    call ice_close_nc(ncid)
    write(il_out,*) '(read_restart_i2a) has read in 18 i2a fields.'
  endif

else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (read_restart_i2a) not found file *** ',fname
  endif
  print *, 'CICE: (read_restart_i2a) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 i2a data file.')
endif
end subroutine read_restart_i2a

!=================================================
subroutine read_restart_i2asum(fname, sec) !'i2a.nc', 0)

! read ice to atm coupling fields from restart file, and send to atm module

implicit none
character*(*), intent(in) :: fname
integer :: sec

integer(kind=int_kind) :: ncid
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(read_restart_i2asum) reading in i2a fields......'
  endif
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, 1, 'maicen1',   maicen(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen2',   maicen(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen3',   maicen(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen4',   maicen(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen5',   maicen(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown1',   msnown(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown2',   msnown(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown3',   msnown(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown4',   msnown(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown5',   msnown(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn1',   mthikn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn2',   mthikn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn3',   mthikn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn4',   mthikn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn5',   mthikn(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'msst',      msst,   dbug)
  call ice_read_nc(ncid, 1, 'mssu',      mssu,   dbug)
  call ice_read_nc(ncid, 1, 'mssv',      mssv,   dbug)
  call ice_read_nc(ncid, 1, 'muvel',     muvel,  dbug)
  call ice_read_nc(ncid, 1, 'mvvel',     mvvel,  dbug)
  call ice_read_nc(ncid, 1, 'maiu',      maiu,   dbug)
  !
  !call ice_read_nc(ncid, 1, 'maice_ia', maice_ia,          dbug)
  call ice_read_nc(ncid, 1, 'mfoifr01', mfoifr(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mfoifr02', mfoifr(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mfoifr03', mfoifr(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mfoifr04', mfoifr(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mfoifr05', mfoifr(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopt01', mitopt(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopt02', mitopt(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopt03', mitopt(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopt04', mitopt(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopt05', mitopt(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopk01', mitopk(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopk02', mitopk(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopk03', mitopk(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopk04', mitopk(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mitopk05', mitopk(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndfn01', mpndfn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndfn02', mpndfn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndfn03', mpndfn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndfn04', mpndfn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndfn05', mpndfn(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndtn01', mpndtn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndtn02', mpndtn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndtn03', mpndtn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndtn04', mpndtn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mpndtn05', mpndtn(:,:,5,:),   dbug)

  if (my_task == master_task) then
    call ice_close_nc(ncid)
    write(il_out,*) '(read_restart_i2asum) has read in 21 i2a fields.'
  endif

else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (read_restart_i2asum) not found file *** ',fname
  endif
  print *, 'CICE: (read_restart_i2asum) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 i2a data file.')
endif
end subroutine read_restart_i2asum

!=================================================
subroutine put_restart_i2a(fname, sec)
! call this subroutine after called get_restart_oi2
! it uses ocn_sst etc to calculate average ocn fields which will be used to send 
!to atm at the 1st step of continue run, because the ocn_sst cannot be sent to ice at the end of last run. 
! average ice fields (done at end of last run) are ready by calling read_restart_i2asum() 
!
implicit none

character*(*), intent(in) :: fname
integer :: sec

  if ( file_exist('i2a.nc') ) then
    write(il_out,*)' calling read_restart_i2a at time_sec = ',sec
    call read_restart_i2a('i2a.nc', sec)  
  endif
  if ( file_exist('i2asum.nc') ) then
    write(il_out,*)' calling read_restart_i2asum at time_sec = ',sec
    call read_restart_i2asum('i2asum.nc', sec) 

    write(il_out,*)' calling ave_ocn_fields_4_i2a at time_sec = ',sec
    call time_average_ocn_fields_4_i2a  !accumulate/average ocn fields needed for IA coupling
    write(il_out,*) ' calling get_i2a_fields at time_sec =', sec
    call get_i2a_fields
  endif

end subroutine put_restart_i2a

!=================================================
subroutine get_restart_o2i(fname)

! To be called at beginning of each run trunk to read in restart o2i fields

implicit none

character*(*), intent(in) :: fname
 
integer(kind=int_kind) :: ncid_o2i
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_restart_o2i) reading in o2i fields......'
  endif
  call ice_open_nc(fname, ncid_o2i)
  call ice_read_nc(ncid_o2i, 1, 'sst_i',    ocn_sst,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'sss_i',    ocn_sss,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssu_i',    ocn_ssu,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssv_i',    ocn_ssv,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'sslx_i',   ocn_sslx,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssly_i',   ocn_ssly,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'pfmice_i', ocn_pfmice, dbug)
  call ice_read_nc(ncid_o2i, 1, 'co2_oi',   ocn_co2, dbug)
  call ice_read_nc(ncid_o2i, 1, 'co2fx_oi',  ocn_co2fx, dbug)
  if (my_task == master_task) then
    call ice_close_nc(ncid_o2i)
    write(il_out,*) '(get_restart_o2i) has read in 7 o2i fields.'
  endif
else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (get_restart_o2i) not found file *** ',fname
  endif
  print *, 'CICE: (get_restart_o2i) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 o2i data file.')
endif

return
end subroutine get_restart_o2i

!=================================================
subroutine get_restart_mice(fname)

! Called at beginning of the run to get 'last' IO cpl int T-M ice variables 
! which are used together with the first received i2a fields to obtain the first 
! i2o fields sent to ocn immediately as the 1st io cpl int forcing there.

implicit none

character*(*), intent(in) :: fname

integer(kind=int_kind) :: ncid_o2i
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_restart_mice) reading in mice variables......'
  endif
  call ice_open_nc(fname, ncid_o2i)
!B: 20170825 ==> need maicen_saved variables
  call ice_read_nc(ncid_o2i, 1, 'maicen1',   maicen_saved(:,:,1,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen2',   maicen_saved(:,:,2,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen3',   maicen_saved(:,:,3,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen4',   maicen_saved(:,:,4,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen5',   maicen_saved(:,:,5,:), dbug)
!b.
  call ice_read_nc(ncid_o2i, 1, 'maice',     maice,     dbug)
  call ice_read_nc(ncid_o2i, 1, 'mstrocnxT', mstrocnxT, dbug)
  call ice_read_nc(ncid_o2i, 1, 'mstrocnyT', mstrocnyT, dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfresh',    mfresh,    dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfsalt',    mfsalt,    dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfhocn',    mfhocn,    dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfswthru',  mfswthru,  dbug)
  call ice_read_nc(ncid_o2i, 1, 'msicemass', msicemass, dbug)

  if (my_task == master_task) then
    call ice_close_nc(ncid_o2i)
    write(il_out,*) '(get_restart_mice) has read in 8 T-M variables.'
  endif
else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (get_restart_mice) not found file *** ',fname
  endif
  print *, 'CICE: (get_restart_mice) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 mice data file.')
endif

return
end subroutine get_restart_mice

!=================================================
subroutine get_restart_i2o(fname)

! To be called at beginning of each run trunk to read in restart i2o fields

implicit none

character*(*), intent(in) :: fname

integer(kind=int_kind) :: ncid_i2o, jf, jfs
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_time0_i2o_fields) reading in i2o fields......'
  endif
  call ice_open_nc(fname, ncid_i2o)
  do jf = nsend_i2a + 1, jpfldout
    call ice_read_nc(ncid_i2o, 1, cl_writ(jf) , vwork, dbug)
    select case(trim(cl_writ(jf)))
    case ('strsu_io'); io_strsu = vwork
    case ('strsv_io'); io_strsv = vwork
    case ('rain_io');  io_rain = vwork
    case ('snow_io');  io_snow = vwork
    case ('stflx_io'); io_stflx = vwork
    case ('htflx_io'); io_htflx = vwork
    case ('swflx_io'); io_swflx = vwork
    case ('qflux_io'); io_qflux = vwork
    case ('shflx_io'); io_shflx = vwork
    case ('lwflx_io'); io_lwflx = vwork
    case ('runof_io'); io_runof = vwork
    case ('press_io'); io_press = vwork
    case ('aice_io');  io_aice  = vwork
    case ('melt_io');  io_melt  = vwork
    case ('form_io');  io_form  = vwork
    case ('co2_i1');  io_co2  = vwork
    case ('wnd_i1');  io_wnd  = vwork
    end select
  enddo
  if (my_task == master_task) then
    call ice_close_nc(ncid_i2o)
    write(il_out,*) '(get_time0_i2o_fields) has read in 11 i2o fields.'
  endif
else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (get_time0_i2o_fields) not found file *** ',fname
  endif
  print *, 'CICE: (get_time0_i2o_fields_old) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 i2o data file.')
endif

return
end subroutine get_restart_i2o

!=================================================
subroutine set_sbc_ice          !!NOTE: This routine is NOT used!! 
!
! Set coupling fields (in units of GMB, from UM and MOM4) needed for CICE
!
! Adapted from "subroutine cice_sbc_in" of HadGem3 Nemo "MODULE sbcice_cice"
! for the "nsbc = 5" case.
!
! It should be called after calling "from_atm" and "from_ocn". 
!-------------------------------------------------------------------------------

implicit none

real :: r1_S0
real, dimension(nx_block,ny_block,nblocks) :: zzs

integer :: i,j,k,cat

!*** Fields from UM (all on T cell center):

!(1) windstress taux:
strax = um_taux * maice      !*tmask   

!(2) windstress tauy:
stray = um_tauy * maice      !*tmask

!(3) surface downward latent heat flux (==> multi-category)
do j = 1, ny_block
do i = 1, nx_block
  do k = 1, nblocks
    if (maice(i,j,k)==0.0) then
      do cat = 1, ncat
        flatn_f(i,j,cat,k) = 0.0
      enddo
      ! This will then be conserved in CICE (done in sfcflux_to_ocn)
      flatn_f(i,j,1,k) = um_lhflx(i,j,k) 
    else
      do cat = 1, ncat
        !!!B: flatn_f(i,j,cat,k) = um_lhflx(i,j,k) * maicen(i,j,cat,k)/maice(i,j,k)
        !???: flatn_f(i,j,cat,k) = um_iceevp(i,j,cat,k) * Lsub 
        flatn_f(i,j,cat,k) = - um_iceevp(i,j,cat,k) * Lsub
      enddo
    endif
  enddo
enddo
enddo 

! GBM conductive flux through ice:
!(4-8) top melting; (9-13) bottom belting ==> surface heatflux
do cat = 1, ncat
  fcondtopn_f(:,:,cat,:) = um_bmlt(:,:,cat,:)
  fsurfn_f   (:,:,cat,:) = um_tmlt(:,:,cat,:) + um_bmlt(:,:,cat,:)
enddo

!(14) snowfall
fsnow = max(maice * um_snow, 0.0)

!(15) rainfall
frain = max(maice * um_rain, 0.0)

!BX-20160718: "save" the ice concentration "maice" used here for scaling-up frain etc in
!ice_step for "consistency"--
!maice_saved = maice 

!*** Fields from MOM4 (SSU/V and sslx/y are on U points): 

!(1) freezing/melting potential
frzmlt = ocn_pfmice
if (limit_icemelt) then
  frzmlt(:,:,:) = max(frzmlt(:,:,:), meltlimit)
endif

!(2) SST
!make sure SST is 'all right' K==>C
sst = ocn_sst
if (maxval(sst).gt.200) then
  sst = sst -273.15
endif

!(3) SSS
sss = ocn_sss

!(4) SSU
uocn = ocn_ssu

!(5) SSV
vocn = ocn_ssv

!(6) surface slope sslx
ss_tltx = ocn_sslx

!(7) surface slope ssly
ss_tlty = ocn_ssly

!(as per S.O.) make sure Tf if properly initialized
Tf (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
!
!B: May use different formula for Tf such as TEOS-10 formulation: 
!
!r1_S0 = 0.875/35.16504
!zzs(:,:,:) = sqrt(abs(sss(:,:,:)) * r1_S0)
!Tf(:,:,:) = ((((1.46873e-03 * zzs(:,:,:) - 9.64972e-03) * zzs(:,:,:) + &
!               2.28348e-02) * zzs(:,:,:) - 3.12775e-02) * zzs(:,:,:) + &
!               2.07679e-02) * zzs(:,:,:) - 5.87701e-02
!Tf(:,:,:) = Tf(:,:,:) * sss(:,:,:) ! - 7.53e-4 * 5.0 !!!5.0 is depth in meters

end subroutine set_sbc_ice

!===============================================================================
subroutine get_sbc_ice
!
! ** Purpose: set GBM coupling fields (from UM and MOM4) needed for CICE
!
! Adapted from "subroutine cice_sbc_in" of HadGem3 Nemo "MODULE sbcice_cice"
!    for the "nsbc = 5" case.
!
! It should be called after calling "from_atm" and "from_ocn". 
!-------------------------------------------------------------------------------

implicit none

real :: r1_S0
real, dimension(nx_block,ny_block,nblocks) :: zzs

integer :: i,j,k,cat

! Fields from UM (all on T cell center):

!(1) windstress taux:
strax = um_taux * aice      !*tmask ?   

!(2) windstress tauy:
stray = um_tauy * aice      !*tmask ?

!(3) surface downward latent heat flux (==> multi_category)
!BX: where is flatn_f "used" in CICE?
do j = 1, ny_block
do i = 1, nx_block
  do k = 1, nblocks
    !BX 20160826: as in NEMO sbccpl.F90, there is no "open water field" um_lhflx involved: 
    !    qla_ice(:,:,1:jpl) = - frcv(jpr_ievp)%z3(:,:,1:jpl) * lsub 
    !-------------------------------------------------------------------------------------
    !if (aice(i,j,k)==0.0) then
    !  do cat = 1, ncat
    !    flatn_f(i,j,cat,k) = 0.0
    !  enddo
    !  ! This will then be conserved in CICE (done in sfcflux_to_ocn)
    !  flatn_f(i,j,1,k) = um_lhflx(i,j,k)
    !else
      do cat = 1, ncat
        !!!BX: flatn_f(i,j,cat,k) = um_lhflx(i,j,k) * aicen(i,j,cat,k)/aice(i,j,k)
        !!!   Double check "Lsub" used here !!! 
        !?! flatn_f(i,j,cat,k) = um_iceevp(i,j,cat,k) * Lsub 
        flatn_f(i,j,cat,k) = - um_iceevp(i,j,cat,k) * Lsub
      enddo
    !endif
  enddo
enddo
enddo

! GBM conductive flux through ice:
!(4-8) top melting; (9-13) bottom belting ==> surface heatflux
do cat = 1, ncat
  fcondtopn_f(:,:,cat,:) = um_bmlt(:,:,cat,:)
  fsurfn_f   (:,:,cat,:) = um_tmlt(:,:,cat,:) + um_bmlt(:,:,cat,:)
enddo

!!! 20130419: Martin Dix's investigation suggests that frain and fsnow should NOT be scaled by 
!!!           aice here. This scaling would caused double-scaling with "fresh" calculation.. 
!(14) snowfall
!!!fsnow = max(aice * um_snow,0.0)
!fsnow = max(um_snow,0.0)            !no more scaling as per M.D.!
!(15) rainfall
!!!frain = max(aice * um_rain,0.0)
!frain = max(um_rain,0.0)            !no more scaling as per M.D.!
!!! 20130420: I dug deeper and checked all the associated steps of "fresh" calculation, found
!!!           the original weighting is CORRECT! so back to *aice:
fsnow = max(aice * um_snow,0.0)
frain = max(aice * um_rain,0.0)  
!
!ice surface skin temperature (from UM)-------------------------------------
!see: tsfc_ice definition in sbccpl.F90 at
!/short/p66/hxy599/fcm_make_ocean_GC3/extract/nemo/NEMOGCM/NEMO/OPA_SRC/SBC
!---------------------------------------------------------------------------
do cat = 1, ncat
  !!!  trcrn(:,:,nt_Tsfc,cat,:) = um_tsfice(:,:,cat,:)
  do k = 1, nblocks
    do j = 1, ny_block
      do i = 1, nx_block
        if (um_tsfice(i,j,cat,k) > 0.0) then
          trcrn(i,j,nt_Tsfc,cat,k) = 0.0 
        else if (um_tsfice(i,j,cat,k) < -60.0) then
          trcrn(i,j,nt_Tsfc,cat,k) = -60.0 
        else
          trcrn(i,j,nt_Tsfc,cat,k) = um_tsfice(i,j,cat,k)
        endif
      enddo
    enddo
  enddo
enddo

! Fields from MOM4 (SSU/V and sslx/y are on U points): 

!(1) freezing/melting potential
frzmlt = ocn_pfmice
!20080312: set maximum melting htflux allowed from ocn, (eg, -200 W/m^2)
!          the artificial "meltlimit = -200 " is read in from input_ice.nml
!20090320: set option 'limit_icemelt' in case no limit needed if cice behaves!
if (limit_icemelt) then
  frzmlt(:,:,:) = max(frzmlt(:,:,:), meltlimit)
endif

!(2) SST
sst = ocn_sst -273.15

!(3) SSS
sss = ocn_sss

!(4) SSU
uocn = ocn_ssu

!(5) SSV
vocn = ocn_ssv
!(6) surface slope sslx
                                                             
ss_tltx = ocn_sslx

!(7) surface slope ssly
ss_tlty = ocn_ssly

! * (as per S. O'Farrel) make sure Tf if properly initialized
!----- should use eos formula to calculate Tf for "consistency" with GCx ----!
!Tf (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
!
!B: May use different formula for Tf such as TEOS-10 formulation: 
!
!r1_S0 = 0.875/35.16504
!zzs(:,:,:) = sqrt(abs(sss(:,:,:)) * r1_S0)
!Tf(:,:,:) = ((((1.46873e-03 * zzs(:,:,:) - 9.64972e-03) * zzs(:,:,:) + &
!               2.28348e-02) * zzs(:,:,:) - 3.12775e-02) * zzs(:,:,:) + &
!               2.07679e-02) * zzs(:,:,:) - 5.87701e-02
!Tf(:,:,:) = Tf(:,:,:) * sss(:,:,:) ! - 7.53e-4 * 5.0 !!!5.0 is depth in meters
!
end subroutine get_sbc_ice

!=================================================
subroutine save_restart_o2i(fname, nstep)

! output the last o2i forcing data received in cice by the end of the run, 
! to be read in at the beginning of next run by cice

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

if (my_task == 0) then
  call create_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
  call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
endif

do jf = nrecv_a2i + 1, jpfldin

  select case (trim(cl_read(jf)))
  case('sst_i'); vwork = ocn_sst
  case('sss_i'); vwork = ocn_sss
  case('ssu_i'); vwork = ocn_ssu
  case('ssv_i'); vwork = ocn_ssv
  case('sslx_i'); vwork = ocn_sslx
  case('ssly_i'); vwork = ocn_ssly
  case('pfmice_i'); vwork = ocn_pfmice
  case('co2_oi'); vwork = ocn_co2
  case('co2fx_oi'); vwork = ocn_co2fx
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)
  if (my_task == 0) then
    call write_nc2D(ncid, cl_read(jf), gwork, 2, il_im, il_jm, 1, ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck( nf_close(ncid) )

return
end subroutine save_restart_o2i

!=================================================
subroutine save_restart_i2asum(fname, nstep)
! output the last i2a forcing data in cice at the end of the run,
! to be read in at the beginning of next run by cice and sent to atm

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

integer(kind=int_kind), parameter :: sumfldin  = 46     !21
character(len=8), dimension(sumfldin) :: sumfld

sumfld(1)='msst'
sumfld(2)='mssu'
sumfld(3)='mssv'
sumfld(4)='muvel'
sumfld(5)='mvvel'
sumfld(6)='maiu'
sumfld(7)='maicen1'
sumfld(8)='maicen2'
sumfld(9)='maicen3'
sumfld(10)='maicen4'
sumfld(11)='maicen5'
sumfld(12)='mthikn1'
sumfld(13)='mthikn2'
sumfld(14)='mthikn3'
sumfld(15)='mthikn4'
sumfld(16)='mthikn5'
sumfld(17)='msnown1'
sumfld(18)='msnown2'
sumfld(19)='msnown3'
sumfld(20)='msnown4'
sumfld(21)='msnown5'
!
sumfld(22)='mfoifr1'
sumfld(23)='mfoifr2'
sumfld(24)='mfoifr3'
sumfld(25)='mfoifr4'
sumfld(26)='mfoifr5'
sumfld(27)='mitopt1'
sumfld(28)='mitopt2'
sumfld(29)='mitopt3'
sumfld(30)='mitopt4'
sumfld(31)='mitopt5'
sumfld(32)='mitopk1'
sumfld(33)='mitopk2'
sumfld(34)='mitopk3'
sumfld(35)='mitopk4'
sumfld(36)='mitopk5'
sumfld(37)='mpndfn1'
sumfld(38)='mpndfn2'
sumfld(39)='mpndfn3'
sumfld(40)='mpndfn4'
sumfld(41)='mpndfn5'
sumfld(42)='mpndtn1'
sumfld(43)='mpndtn2'
sumfld(44)='mpndtn3'
sumfld(45)='mpndtn4'
sumfld(46)='mpndtn5'

if (my_task == 0) then
  call create_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
endif

do jf = 1, sumfldin
    select case (trim(sumfld(jf)))
    case('msst'); vwork = msst
    case('mssu'); vwork = mssu
    case('mssv'); vwork = mssv
    case('muvel'); vwork = muvel
    case('mvvel'); vwork = mvvel
    case('maiu'); vwork = maiu
    case('maicen1'); vwork = maicen(:,:,1,:)
    case('maicen2'); vwork = maicen(:,:,2,:)
    case('maicen3'); vwork = maicen(:,:,3,:)
    case('maicen4'); vwork = maicen(:,:,4,:)
    case('maicen5'); vwork = maicen(:,:,5,:)
    case('mthikn1'); vwork = mthikn(:,:,1,:)
    case('mthikn2'); vwork = mthikn(:,:,2,:)
    case('mthikn3'); vwork = mthikn(:,:,3,:)
    case('mthikn4'); vwork = mthikn(:,:,4,:)
    case('mthikn5'); vwork = mthikn(:,:,5,:)
    case('msnown1'); vwork = msnown(:,:,1,:)
    case('msnown2'); vwork = msnown(:,:,2,:)
    case('msnown3'); vwork = msnown(:,:,3,:)
    case('msnown4'); vwork = msnown(:,:,4,:)
    case('msnown5'); vwork = msnown(:,:,5,:)
    case('mfoifr1'); vwork = mfoifr(:,:,1,:)
    case('mfoifr2'); vwork = mfoifr(:,:,2,:)
    case('mfoifr3'); vwork = mfoifr(:,:,3,:)
    case('mfoifr4'); vwork = mfoifr(:,:,4,:)
    case('mfoifr5'); vwork = mfoifr(:,:,5,:)
    case('mitopt1'); vwork = mitopt(:,:,1,:)
    case('mitopt2'); vwork = mitopt(:,:,2,:)
    case('mitopt3'); vwork = mitopt(:,:,3,:)
    case('mitopt4'); vwork = mitopt(:,:,4,:)
    case('mitopt5'); vwork = mitopt(:,:,5,:)
    case('mitopk1'); vwork = mitopk(:,:,1,:)
    case('mitopk2'); vwork = mitopk(:,:,2,:)
    case('mitopk3'); vwork = mitopk(:,:,3,:)
    case('mitopk4'); vwork = mitopk(:,:,4,:)
    case('mitopk5'); vwork = mitopk(:,:,5,:)
    case('mpndfn1'); vwork = mpndfn(:,:,1,:)
    case('mpndfn2'); vwork = mpndfn(:,:,2,:)
    case('mpndfn3'); vwork = mpndfn(:,:,3,:)
    case('mpndfn4'); vwork = mpndfn(:,:,4,:)
    case('mpndfn5'); vwork = mpndfn(:,:,5,:)
    case('mpndtn1'); vwork = mpndtn(:,:,1,:)
    case('mpndtn2'); vwork = mpndtn(:,:,2,:)
    case('mpndtn3'); vwork = mpndtn(:,:,3,:)
    case('mpndtn4'); vwork = mpndtn(:,:,4,:)
    case('mpndtn5'); vwork = mpndtn(:,:,5,:)

    end select
    call gather_global(gwork, vwork, master_task, distrb_info)
  if (my_task == 0) then
    call write_nc2D(ncid, sumfld(jf), gwork, 2, il_im, il_jm, 1, ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck( nf_close(ncid) )

end subroutine save_restart_i2asum

!=================================================
subroutine save_restart_mice(fname, nstep)

! output ice variable averaged over the last IO cpl int of this run, 
! cice reads in these vars at the beginning of next run, uses them with the first 
! received a2i fields to obtain the first i2o fields to be sent to ocn  

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

if (my_task == 0) then
  call create_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
  call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
endif

!B: 20170825 ==> add maicen_saved for atm_icefluxes_back2GBM calculation!
!        note maicen_saved is the last ia interval mean.  
vwork(:,:,:) = maicen_saved(:,:,1,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen1', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,2,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen2', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,3,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen3', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,4,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen4', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,5,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen5', gwork, 2, il_im, il_jm, 1, ilout=il_out)
!b.

!The following fields are actually the ice state of last timestep 
!(no time-averaging is required in each timestep io coupling, see time_average_fields_4_i2o)
vwork = maice
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maice', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mstrocnxT
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnxT', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mstrocnyT
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnyT', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfresh
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfresh', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfsalt
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfsalt', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfhocn
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfhocn', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfswthru
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfswthru', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = msicemass
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'msicemass', gwork, 2, il_im, il_jm, 1, ilout=il_out)

if (my_task == 0) call ncheck( nf_close(ncid) )

return
end subroutine save_restart_mice

!=================================================
subroutine get_i2a_fields

implicit none

! all fields (except for vector) obtained here are all on T cell center

!(1) ocean surface temperature 
ia_sst(:,:,:) = msst(:,:,:)

!(2-3) ice/ocn combined surface velocity 
!CH: should use "aiu", not aice!
ia_uvel(:,:,:) = mssu(:,:,:) * (1. - maiu(:,:,:)) + muvel(:,:,:) * maiu(:,:,:) 
ia_vvel(:,:,:) = mssv(:,:,:) * (1. - maiu(:,:,:)) + mvvel(:,:,:) * maiu(:,:,:)

!(4-8) ice concentration
ia_aicen(:,:,:,:) = maicen(:,:,:,:)
!BX: save it for use in atm_icefluxes_back2GBM ---
maicen_saved = maicen

!XXX -- As per Alex West, only two of the ice vaiables below need to be scaled down 
!       by "* aice": ice top layer "temperature" and "effective conductivity"!

!(9-13) ice thickness
ia_thikn(:,:,:,:) = mthikn(:,:,:,:)
!ia_thikn(:,:,:,:) = mthikn(:,:,:,:) * mfoifr(:,:,:,:)  !X

!(14-18) snow thickness
ia_snown(:,:,:,:) = msnown(:,:,:,:)
!ia_snown(:,:,:,:) = msnown(:,:,:,:) * mfoifr(:,:,:,:)  !X

!(19-20) co2 flux stuff
ia_co2 = mco2
ia_co2fx = mco2fx

!(21) ocean surface freezing temperature
ia_sstfz(:,:,:) = msstfz(:,:,:) + 273.15

!(22-26) first order ice concentration
ia_foifr(:,:,:,:) = mfoifr(:,:,:,:)

!(27-31) ice top layer temperature
!XXX ia_itopt(:,:,:,:) = mitopt(:,:,:,:) + 273.15
ia_itopt(:,:,:,:) = (mitopt(:,:,:,:) + 273.15) * mfoifr(:,:,:,:)        !Y

!(32-36) ice top layer effective conductivity
!XXX ia_itopk(:,:,:,:) = mitopk(:,:,:,:)
ia_itopk(:,:,:,:) = mitopk(:,:,:,:) * mfoifr(:,:,:,:)                   !Y

!(37-41) ice melt pond concentration
ia_pndfn(:,:,:,:) = mpndfn(:,:,:,:)
!ia_pndfn(:,:,:,:) = mpndfn(:,:,:,:) * mfoifr(:,:,:,:)  !X

!(42-46) ice melt pond thickness
ia_pndtn(:,:,:,:) = mpndtn(:,:,:,:)
!ia_pndtn(:,:,:,:) = mpndtn(:,:,:,:) * mfoifr(:,:,:,:)  !X

return
end subroutine get_i2a_fields

!=================================================
subroutine get_i2o_fields

! All fluxes should be in GBM units before passing into coupler.
! e.g.,  io_htflx(:,:,:) = fhocn_gbm(:,:,:)	

implicit none

real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: pice

! Fields obtained here are all at T cell center. before being sent to MOM4, vector 
! (Taux, Tauy) should be shifted on to U point as required
!------------------------------------------------------------------------------- 

!(1-2) air/ice - sea stress TAUX/TAUY
!    Caution: in nemo, "strocnx/y" are NOT weighted by aice here, 'cos strocnx/y
!    have already been 'weighted' using aice (when calculated in "evp_finish". 
!    But, this weight has been removed in strocnx/yT (see "evp_finish"), therfore 
!    we need put it on again here. 
io_strsu = um_taux * (1. - maice) - mstrocnxT * maice
io_strsv = um_tauy * (1. - maice) - mstrocnyT * maice

!(3) freshwater flux to ocean: rainfall (+ ice melting water flux ?)
io_rain = um_rain * (1. - maice)
!!CH: confirmed:
!!if (ice_fwflux) io_rain = io_rain + mfresh	!always .t.
!!NOTE mfresh is now splitted into melt (14) and form (15) and passed into ocn seperately.

!(4) freshwater flux to ocean: snowfall
io_snow = um_snow * (1. - maice)

!(5) salt flux to ocean
io_stflx = mfsalt

!(6) ice/snow melting heatflux into ocean 
io_htflx = mfhocn

!(7) short wave radiation 
!(CH: the (1-aice) weight should not be here 'cos all fluxes passed in from
!     UM have already been aice-weighted when they are calculated there!!!) 
!io_swflx = um_swflx * (1. - maice) + mfswthru
io_swflx = um_swflx + mfswthru
!!!20100616: test for more swflx
!!!io_swflx = 1.064 * um_swflx + mfswthru

!(8) latent heat flux (positive out of ocean as required by MOM4)
io_qflux = um_evap * Lvap        !Note it's already weighted in UM for open sea.

!(9) sensible heat flux (positive out of ocean as required by MOM4)
io_shflx = um_shflx

!(10) net long wave radiation positive down
io_lwflx = um_lwflx

!(11) runoff (!check the incoming field! pattern? remapping ok? conserved? ...) 
io_runof = um_runoff 
! CHECK with SM about the annual cycle of core-runoff! (we only have annual mean)

!(12) pressure
pice = gravit * msicemass
!----------------------------------------------------------------------------
!sicemass = rho_ice x hi + rho_snow x hs (in m)
!
! Should we set limit to the ovelying ice pressure as suggested in MOM4 code?
!(see ocean_sbc.F90) if yes, we may use following 
!pice(i,j) = min(pice(i,j), gravit*rhow*max_ice_thickness) 
! (note  rhow = 1026 kg/m^3 here, but mom4 instead uses rho0 = 1035 kg/m^3)
! No, let mom4 handle it (see ocean_sbc.F90)
!
!as GFDL SIS, we use patm 'anormaly' and then add in the ice/snow pressure! 
!29/11/2007
!----------------------------------------------------------------------------
if (ice_pressure_on) then
  io_press = pice * maice
endif
if (air_pressure_on) then
   !as GFDL SIS, we use patm anormaly, i.e., taking off 1.e5 Pa !
   io_press(:,:,:) = io_press(:,:,:) + um_press(:,:,:) - 1.0e5 
endif
!(13) ice concentration
io_aice = maice
!(14) ice melt fwflux 
io_melt = max(0.0,mfresh(:,:,:)) 
!(15) ice form fwflux
io_form = min(0.0,mfresh(:,:,:))

io_co2 = um_co2
io_wnd = um_wnd

return
end subroutine get_i2o_fields

!=================================================
subroutine initialize_mice_fields_4_i2o

implicit none

maice = 0.
mstrocnxT = 0.
mstrocnyT = 0.
mfresh = 0.
mfsalt = 0.
mfhocn = 0. 
mfswthru = 0.
msicemass = 0.

return
end subroutine initialize_mice_fields_4_i2o

!=================================================
subroutine initialize_mice_fields_4_i2a

implicit none

muvel = 0.
mvvel = 0.

maiu = 0.
maicen = 0.
mthikn = 0.
msnown = 0. 

mfoifr = 0.
mitopt = 0.
mitopk = 0.
mpndfn = 0.
mpndtn = 0.

return
end subroutine initialize_mice_fields_4_i2a

!=================================================
subroutine initialize_mocn_fields_4_i2a

implicit none

msst = 0.
mssu = 0.
mssv = 0.
mco2 = 0.
mco2fx = 0.
msstfz = 0.

return
end subroutine initialize_mocn_fields_4_i2a

!=================================================
subroutine time_average_ocn_fields_4_i2a

implicit none

msst(:,:,:) = msst(:,:,:) + ocn_sst(:,:,:) * coef_ai
mssu(:,:,:) = mssu(:,:,:) + ocn_ssu(:,:,:) * coef_ai
mssv(:,:,:) = mssv(:,:,:) + ocn_ssv(:,:,:) * coef_ai
mco2(:,:,:) = mco2(:,:,:) + ocn_co2(:,:,:) * coef_ai
mco2fx(:,:,:) = mco2fx(:,:,:) + ocn_co2fx(:,:,:) * coef_ai
msstfz(:,:,:) = msstfz(:,:,:) + Tf(:,:,:) * coef_ai

return
end subroutine time_average_ocn_fields_4_i2a

subroutine time_average_fields_4_i2o
!now for each timestep io coupling, so no time-averaging is required. 
implicit none

maice(:,:,:)     = aice(:,:,:)
mstrocnxT(:,:,:) = strocnxT(:,:,:) 
mstrocnyT(:,:,:) = strocnyT(:,:,:) 
mfresh(:,:,:)    = fresh(:,:,:) 
mfsalt(:,:,:)    = fsalt(:,:,:) 
mfhocn(:,:,:)    = fhocn(:,:,:) 
mfswthru(:,:,:)  = fswthru(:,:,:)
msicemass(:,:,:) = sicemass(:,:,:)

return
end subroutine time_average_fields_4_i2o

!=================================================
subroutine time_average_fields_4_i2a

implicit none

! ice fields:
muvel(:,:,:) = muvel(:,:,:) + uvel(:,:,:) * coef_ai
mvvel(:,:,:) = mvvel(:,:,:) + vvel(:,:,:) * coef_ai
maicen(:,:,:,:) = maicen(:,:,:,:) + aicen(:,:,:,:) * coef_ai  !T cat. ice concentration
mthikn(:,:,:,:) = mthikn(:,:,:,:) + vicen(:,:,:,:) * coef_ai  !T cat. ice thickness
msnown(:,:,:,:) = msnown(:,:,:,:) + vsnon(:,:,:,:) * coef_ai  !T cat. snow thickness

call to_ugrid(aice, aiiu)
maiu(:,:,:)  = maiu(:,:,:) + aiiu(:,:,:) * coef_ai            !U cell ice concentraction

!BX: "First order" ice fraction (mfoifr, below) is required for GSI8 "Time-Travelling Ice" (TTI)
!    coupling approach. It may be different than the "normal" ice fraction (maicen, above) if  
!    maicen is regridded with second order conservation scheme (as "proposed" in GC3).
!    BUT, GC3 actually uses 1st order remapping for both of them, so they are identical! 
!    In ACCESS practice, no second order remapping has been appllied to any coupling field, and 
!    maicen and mfoifr are ALWAYS the same thing.
!    We pass both of them to UM for "concictency" (thus keeping UM coupling code intact)!
mfoifr(:,:,:,:) = mfoifr(:,:,:,:) + aicen(:,:,:,:)* coef_ai    !==maicen 
mitopt(:,:,:,:) = mitopt(:,:,:,:) + Tn_top(:,:,:,:) * coef_ai
mitopk(:,:,:,:) = mitopk(:,:,:,:) + keffn_top(:,:,:,:) * coef_ai
mpndfn(:,:,:,:) = mpndfn(:,:,:,:) + apeffn(:,:,:,:) * coef_ai
mpndtn(:,:,:,:) = mpndtn(:,:,:,:) + trcrn(:,:,nt_hpnd,:,:) * coef_ai 

!add one more a-i interval mean field (integrated ice concentration), which, togthere with maicen, 
!should be saved at the end of current run for use at the beginning of the continue run (e.g., 
!converting ice fluxes into GBM. see routines "atm_icefluxes_back2GBM", and "get_sbc_ice")...... 
!maice_ia(:,:,:) = maice_ia(:,:,:) + aice(:,:,:) * coef_ai

!ocn fields:
!must be done after calling from_ocn so as to get the most recently updated ocn fields,
!therefore a separate call to "time_average_ocn_fields_4_i2a" is done for this purpose.

return
end subroutine time_average_fields_4_i2a

!=================================================
subroutine check_i2a_fields(nstep)

implicit none

integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ilout, ll, jf
integer(kind=int_kind), save :: ncid,currstep 
data currstep/0/ 

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist('fields_i2a_in_ice.nc') ) then
  call create_ncfile('fields_i2a_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file fields_i2a_in_ice.nc at nstep = ', nstep
  call ncheck( nf_open('fields_i2a_in_ice.nc',nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid) 
end if

do jf = 1, nsend_i2a
   
  select case(trim(cl_writ(jf)))
    case('isst_ia');  vwork = ia_sst
    case('icecon01'); vwork(:,:,:) = ia_aicen(:,:,1,:) 
    case('icecon02'); vwork(:,:,:) = ia_aicen(:,:,2,:) 
    case('icecon03'); vwork(:,:,:) = ia_aicen(:,:,3,:) 
    case('icecon04'); vwork(:,:,:) = ia_aicen(:,:,4,:) 
    case('icecon05'); vwork(:,:,:) = ia_aicen(:,:,5,:) 
    case('snwthk01'); vwork(:,:,:) = ia_snown(:,:,1,:)
    case('snwthk02'); vwork(:,:,:) = ia_snown(:,:,2,:)  
    case('snwthk03'); vwork(:,:,:) = ia_snown(:,:,3,:)
    case('snwthk04'); vwork(:,:,:) = ia_snown(:,:,4,:)
    case('snwthk05'); vwork(:,:,:) = ia_snown(:,:,5,:)
    case('icethk01'); vwork(:,:,:) = ia_thikn(:,:,1,:)
    case('icethk02'); vwork(:,:,:) = ia_thikn(:,:,2,:)
    case('icethk03'); vwork(:,:,:) = ia_thikn(:,:,3,:)
    case('icethk04'); vwork(:,:,:) = ia_thikn(:,:,5,:)
    case('icethk05'); vwork(:,:,:) = ia_thikn(:,:,5,:)
    case('uvel_ia');  vwork = ia_uvel 
    case('vvel_ia');  vwork = ia_vvel 
    case('co2_i2');  vwork = ia_co2
    case('co2fx_i2');  vwork = ia_co2fx
    case('sstfz_ia'); vwork = ia_sstfz
    case('foifr01'); vwork(:,:,:) = ia_foifr(:,:,1,:)
    case('foifr02'); vwork(:,:,:) = ia_foifr(:,:,2,:)
    case('foifr03'); vwork(:,:,:) = ia_foifr(:,:,3,:)
    case('foifr04'); vwork(:,:,:) = ia_foifr(:,:,4,:)
    case('foifr05'); vwork(:,:,:) = ia_foifr(:,:,5,:)
    case('itopt01'); vwork(:,:,:) = ia_itopt(:,:,1,:)
    case('itopt02'); vwork(:,:,:) = ia_itopt(:,:,2,:)
    case('itopt03'); vwork(:,:,:) = ia_itopt(:,:,3,:)
    case('itopt04'); vwork(:,:,:) = ia_itopt(:,:,4,:)
    case('itopt05'); vwork(:,:,:) = ia_itopt(:,:,5,:)
    case('itopk01'); vwork(:,:,:) = ia_itopk(:,:,1,:)
    case('itopk02'); vwork(:,:,:) = ia_itopk(:,:,2,:)
    case('itopk03'); vwork(:,:,:) = ia_itopk(:,:,3,:)
    case('itopk04'); vwork(:,:,:) = ia_itopk(:,:,4,:)
    case('itopk05'); vwork(:,:,:) = ia_itopk(:,:,5,:)
    case('pndfn01'); vwork(:,:,:) = ia_pndfn(:,:,1,:)
    case('pndfn02'); vwork(:,:,:) = ia_pndfn(:,:,2,:)
    case('pndfn03'); vwork(:,:,:) = ia_pndfn(:,:,3,:)
    case('pndfn04'); vwork(:,:,:) = ia_pndfn(:,:,4,:)
    case('pndfn05'); vwork(:,:,:) = ia_pndfn(:,:,5,:)
    case('pndtn01'); vwork(:,:,:) = ia_pndtn(:,:,1,:)
    case('pndtn02'); vwork(:,:,:) = ia_pndtn(:,:,2,:)
    case('pndtn03'); vwork(:,:,:) = ia_pndtn(:,:,3,:)
    case('pndtn04'); vwork(:,:,:) = ia_pndtn(:,:,4,:)
    case('pndtn05'); vwork(:,:,:) = ia_pndtn(:,:,5,:)
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)
    
  if (my_task == 0 ) then
    call write_nc2D(ncid, trim(cl_writ(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_i2a_fields

!=================================================
subroutine check_a2i_fields(nstep)

implicit none

integer(kind=int_kind), intent(in) :: nstep
character*80 :: ncfile='fields_a2i_in_ice_2.nc'
integer(kind=int_kind) :: ncid, currstep, ll, ilout, jf
data currstep/0/
save currstep

currstep=currstep+1

if ( my_task == 0 .and. .not. file_exist(trim(ncfile)) ) then
  call create_ncfile(trim(ncfile),ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file ', trim(ncfile), ' at nstep = ', nstep
  call ncheck( nf_open(trim(ncfile),nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

do jf = 1, nrecv_a2i

  select case (trim(cl_read(jf)))
    case ('thflx_i');  vwork = um_thflx
    case ('pswflx_i'); vwork = um_pswflx
    case ('runoff_i'); vwork = um_runoff
    case ('wme_i');    vwork = um_wme
    case ('rain_i');   vwork = um_rain
    case ('snow_i');   vwork = um_snow
    case ('evap_i');   vwork = um_evap
    case ('lhflx_i');  vwork = um_lhflx
    case ('tmlt01'); vwork(:,:,:) = um_tmlt(:,:,1,:)
    case ('tmlt02'); vwork(:,:,:) = um_tmlt(:,:,2,:)
    case ('tmlt03'); vwork(:,:,:) = um_tmlt(:,:,3,:)
    case ('tmlt04'); vwork(:,:,:) = um_tmlt(:,:,4,:)
    case ('tmlt05'); vwork(:,:,:) = um_tmlt(:,:,5,:)
    case ('bmlt01'); vwork(:,:,:) = um_tmlt(:,:,1,:)
    case ('bmlt02'); vwork(:,:,:) = um_tmlt(:,:,2,:)
    case ('bmlt03'); vwork(:,:,:) = um_tmlt(:,:,3,:)
    case ('bmlt04'); vwork(:,:,:) = um_tmlt(:,:,4,:)
    case ('bmlt05'); vwork(:,:,:) = um_tmlt(:,:,5,:)
    case ('taux_i'); vwork = um_taux
    case ('tauy_i'); vwork = um_tauy
    case ('swflx_i'); vwork = um_swflx
    case ('lwflx_i'); vwork = um_lwflx
    case ('shflx_i'); vwork = um_shflx
    case ('press_i'); vwork = um_press
    case ('co2_ai'); vwork = um_co2
    case ('wnd_ai'); vwork = um_wnd
    case ('icenth_i'); vwork = um_icenth
    case ('icesth_i'); vwork = um_icesth
    case ('tsfice01'); vwork = um_tsfice(:,:,1,:)
    case ('tsfice02'); vwork = um_tsfice(:,:,2,:)
    case ('tsfice03'); vwork = um_tsfice(:,:,3,:)
    case ('tsfice04'); vwork = um_tsfice(:,:,4,:)
    case ('tsfice05'); vwork = um_tsfice(:,:,5,:)
    case ('iceevp01'); vwork = um_iceevp(:,:,1,:)
    case ('iceevp02'); vwork = um_iceevp(:,:,2,:)
    case ('iceevp03'); vwork = um_iceevp(:,:,3,:)
    case ('iceevp04'); vwork = um_iceevp(:,:,4,:)
    case ('iceevp05'); vwork = um_iceevp(:,:,5,:)
  end select 

  call gather_global(gwork, vwork, master_task, distrb_info)

  if (my_task == 0) then
    call write_nc2D(ncid, trim(cl_read(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_a2i_fields

!=================================================
subroutine check_i2o_fields(nstep, scale)

implicit none

integer(kind=int_kind), intent(in) :: nstep
real, intent(in) :: scale
integer(kind=int_kind) :: ncid, currstep, ll, ilout, jf
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist('fields_i2o_in_ice.nc') ) then
  call create_ncfile('fields_i2o_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file fields_i2o_in_ice.nc at nstep = ', nstep
  call ncheck( nf_open('fields_i2o_in_ice.nc',nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

do jf = nsend_i2a + 1, jpfldout

  select case(trim(cl_writ(jf)))
    case('strsu_io')
      vwork = scale * io_strsu
    case('strsv_io')
      vwork = scale * io_strsv
    case('rain_io')
      vwork = scale * io_rain
    case('snow_io')
      vwork = scale * io_snow
    case('stflx_io')
      vwork = scale * io_stflx
    case('htflx_io')
      vwork = scale * io_htflx
    case('swflx_io')
      vwork = scale * io_swflx
    case('qflux_io')
      vwork = scale * io_qflux
    case('shflx_io')
      vwork = scale * io_shflx
    case('lwflx_io')
      vwork = scale * io_lwflx
    case('runof_io')
      vwork = scale * io_runof
    case('press_io')
      vwork = scale * io_press
    case('aice_io')
      vwork = scale * io_aice
    case('form_io')
      vwork = scale * io_form
    case('melt_io')
      vwork = scale * io_melt
    case('co2_i1')
      vwork = scale * io_co2
    case('wnd_i1')
      vwork = scale * io_wnd
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)

  if (my_task == 0 ) then
    call write_nc2D(ncid, trim(cl_writ(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_i2o_fields

!=================================================
subroutine check_o2i_fields(nstep)

implicit none

integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid, currstep, ilout, ll, jf
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist('fields_o2i_in_ice.nc') ) then
  call create_ncfile('fields_o2i_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file fields_o2i_in_ice.nc at nstep = ', nstep
  call ncheck( nf_open('fields_o2i_in_ice.nc',nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

do jf = nrecv_a2i + 1, jpfldin

  select case (trim(cl_read(jf)))
    case ('sst_i'); vwork = ocn_sst
    case ('sss_i'); vwork = ocn_sss
    case ('ssu_i'); vwork = ocn_ssu
    case ('ssv_i'); vwork = ocn_ssv
    case ('sslx_i'); vwork = ocn_sslx
    case ('ssly_i'); vwork = ocn_ssly
    case ('pfmice_i'); vwork = ocn_pfmice
    case ('co2_oi'); vwork = ocn_co2
    case ('co2fx_oi'); vwork = ocn_co2fx
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)

  if (my_task == 0) then
    call write_nc2D(ncid, trim(cl_read(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_o2i_fields

!============================================================================
subroutine check_frzmlt_sst(ncfilenm)

!this is (mainly) used to check cice solo run frzmlt and sst !
! (for comparison against a coupled run forcing into cice)

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, sst, master_task, distrb_info) 
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, frzmlt, master_task, distrb_info) 
if (my_task == 0) call write_nc2D(ncid, 'frzmlt', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_frzmlt_sst

!=================================================
subroutine check_i2o_uvfluxes(ncfilenm)

!this is temporarily used to check i2o fields (uflux, vflux and maice) 
!for debug purpose

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, maice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maice', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, um_taux, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'um_taux', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, um_tauy, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'um_tauy', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, mstrocnxT, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnxT', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, mstrocnyT, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnyT', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_i2o_uvfluxes

!=================================================
subroutine check_ice_fields(ncfilenm)
!this is temporarily used to check ice fields immediately after ice_step 
!for debug purpose

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, aice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aice', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,1,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,2,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,3,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,4,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,5,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_ice_fields

!=================================================
subroutine check_ice_sbc_fields(ncfilenm)

!this is temporarily used to check ice_sbc fields got from get_ice_sbc 
!for debug purpose

implicit none

real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: v3d

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

v3d = 0.0

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, aice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aice', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, strax, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'strax', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, stray, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'stray', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = flatn_f(:,:,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = fcondtopn_f(:,:,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = fsurfn_f(:,:,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

call gather_global(gwork, fsnow, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsnow', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, frain, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'frain', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = trcrn(:,:,nt_Tsfc,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'trcrn1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = trcrn(:,:,nt_Tsfc,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'trcrn2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = trcrn(:,:,nt_Tsfc,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'trcrn3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = trcrn(:,:,nt_Tsfc,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'trcrn4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = trcrn(:,:,nt_Tsfc,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'trcrn5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

!from ocen:
call gather_global(gwork, frzmlt, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'frzmlt', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, sst, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, sss, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sss', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, uocn, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'uocn', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, vocn, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'vocn', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

call gather_global(gwork, ss_tltx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ss_tltx', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, ss_tlty, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ss_tlty', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

call gather_global(gwork, Tf, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'Tf', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_ice_sbc_fields

!=================================================
subroutine check_sstsss(ncfilenm)

!this is used to check cice sst/sss : temporary use (20091019)

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, sst, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, sss, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sss', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_sstsss

!=================================================
function file_exist (file_name)
!
character(len=*), intent(in) :: file_name
logical  file_exist

file_exist = .false.
if (len_trim(file_name) == 0) return
if (file_name(1:1) == ' ')    return

inquire (file=trim(file_name), exist=file_exist)

end function file_exist

!=================================================

end module cpl_forcing_handler
