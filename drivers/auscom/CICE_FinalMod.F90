!  SVN:$Id: CICE_FinalMod.F90 744 2013-09-27 22:53:24Z eclare $
!=======================================================================
!
!  This module contains routines for the final exit of the CICE model,
!  including final output and clean exit from any message passing
!  environments and frameworks.
!
!  authors: Philip W. Jones, LANL
!  2006: Converted to free source form (F90) by Elizabeth Hunke
!  2008: E. Hunke moved ESMF code to its own driver

      module CICE_FinalMod

      use ice_kinds_mod

      use accessom2_mod, only : accessom2_type => accessom2
      use cpl_interface, only : coupler_termination

      implicit none
      private
      public :: CICE_Finalize
      save

!=======================================================================

      contains

!=======================================================================
!
!  This routine shuts down CICE by exiting all relevent environments.

      subroutine CICE_Finalize(accessom2)

      use ice_atmo, only: dealloc_atmo
      use ice_brine, only: dealloc_brine
      use ice_dyn_eap, only: dealloc_dyn_eap
      use ice_dyn_shared, only: dealloc_dyn_shared
      use ice_flux, only: dealloc_flux
      use ice_forcing, only: dealloc_forcing
      use ice_grid, only: dealloc_grid
      use ice_meltpond_lvl, only: dealloc_meltpond_lvl
      use ice_shortwave, only: dealloc_shortwave
      use ice_state, only: dealloc_state
      use ice_therm_shared, only: dealloc_therm_shared
      use ice_zbgc_shared, only: dealloc_zbgc_shared

      use ice_exit, only: end_run
      use ice_fileunits, only: nu_diag, release_all_fileunits
      use ice_restart_shared, only: runid
      use ice_timers, only: ice_timer_stop, ice_timer_print_all, timer_total
      use ice_communicate, only: my_task, master_task
      use ice_calendar, only : year_init, nyr, month, mday, hour, sec
      use ice_calendar, only : calendar, time, dt

      type(accessom2_type), intent(inout) :: accessom2
      integer, dimension(6) :: date_array

   !-------------------------------------------------------------------
   ! stop timers and print timer info
   !-------------------------------------------------------------------

      call ice_timer_stop(timer_total)        ! stop timing entire run
      call ice_timer_print_all(stats=.true.) ! print timing information

!echmod      if (nu_diag /= 6) close (nu_diag) ! diagnostic output

      !-----------------------------------------------------------------
      ! Release the run-time-sized module arrays allocated in cice_init.
      ! Each dealloc_* is a no-op if the matching alloc_* never ran, so
      ! this is safe on an aborted or partial initialisation.
      !-----------------------------------------------------------------
      call dealloc_grid
      call dealloc_state
      call dealloc_flux
      call dealloc_forcing
      call dealloc_shortwave
      call dealloc_atmo
      call dealloc_brine
      call dealloc_dyn_shared
      call dealloc_dyn_eap
      call dealloc_meltpond_lvl
      call dealloc_therm_shared
      call dealloc_zbgc_shared

      call release_all_fileunits

   !-------------------------------------------------------------------
   ! write 'finished' file if needed
   !-------------------------------------------------------------------

      if (runid == 'bering') call writeout_finished_file()

   !-------------------------------------------------------------------
   ! quit MPI
   !-------------------------------------------------------------------

    ! Allow libaccessom2 to check that datetime of all models is synchronised at
    ! the end of the run.
    call calendar(time-dt)
    date_array(1) = nyr + year_init - 1
    date_array(2) = month
    date_array(3) = mday
    date_array(4) = int(sec / 3600)
    date_array(5) = int(mod(sec, 3600) / 60)
    date_array(6) = mod(sec, 60)
    call accessom2%deinit(cur_date_array=date_array)

   call coupler_termination

end subroutine CICE_Finalize

!=======================================================================
!
! Write a file indicating that this run finished cleanly.  This is
! needed only for runs on machine 'bering' (set using runid = 'bering').
!
!  author: Adrian Turner, LANL

      subroutine writeout_finished_file()
      
      use ice_restart_shared, only: restart_dir
      use ice_communicate, only: my_task, master_task

      character(len=char_len_long) :: filename

      if (my_task == master_task) then
           
         filename = trim(restart_dir)//"finished"
         open(11,file=filename)
         write(11,*) "finished"
         close(11)

      endif

      end subroutine writeout_finished_file

!=======================================================================

      end module CICE_FinalMod

!=======================================================================
