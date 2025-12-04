!  SVN:$Id: ice_history_write.F90 567 2013-01-07 02:57:36Z eclare $
!=======================================================================
!
! Writes history in using PIO
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
      module ice_history_write

      use ice_pio
      use pio
      use netcdf, only: NF90_UNLIMITED, NF90_CHUNKED, nf90_global, nf90_noerr

      use ice_kinds_mod
      use ice_constants, only: c0, c360, secday, spval, rad_to_deg
      use ice_blocks, only: nx_block, ny_block, block, get_block
      use ice_exit, only: abort_ice
      use ice_domain, only: distrb_info, nblocks, blocks_ice
      use ice_domain, only: equal_num_blocks_per_cpu
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_broadcast, only: broadcast_scalar
      use ice_gather_scatter, only: gather_global
      use ice_domain_size, only: nx_global, ny_global, max_nstrm, max_blocks
      use ice_grid, only: TLON, TLAT, ULON, ULAT, hm, bm, tarea, uarea, &
         dxu, dxt, dyu, dyt, HTN, HTE, ANGLE, ANGLET, &
         lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      use ice_history_shared
      use ice_itd, only: hin_max
      use ice_calendar, only: write_ic, histfreq
      use ice_fileunits, only: nu_diag, flush_fileunit

      implicit none
      private
      public :: ice_write_hist
      save

      type coord_attributes         ! netcdf coordinate attributes
         character (len=11)   :: short_name
         character (len=45)   :: long_name
         character (len=20)   :: units
      end type coord_attributes

      type req_attributes         ! req'd netcdf attributes
         type (coord_attributes) :: req
         character (len=20)   :: coordinates
      end type req_attributes

      type(io_desc_t)       :: iodesc2d, &
                               iodesc3dc, iodesc3dv, iodesc3di, iodesc3db, &
                               iodesc4di, iodesc4ds

      ! 4 coordinate variables: TLON, TLAT, ULON, ULAT
      INTEGER (kind=int_kind), PARAMETER :: ncoord = 4

      ! 4 vertices in each grid cell
      INTEGER (kind=int_kind), PARAMETER :: nverts = 4

      ! 4 variables describe T, U grid boundaries:
      ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 4

      contains

!=======================================================================

!
! write average ice quantities or snapshots
!
! author:   Elizabeth C. Hunke, LANL

      subroutine ice_write_hist (ns)

      use ice_calendar, only: time, sec, idate, idate0, &
        month, daymo, dayyr, days_per_year, use_leap_years

      integer (kind=int_kind), intent(in) :: ns !history stream number

      ! local variables

      real (kind=real_kind) :: ltime                 !history timestamp in seconds
      character (char_len_long) :: ncfile(max_nstrm), filename !filenames
      character (char_len) :: time_string            !model time for logging
      logical :: file_exists
      integer (kind=int_kind) :: & 
        i_time, &    ! time index
        timid       ! time var id
      type(file_desc_t)     :: File

      type(var_desc_t)      :: varid
      TYPE(req_attributes), dimension(nvar) :: var
      TYPE(coord_attributes), dimension(ncoord) :: coord_var
      TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
      TYPE(coord_attributes), dimension(nvarz) :: var_nz

      if (my_task == master_task) then
        ! set timestamp in middle of time interval
        if (histfreq(ns) == 'm' .or. histfreq(ns) == 'M') then
            if (month /= 1) then
                ltime=time/int(secday)-real(daymo(month-1))/2.0
            else
                ltime=time/int(secday)-real(daymo(12))/2.0
            endif
        else if(histfreq(ns) == 'd' .or. histfreq(ns) == 'D') then 
            ltime=time/int(secday) - 0.5
        else
            ltime=time/int(secday)
        endif

        call construct_filename(ncfile(ns),'nc',ns,time_string)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif

        ! run inquire only on master_task, to avoid possible race condition with
        ! multiple procs creating the file
        inquire(file=trim(ncfile(ns)),exist=file_exists)
      endif

      call broadcast_scalar(ncfile(ns), master_task)
      call broadcast_scalar(file_exists, master_task)

      File%fh=-1

      call ice_pio_initfile(mode='write', filename=trim(ncfile(ns)), File=File)

      if (.not. file_exists) then
         write(nu_diag,*) 'Writing dimensions and metadata: '//trim(ncfile(ns))
         call ice_hist_create(ns, ncfile(ns), File, var, coord_var, var_nverts, var_nz)
      endif

      !-----------------------------------------------------------------
      ! extent time coord by 1 and write time variable
      !-----------------------------------------------------------------
      call ice_pio_check(pio_inq_dimid(File, 'time', timid), &
                  'inq dimid time')
      call ice_pio_check(pio_inquire_dimension(File, timid, len=i_time), &
                  'inquire dim time')
      call ice_pio_check(pio_inq_varid(File,'time',varid), &
                  'inq varid time')
      i_time = i_time + 1 ! index of the current history time
      call ice_pio_check(pio_put_var(File,varid,(/i_time/),ltime), &
                  'put var time')

      !-----------------------------------------------------------------
      ! write time_bounds info
      !-----------------------------------------------------------------

      if (hist_avg) then
         time_bounds=(/time_beg(ns),time_end(ns)/)
         call ice_pio_check(pio_inq_varid(File,'time_bounds',varid), &
                     'inq varid time_bounds')
         call ice_pio_check(pio_put_var(File,varid,ival=(/time_beg(ns),time_end(ns)/), &
                     start=(/1,i_time/), count=(/2,1/)), &
                     'put var time_bounds')
      endif


      ! these iodesc variables describe how to many from the local variable in blocks, 
      ! to the output variable
      call ice_pio_initdecomp(iodesc=iodesc2d)
      call ice_pio_initdecomp(ndim3=nverts,    iodesc=iodesc3dv, inner_dim=.true.)
      call ice_pio_initdecomp(ndim3=ncat_hist, iodesc=iodesc3dc)
      call ice_pio_initdecomp(ndim3=nzilyr,     iodesc=iodesc3di)
      call ice_pio_initdecomp(ndim3=nzlyrb,    iodesc=iodesc3db)
      call ice_pio_initdecomp(ndim3=nzilyr,  ndim4=ncat_hist,  iodesc=iodesc4di)
      call ice_pio_initdecomp(ndim3=nzslyr,  ndim4=ncat_hist,  iodesc=iodesc4ds)

      if (i_time == 1) then
         ! these variables are time-invariant, only write once per file
         !-----------------------------------------------------------------
         ! write coordinate variables
         !-----------------------------------------------------------------
         call write_coordinate_variables(File, coord_var, var_nz)

         !-----------------------------------------------------------------
         ! write grid masks, area and rotation angle
         !-----------------------------------------------------------------
         call write_grid_variables(File, var, var_nverts)

      endif

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      call write_2d_variables(ns, File, i_time)
      call write_3d_and_4d_variables(ns, File, i_time)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      call pio_closefile(File)
      if (my_task == master_task) then
        write(nu_diag,*) 'Wrote ',trim(ncfile(ns)),' at time ',trim(time_string)
      endif

      !-----------------------------------------------------------------
      ! clean-up PIO descriptors
      !-----------------------------------------------------------------
      call pio_freedecomp(ice_pio_subsystem, iodesc2d)
      call pio_freedecomp(ice_pio_subsystem, iodesc3dv)
      call pio_freedecomp(ice_pio_subsystem, iodesc3dc)
      call pio_freedecomp(ice_pio_subsystem, iodesc3di)
      call pio_freedecomp(ice_pio_subsystem, iodesc3db)
      call pio_freedecomp(ice_pio_subsystem, iodesc4di)
      call pio_freedecomp(ice_pio_subsystem, iodesc4ds)

end subroutine ice_write_hist


subroutine ice_hist_create(ns, ncfile, File, var, coord_var, var_nverts, var_nz)

      use ice_calendar, only: idate, idate0, &
        dayyr, days_per_year, use_leap_years, histfreq_n, sec
      use ice_restart_shared, only: runid

      integer (kind=int_kind), intent(in) :: ns
      character (char_len_long), intent(in) :: ncfile
      type(file_desc_t), intent(inout)     :: File

      ! local variables

      integer (kind=int_kind) :: i,k,ic,n,nn, &
      status,imtid,jmtid,kmtidi,kmtids,kmtidb, cmtid,timid, &
      nvertexid,ivertex
      type(var_desc_t)      :: varid
      integer (kind=int_kind), dimension(3) :: dimid, dimid_nverts
      integer (kind=int_kind), dimension(4) :: dimidz, dimidex
      integer (kind=int_kind), dimension(5) :: dimidcz

      integer (kind=int_kind) :: deflate, deflate_level ! compression settings
      ! We leave shuffle at 0, this is only useful for integer data.
      integer (kind=int_kind), parameter :: shuffle = 0

      integer (kind=int_kind) :: ind,boundid

      character (char_len) :: title, start_time,current_date,current_time,time_period_freq
      character (len=8) :: cdate
      CHARACTER (char_len), dimension(ncoord) :: coord_bounds

      TYPE(req_attributes), dimension(nvar), intent(inout) :: var
      TYPE(coord_attributes), dimension(ncoord), intent(inout) :: coord_var
      TYPE(coord_attributes), dimension(nvar_verts), intent(inout) :: var_nverts
      TYPE(coord_attributes), dimension(nvarz), intent(inout) :: var_nz

      ! If history_deflate_level < 0 then don't do deflation,
      ! otherwise it sets the deflate level
      if (history_deflate_level < 0) then
         deflate = 0
         deflate_level = 0
      else
         deflate = 1
         deflate_level = history_deflate_level
      endif

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------
      call ice_pio_check(pio_def_dim(File, 'time', NF90_UNLIMITED, timid), &
                  'def dim time')

      if (hist_avg .and. histfreq(ns) /= '1') then
            call ice_pio_check(pio_def_dim(File,'d2',2,boundid), 'def dim d2')
      endif

      call ice_pio_check(pio_def_dim(File, 'ni', nx_global, imtid), &
                  'def dim ni')
      call ice_pio_check(pio_def_dim(File, 'nj', ny_global, jmtid), &
                  'def dim nj')
      call ice_pio_check(pio_def_dim(File, 'nc', ncat_hist, cmtid), &
                  'def dim nc')
      call ice_pio_check(pio_def_dim(File, 'nkice', nzilyr, kmtidi), &
                  'def dim nkice')
      call ice_pio_check(pio_def_dim(File, 'nksnow', nzslyr, kmtids), &
                  'def dim nksnow')
      call ice_pio_check(pio_def_dim(File, 'nkbio', nzblyr, kmtidb), &
                  'def dim nkbio')
      call ice_pio_check(pio_def_dim(File, 'nvertices', nverts, nvertexid), &
                  'def dim nverts')

      !-----------------------------------------------------------------
      ! define coordinate variables
      !-----------------------------------------------------------------

      call ice_pio_check(pio_def_var(File,'time',pio_real,(/timid/),varid), &
                  'def var time')
      call ice_pio_check(pio_put_att(File,varid,'long_name','model time'), &
                  'put_att long_name')

      write(cdate,'(i8.8)') idate0
      write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
            cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
      call ice_pio_check(pio_put_att(File,varid,'units',title), &
                  'put_att time units')

      if (days_per_year == 360) then
            call ice_pio_check(pio_put_att(File,varid,'calendar','360_day'), &
                        'att time calendar')
      elseif (days_per_year == 365 .and. .not.use_leap_years ) then
            call ice_pio_check(pio_put_att(File,varid,'calendar','NoLeap'), &
                        'att time calendar')
      elseif (use_leap_years) then
            call ice_pio_check(pio_put_att(File,varid,'calendar','proleptic_gregorian'), &
                        'att time calendar')
      else
            call abort_ice( 'ice Error: invalid calendar settings')
      endif

      if (hist_avg .and. histfreq(ns) /= '1' ) then
            call ice_pio_check(pio_put_att(File,varid,'bounds','time_bounds'), &
                        'att time bounds')
      endif

      !-----------------------------------------------------------------
      ! Define attributes for time bounds if hist_avg is true
      !-----------------------------------------------------------------

      if (hist_avg .and. histfreq(ns) /= '1') then
            dimid(1) = boundid
            dimid(2) = timid
            call ice_pio_check(pio_def_var(File, 'time_bounds', &
                                 pio_real,dimid(1:2),varid), &
                        'def var time_bounds')

            call ice_pio_check(pio_put_att(File,varid,'long_name', &
                              'boundaries for time-averaging interval'), &
                        'att time_bounds long_name')
            write(cdate,'(i8.8)') idate0
            write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
                  cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
            call ice_pio_check(pio_put_att(File,varid,'units',title), &
                                 'att time_bounds units')
      endif

      !-----------------------------------------------------------------
      ! define information for required time-invariant variables
      !-----------------------------------------------------------------

      ind = 0
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLON', &
                        'T grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lont_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLAT', &
                        'T grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latt_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULON', &
                        'U grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lonu_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULAT', &
                        'U grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latu_bounds'

      var_nz(1) = coord_attributes('NCAT', 'category maximum thickness', 'm')
      var_nz(2) = coord_attributes('VGRDi', 'vertical ice levels', '1')
      var_nz(3) = coord_attributes('VGRDs', 'vertical snow levels', '1')
      var_nz(4) = coord_attributes('VGRDb', 'vertical ice-bio levels', '1')

      !-----------------------------------------------------------------
      ! define information for optional time-invariant variables
      !-----------------------------------------------------------------

      var(n_tmask)%req = coord_attributes('tmask', &
                  'ocean grid mask', ' ')
      var(n_tmask)%coordinates = 'TLON TLAT'

      var(n_blkmask)%req = coord_attributes('blkmask', &
                  'ice grid block mask', ' ')
      var(n_blkmask)%coordinates = 'TLON TLAT'

      var(n_tarea)%req = coord_attributes('tarea', &
                  'area of T grid cells', 'm^2')
      var(n_tarea)%coordinates = 'TLON TLAT'

      var(n_uarea)%req = coord_attributes('uarea', &
                  'area of U grid cells', 'm^2')
      var(n_uarea)%coordinates = 'ULON ULAT'
      var(n_dxt)%req = coord_attributes('dxt', &
                  'T cell width through middle', 'm')
      var(n_dxt)%coordinates = 'TLON TLAT'
      var(n_dyt)%req = coord_attributes('dyt', &
                  'T cell height through middle', 'm')
      var(n_dyt)%coordinates = 'TLON TLAT'
      var(n_dxu)%req = coord_attributes('dxu', &
                  'U cell width through middle', 'm')
      var(n_dxu)%coordinates = 'ULON ULAT'
      var(n_dyu)%req = coord_attributes('dyu', &
                  'U cell height through middle', 'm')
      var(n_dyu)%coordinates = 'ULON ULAT'
      var(n_HTN)%req = coord_attributes('HTN', &
                  'T cell width on North side','m')
      var(n_HTN)%coordinates = 'TLON TLAT'
      var(n_HTE)%req = coord_attributes('HTE', &
                  'T cell width on East side', 'm')
      var(n_HTE)%coordinates = 'TLON TLAT'
      var(n_ANGLE)%req = coord_attributes('ANGLE', &
                  'angle grid makes with latitude line on U grid', &
                  'radians')
      var(n_ANGLE)%coordinates = 'ULON ULAT'
      var(n_ANGLET)%req = coord_attributes('ANGLET', &
                  'angle grid makes with latitude line on T grid', &
                  'radians')
      var(n_ANGLET)%coordinates = 'TLON TLAT'

      ! These fields are required for CF compliance
      ! dimensions (nx,ny,nverts)
      var_nverts(n_lont_bnds) = coord_attributes('lont_bounds', &
                  'longitude boundaries of T cells', 'degrees_east')
      var_nverts(n_latt_bnds) = coord_attributes('latt_bounds', &
                  'latitude boundaries of T cells', 'degrees_north')
      var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds', &
                  'longitude boundaries of U cells', 'degrees_east')
      var_nverts(n_latu_bnds) = coord_attributes('latu_bounds', &
                  'latitude boundaries of U cells', 'degrees_north')

      !-----------------------------------------------------------------
      ! define attributes for time-invariant variables
      !-----------------------------------------------------------------

      dimid(1) = imtid
      dimid(2) = jmtid
      dimid(3) = timid

      do i = 1, ncoord
         call ice_pio_check(pio_def_var(File, coord_var(i)%short_name, pio_real, &
                                 dimid(1:2), varid), &
                     'def var '//coord_var(i)%short_name)

         if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
            call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                        (/ history_chunksize_x, history_chunksize_y /)), &
                     'def var chunking '//coord_var(i)%short_name)
         endif

         call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, deflate_level), &
                     'deflate '//coord_var(i)%short_name)

         call ice_pio_check(pio_put_att(File, varid,'long_name',coord_var(i)%long_name), &
                     'put att long_name '//coord_var(i)%short_name)
         call ice_pio_check(pio_put_att(File, varid, 'units', coord_var(i)%units), &
                     'put att units '//coord_var(i)%short_name)
         call ice_pio_check(pio_put_att(File, varid,'missing_value',spval), &
                     'put att missing_value '//coord_var(i)%short_name)

         call ice_pio_check(pio_put_att(File, varid, '_FillValue', spval), &
                     'put att _FillValue '//coord_var(i)%short_name)

         if (coord_var(i)%short_name == 'ULAT') then
            call ice_pio_check(pio_put_att(File,varid,'comment', &
                                    'Latitude of NE corner of T grid cell'), &
                        'put att comment for '//coord_var(i)%short_name)
         endif
         if (f_bounds) then
            call ice_pio_check(pio_put_att(File, varid, 'bounds', coord_bounds(i)), &
                        'put att bounds '//coord_var(i)%short_name)
         endif
      enddo

      ! Extra dimensions (NCAT, NZILYR, NZSLYR, NZBLYR)
      dimidex(1)=cmtid
      dimidex(2)=kmtidi
      dimidex(3)=kmtids
      dimidex(4)=kmtidb

      do i = 1, nvarz
            if (igrdz(i)) then
               call ice_pio_check(pio_def_var(File, var_nz(i)%short_name, &
                                    pio_real, (/dimidex(i)/), varid), &
                           'def var '//var_nz(i)%short_name)

               call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, &
                                                deflate_level), &
                           'deflate '//var_nz(i)%short_name)

               call ice_pio_check(pio_put_att(File, varid,'long_name',var_nz(i)%long_name), &
                           'put att long_name '//var_nz(i)%short_name)
               call ice_pio_check(pio_put_att(File, varid, 'units', var_nz(i)%units), &
                           'for att units '//var_nz(i)%short_name)
            endif
      enddo

      ! Attributes for tmask, blkmask defined separately, since they have no units
      if (igrd(n_tmask)) then
            call ice_pio_check(pio_def_var(File, 'tmask', pio_real, dimid(1:2), varid), &
                     'def var tmask')

            if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
               call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                           (/ history_chunksize_x, history_chunksize_y /)), &
                           'def var chunking tmask')
            endif

            call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, &
                                       deflate_level), 'deflating var tmask')

            call ice_pio_check(pio_put_att(File,varid, 'long_name', 'ocean grid mask'), &
                     'put att tmask long_name')
            call ice_pio_check(pio_put_att(File, varid, 'coordinates', 'TLON TLAT'), &
                     'put att tmask units')
            call ice_pio_check(pio_put_att(File,varid,'comment', '0 = land, 1 = ocean'), &
                     'put att tmask comment')
            call ice_pio_check(pio_put_att(File,varid,'missing_value',spval), &
                     'put att missing_value for tmask')
            call ice_pio_check(pio_put_att(File,varid,'_FillValue',spval), &
                     'put att _FillValue for tmask')
      endif

      if (igrd(n_blkmask)) then
            call ice_pio_check(pio_def_var(File, 'blkmask', pio_real, dimid(1:2), varid), &
                     'def var blkmask')

            if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
               call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                           (/ history_chunksize_x, history_chunksize_y /)), &
                           'def var chunking blkmask')
            endif

            call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, &
                                          deflate_level), &
                     'deflating var blkmask')

            call ice_pio_check(pio_put_att(File,varid, 'long_name', 'ice grid block mask'), &
                        'put att blkmask long_name')
            call ice_pio_check(pio_put_att(File, varid, 'coordinates', 'TLON TLAT'), &
                        'put att blkmask coordinates')
            call ice_pio_check(pio_put_att(File,varid,'comment', 'mytask + iblk/100'), &
                        'put att blkmask comment')
            call ice_pio_check(pio_put_att(File,varid,'missing_value',spval), &
                        'put att blkmask missing_value')
            call ice_pio_check(pio_put_att(File,varid,'_FillValue',spval), &
                        'put att blkmask _FillValue')
      endif

      do i = 3, nvar      ! note n_tmask=1, n_blkmask=2
            if (igrd(i)) then
               call ice_pio_check(pio_def_var(File, var(i)%req%short_name, &
                                       pio_real, dimid(1:2), varid), &
                           'def variable '//var(i)%req%short_name)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y /)), &
                              'def var chunking '//var(i)%req%short_name)
               endif

               call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level), &
                           'deflate var '//var(i)%req%short_name)

               call ice_pio_check(pio_put_att(File,varid, 'long_name', var(i)%req%long_name), &
                           'put att long_name '//var(i)%req%short_name)
               call ice_pio_check(pio_put_att(File, varid, 'units', var(i)%req%units), &
                           'put att units '//var(i)%req%short_name)
               call ice_pio_check(pio_put_att(File, varid, 'coordinates', var(i)%coordinates), &
                           'put att coordinates '//var(i)%req%short_name)
               call ice_pio_check(pio_put_att(File,varid,'missing_value',spval), &
                           'put att missing_value '//var(i)%req%short_name)
               call ice_pio_check(pio_put_att(File,varid,'_FillValue',spval), &
                           'put att _FillValue '//var(i)%req%short_name)
            endif
      enddo

      ! Fields with dimensions (nverts,nx,ny)
      dimid_nverts(1) = nvertexid
      dimid_nverts(2) = imtid
      dimid_nverts(3) = jmtid
      do i = 1, nvar_verts
            if (f_bounds) then
               call ice_pio_check(pio_def_var(File, var_nverts(i)%short_name, &
                                    pio_real,dimid_nverts, varid), &
                           'def var '//var_nverts(i)%short_name)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ 1, history_chunksize_x, history_chunksize_y /)), &
                              'def var chunking '//var_nverts(i)%short_name)
               endif

               call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level), &
                           'deflate var '//var_nverts(i)%short_name)

               call ice_pio_check(pio_put_att(File,varid, 'long_name', &
                           var_nverts(i)%long_name), &
                           'put att long_name '//var_nverts(i)%short_name)
               call ice_pio_check(pio_put_att(File, varid, 'units', var_nverts(i)%units), &
                           'put att units '//var_nverts(i)%short_name)
               call ice_pio_check(pio_put_att(File,varid,'missing_value',spval), &
                           'put att missing_value '//var_nverts(i)%short_name)
               call ice_pio_check(pio_put_att(File,varid,'_FillValue',spval), &
                           'put att _FillValue '//var_nverts(i)%short_name)
            endif
      enddo

      do n=1,num_avail_hist_fields_2D
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_pio_check(pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimid, varid), &
                           'def var '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               call ice_pio_check(pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level), &
                           'deflate '//avail_hist_fields(n)%vname)

               call ice_pio_check(pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit), &
                           'put att units '//avail_hist_fields(n)%vname)
               call ice_pio_check(pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc), &
                           'put att long_name '//avail_hist_fields(n)%vname)
               call ice_pio_check(pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord), &
                           'put att coordinates '//avail_hist_fields(n)%vname)
               call ice_pio_check(pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas), &
                           'put att cell_measures '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                            avail_hist_fields(n)%vcomment), &
                            'put att comment '//avail_hist_fields(n)%vname)
               endif
               call ice_pio_check(pio_put_att(File,varid,'missing_value',spval), &
                           'put att missing_value '//avail_hist_fields(n)%vname)
               call ice_pio_check(pio_put_att(File,varid,'_FillValue',spval), &
                           'put att _FillValue '//avail_hist_fields(n)%vname)

               !---------------------------------------------------------------
               ! Add cell_methods attribute to variables if averaged
               !---------------------------------------------------------------
               if (hist_avg) then
                  if (TRIM(avail_hist_fields(n)%vname)/='sig1' .or. &
                     TRIM(avail_hist_fields(n)%vname)/='sig2') then
                     call ice_pio_check(pio_put_att(File,varid,'cell_methods','time: mean'), &
                              'put att cell methods time mean '//avail_hist_fields(n)%vname)
                  endif
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg         &
                  .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                  .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                  .or. n==n_trsig(ns)                             &
                  .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                  .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
                  call ice_pio_check(pio_put_att(File,varid,'time_rep','instantaneous'), &
                              'put att time_rep instantaneous')
               else
                  call ice_pio_check(pio_put_att(File,varid,'time_rep','averaged'), &
                              'put att time_rep averaged')
               endif
            endif
      enddo  ! num_avail_hist_fields_2D

         ! 3D (category)

      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = cmtid
      dimidz(4) = timid

      do n = n2D + 1, n3Dccum
         ! to-do: use check subroutine
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimidz, varid)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1, 1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               status = pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error deflating variable '//avail_hist_fields(n)%vname)

               status = pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining cell measures for '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                     avail_hist_fields(n)%vcomment), &
                     'put att comment '//avail_hist_fields(n)%vname)
               endif
               status = pio_put_att(File,varid,'missing_value',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining missing_value for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'_FillValue',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining _FillValue for '//avail_hist_fields(n)%vname)

               !---------------------------------------------------------------
               ! Add cell_methods attribute to variables if averaged
               !---------------------------------------------------------------
               if (hist_avg .and. histfreq(ns) /= '1') then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
                  if (status /= nf90_noerr) call abort_ice( &
                     'Error defining cell methods for '//avail_hist_fields(n)%vname)
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg) then
                  status = pio_put_att(File,varid,'time_rep','instantaneous')
               else
                  status = pio_put_att(File,varid,'time_rep','averaged')
               endif
            endif
      enddo  ! num_avail_hist_fields_3Dc

         ! 3D (ice layers)

      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = kmtidi
      dimidz(4) = timid

      do n = n3Dccum + 1, n3Dzcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimidz, varid)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1, 1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               status = pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error deflating variable '//avail_hist_fields(n)%vname)

               status = pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining cell measures for '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                     avail_hist_fields(n)%vcomment), &
                     'put att comment '//avail_hist_fields(n)%vname)
               endif
               status = pio_put_att(File,varid,'missing_value',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining missing_value for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'_FillValue',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining _FillValue for '//avail_hist_fields(n)%vname)

               ! Add cell_methods attribute to variables if averaged
               if (hist_avg .and. histfreq(ns) /= '1') then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg) then
                  status = pio_put_att(File,varid,'time_rep','instantaneous')
               else
                  status = pio_put_att(File,varid,'time_rep','averaged')
               endif

            endif
      enddo  ! num_avail_hist_fields_3Dz

         ! 3D (biology layers)

      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = kmtidb
      dimidz(4) = timid

      do n = n3Dzcum + 1, n3Dbcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimidz, varid)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1, 1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               status = pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error deflating variable '//avail_hist_fields(n)%vname)

               status = pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining cell measures for '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                     avail_hist_fields(n)%vcomment), &
                     'put att comment '//avail_hist_fields(n)%vname)
               endif
               status = pio_put_att(File,varid,'missing_value',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining missing_value for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'_FillValue',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining _FillValue for '//avail_hist_fields(n)%vname)

               ! Add cell_methods attribute to variables if averaged
               if (hist_avg .and. histfreq(ns) /= '1') then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg) then
                  status = pio_put_att(File,varid,'time_rep','instantaneous')
               else
                  status = pio_put_att(File,varid,'time_rep','averaged')
               endif

            endif
      enddo  ! num_avail_hist_fields_3Db

         ! 4D (ice categories and layers)

      dimidcz(1) = imtid
      dimidcz(2) = jmtid
      dimidcz(3) = kmtidi
      dimidcz(4) = cmtid
      dimidcz(5) = timid

      do n = n3Dbcum + 1, n4Dicum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimidcz, varid)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1, 1, 1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               status = pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error deflating variable '//avail_hist_fields(n)%vname)

               status = pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining cell measures for '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                     avail_hist_fields(n)%vcomment), &
                     'put att comment '//avail_hist_fields(n)%vname)
               endif
               status = pio_put_att(File,varid,'missing_value',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining missing_value for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'_FillValue',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining _FillValue for '//avail_hist_fields(n)%vname)

               !---------------------------------------------------------------
               ! Add cell_methods attribute to variables if averaged
               !---------------------------------------------------------------
               if (hist_avg .and. histfreq(ns) /= '1') then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
                  if (status /= nf90_noerr) call abort_ice( &
                        'Error defining cell methods for '//avail_hist_fields(n)%vname)
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg) then
                  status = pio_put_att(File,varid,'time_rep','instantaneous')
               else
                  status = pio_put_att(File,varid,'time_rep','averaged')
               endif
            endif
      enddo  ! num_avail_hist_fields_4Di


      ! 4D (ice categories and snow layers)

      dimidcz(1) = imtid
      dimidcz(2) = jmtid
      dimidcz(3) = kmtids
      dimidcz(4) = cmtid
      dimidcz(5) = timid

      do n = n4Dicum + 1, n4Dscum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimidcz, varid) 
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1, 1 ,1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               status = pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error deflating variable '//avail_hist_fields(n)%vname)

               status = pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining cell measures for '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                     avail_hist_fields(n)%vcomment), &
                     'put att comment '//avail_hist_fields(n)%vname)
               endif
               status = pio_put_att(File,varid,'missing_value',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining missing_value for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'_FillValue',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining _FillValue for '//avail_hist_fields(n)%vname)

               !---------------------------------------------------------------
               ! Add cell_methods attribute to variables if averaged
               !---------------------------------------------------------------
               if (hist_avg .and. histfreq(ns) /= '1') then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
                  if (status /= nf90_noerr) call abort_ice( &
                     'Error defining cell methods for '//avail_hist_fields(n)%vname)
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg) then
                  status = pio_put_att(File,varid,'time_rep','instantaneous')
               else
                  status = pio_put_att(File,varid,'time_rep','averaged')
               endif
            endif
      enddo  ! num_avail_hist_fields_4Ds

      ! 4D (ice categories and biology layers)

      dimidcz(1) = imtid
      dimidcz(2) = jmtid
      dimidcz(3) = kmtidb
      dimidcz(4) = cmtid
      dimidcz(5) = timid

      do n = n4Dscum + 1, n4Dbcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                                       pio_real, dimidcz, varid)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//avail_hist_fields(n)%vname)

               if (history_chunksize_x > 0 .and. history_chunksize_y > 0) then
                  call ice_pio_check(pio_def_var_chunking(File, varid, NF90_CHUNKED, &
                              (/ history_chunksize_x, history_chunksize_y, 1, 1 /)), &
                              'def var chunking '//avail_hist_fields(n)%vname)
               endif

               status = pio_def_var_deflate(File, varid, shuffle, deflate, &
                                             deflate_level)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error deflating variable '//avail_hist_fields(n)%vname)

               status = pio_put_att(File,varid,'units', &
                           avail_hist_fields(n)%vunit)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid, 'long_name', &
                           avail_hist_fields(n)%vdesc)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'coordinates', &
                           avail_hist_fields(n)%vcoord)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'cell_measures', &
                           avail_hist_fields(n)%vcellmeas)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining cell measures for '//avail_hist_fields(n)%vname)
               if (avail_hist_fields(n)%vcomment /= "none") then
                  call ice_pio_check(pio_put_att(File,varid,'comment', &
                     avail_hist_fields(n)%vcomment), &
                     'put att comment '//avail_hist_fields(n)%vname)
               endif
               status = pio_put_att(File,varid,'missing_value',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining missing_value for '//avail_hist_fields(n)%vname)
               status = pio_put_att(File,varid,'_FillValue',spval)
               if (status /= nf90_noerr) call abort_ice( &
                  'Error defining _FillValue for '//avail_hist_fields(n)%vname)

               !---------------------------------------------------------------
               ! Add cell_methods attribute to variables if averaged
               !---------------------------------------------------------------
               if (hist_avg .and. histfreq(ns) /= '1' ) then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
                  if (status /= nf90_noerr) call abort_ice( &
                     'Error defining cell methods for '//avail_hist_fields(n)%vname)
               endif

               if (histfreq(ns) == '1' .or. .not. hist_avg) then
                  call ice_pio_check(pio_put_att(File,varid,'time_rep','instantaneous'), &
                              'put att time_rep instantaneous')
               else
                  call ice_pio_check(pio_put_att(File,varid,'time_rep','averaged'), &
                              'put att time_rep averaged')
               endif
            endif
      enddo  ! num_avail_hist_fields_4Db

      title  = 'sea ice model output for CICE'
      call ice_pio_check(pio_put_att(File,nf90_global,'title',title), &
                  'global attribute title')

      title = 'Diagnostic and Prognostic Variables'
      call ice_pio_check(pio_put_att(File,nf90_global,'contents',title), &
                  'global attribute contents')

      title  = 'Los Alamos Sea Ice Model (CICE) Version 5'
      call ice_pio_check(pio_put_att(File,nf90_global,'source',title), &
                  'global attribute source')

      select case (histfreq(ns))
         case ("y", "Y")
            write(time_period_freq,'(a,i0)') 'year_',histfreq_n(ns)
         case ("m", "M")
            write(time_period_freq,'(a,i0)') 'month_',histfreq_n(ns)
         case ("d", "D")
            write(time_period_freq,'(a,i0)') 'day_',histfreq_n(ns)
         case ("h", "H")
            write(time_period_freq,'(a,i0)') 'hour_',histfreq_n(ns)
         case ("1")
            write(time_period_freq,'(a,i0)') 'step_',histfreq_n(ns)
      end select

      status = pio_put_att(File,pio_global,'time_period_freq',trim(time_period_freq))

#ifdef AusCOM
      write(title,'(a,i3,a)') 'This Year Has ',int(dayyr),' days'
#else
      if (use_leap_years) then
            write(title,'(a,i3,a)') 'This year has ',int(dayyr),' days'
      else
            write(title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
      endif
#endif
      call ice_pio_check(pio_put_att(File,nf90_global,'comment',title), &
                  'global attribute comment')

      write(title,'(a,i8.8)') 'File started on model date ',idate
      call ice_pio_check(pio_put_att(File,nf90_global,'comment2',title), &
                  'global attribute date1')

      ! TO-DO: Update output for CF compliance !
      ! title = 'CF-1.0'
      ! call ice_pio_check(pio_put_att(File,nf90_global,'conventions',title), &
                  ! 'global attribute conventions')

      call date_and_time(date=current_date, time=current_time)
      write(start_time,1000) current_date(1:4), current_date(5:6), &
                              current_date(7:8), current_time(1:2), &
                              current_time(3:4), current_time(5:8)
1000    format('This dataset was created on ', &
               a,'-',a,'-',a,' at ',a,':',a,':',a)

      call ice_pio_check(pio_put_att(File,nf90_global,'history',start_time), &
                  'global attribute history')

      call ice_pio_check(pio_put_att(File,nf90_global,'io_flavor','io_pio'), &
                  'global attribute io_flavor')

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

      call ice_pio_check(pio_enddef(File), 'enddef')

end subroutine ice_hist_create

subroutine write_coordinate_variables(File, coord_var, var_nz)

      type(file_desc_t), intent(inout)     :: File
      type(coord_attributes), dimension(ncoord), intent(in) :: coord_var
      type(coord_attributes), dimension(nvarz) :: var_nz

      real (kind=real_kind), allocatable :: workr2(:,:,:)

      integer :: i, k, status
      type(var_desc_t)      :: varid
      character (len=len(coord_var(1)%short_name)) :: coord_var_name

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------
        call pio_seterrorhandling(File, PIO_RETURN_ERROR)
        allocate(workr2(nx_block,ny_block,nblocks))

        do i = 1,ncoord
          call ice_pio_check(pio_inq_varid(File, coord_var(i)%short_name, varid), &
                     'inquire varid '//coord_var(i)%short_name)
          SELECT CASE (coord_var(i)%short_name)
            CASE ('TLON')
              ! Convert T grid longitude from -180 -> 180 to 0 to 360
                 workr2(:,:,:) = mod(tlon(:,:,1:nblocks)*rad_to_deg + c360, c360)
            CASE ('TLAT')
              workr2(:,:,:) = tlat(:,:,1:nblocks)*rad_to_deg
            CASE ('ULON')
              workr2(:,:,:) = ulon(:,:,1:nblocks)*rad_to_deg
            CASE ('ULAT')
              workr2(:,:,:) = ulat(:,:,1:nblocks)*rad_to_deg
          END SELECT
          call pio_write_darray(File, varid, iodesc2d, &
               workr2, status, fillval=spval)
          call ice_pio_check(status, 'write darray '//coord_var(i)%short_name)
        enddo

        ! Extra dimensions (NCAT, VGRD*)

        do i = 1, nvarz
          if (igrdz(i)) then
            status = pio_inq_varid(File, var_nz(i)%short_name, varid)
            SELECT CASE (var_nz(i)%short_name)
              CASE ('NCAT')
                status = pio_put_var(File, varid, hin_max(1:ncat_hist)) 
              CASE ('VGRDi')
                status = pio_put_var(File, varid, (/(k, k=1,nzilyr)/))
              CASE ('VGRDs')
                status = pio_put_var(File, varid, (/(k, k=1,nzslyr)/))
              CASE ('VGRDb')
                status = pio_put_var(File, varid, (/(k, k=1,nzblyr)/))
             END SELECT
           endif
        enddo
        call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)


end subroutine write_coordinate_variables



subroutine write_grid_variables(File, var, var_nverts)

   type(file_desc_t), intent(inout)    :: File
   type(req_attributes), dimension(nvar), intent(in) :: var
   type(coord_attributes), dimension(nvar_verts), intent(in) :: var_nverts

   real (kind=real_kind), allocatable :: workr2(:,:,:)
   real (kind=real_kind), allocatable :: workr3v(:,:,:,:)

   integer :: ivertex, i, status
   type(var_desc_t)      :: varid
   character (len=len(var(1)%req%short_name)) :: var_name
   character (len=len(var_nverts(1)%short_name)) :: var_nverts_name

   allocate(workr2(nx_block,ny_block,nblocks))

   do i = 1, nvar       ! note: n_tmask=1, n_blkmask=2
         if (igrd(i)) then
            SELECT CASE (var(i)%req%short_name)
            CASE ('tmask')
               workr2 = hm(:,:,1:nblocks)
            CASE ('blkmask')
               workr2 = bm(:,:,1:nblocks)
            CASE ('tarea')
               workr2 = tarea(:,:,1:nblocks)
            CASE ('uarea')
               workr2 = uarea(:,:,1:nblocks)
            CASE ('dxu')
               workr2 = dxu(:,:,1:nblocks)
            CASE ('dyu')
               workr2 = dyu(:,:,1:nblocks)
            CASE ('dxt')
               workr2 = dxt(:,:,1:nblocks)
            CASE ('dyt')
               workr2 = dyt(:,:,1:nblocks)
            CASE ('HTN')
               workr2 = HTN(:,:,1:nblocks)
            CASE ('HTE')
               workr2 = HTE(:,:,1:nblocks)
            CASE ('ANGLE')
               workr2 = ANGLE(:,:,1:nblocks)
            CASE ('ANGLET')
               workr2 = ANGLET(:,:,1:nblocks)
            END SELECT
            status = pio_inq_varid(File, var(i)%req%short_name, varid)
            call pio_write_darray(File, varid, iodesc2d, &
                 workr2, status, fillval=spval)
         endif
      enddo

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds) then
      allocate(workr3v(nverts,nx_block,ny_block,nblocks))
      workr3v (:,:,:,:) = c0
      do i = 1, nvar_verts
        SELECT CASE (var_nverts(i)%short_name)
        CASE ('lont_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = lont_bounds(ivertex,:,:,1:nblocks)
           enddo
        CASE ('latt_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = latt_bounds(ivertex,:,:,1:nblocks)
           enddo
        CASE ('lonu_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = lonu_bounds(ivertex,:,:,1:nblocks)
           enddo
        CASE ('latu_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = latu_bounds(ivertex,:,:,1:nblocks)
           enddo
        END SELECT

          status = pio_inq_varid(File, var_nverts(i)%short_name, varid)
          call pio_write_darray(File, varid, iodesc3dv, &
                                workr3v, status, fillval=spval)
      enddo
      deallocate(workr3v)
      endif  ! f_bounds

      deallocate(workr2)

end subroutine write_grid_variables


subroutine write_2d_variables(ns, File, i_time)

   integer, intent(in) :: ns, i_time
   type(file_desc_t), intent(inout)     :: File

   real (kind=real_kind), allocatable :: workr2(:,:,:)

   integer :: n, status
   type(var_desc_t)      :: varid
   integer (kind=PIO_OFFSET_KIND) :: FRAME_1

      FRAME_1 = i_time

      allocate(workr2(nx_block,ny_block,nblocks))

      do n=1,num_avail_hist_fields_2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            workr2(:,:,:) = a2D(:,:,n,1:nblocks)
            call pio_setframe(File, varid, FRAME_1)
            call pio_write_darray(File, varid, iodesc2d,&
                                  workr2, status, fillval=spval)
         endif
      enddo ! num_avail_hist_fields_2D

      deallocate(workr2)

end subroutine write_2d_variables


subroutine write_3d_and_4d_variables(ns, File, i_time)

   integer, intent(in) :: ns, i_time
   type(file_desc_t), intent(inout)   :: File

   real (kind=dbl_kind),  dimension(:,:), allocatable :: work_g1
   real (kind=real_kind), dimension(:,:), allocatable :: work_gr
   real (kind=real_kind), allocatable :: workr3(:,:,:,:)
   real (kind=real_kind), allocatable :: workr4(:,:,:,:,:)

   type(var_desc_t)      :: varid
   integer :: n, nn, i, j, k, ic, status
   integer (kind=PIO_OFFSET_KIND) :: FRAME_1

      FRAME_1 = i_time

          ! 3D (category)
      allocate(workr3(nx_block,ny_block,nblocks,ncat_hist))
      do n = n2D + 1, n3Dccum
         nn = n - n2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
               workr3(:,:,j,i) = a3Dc(:,:,i,nn,j)
            enddo
            enddo
            call pio_setframe(File, varid, FRAME_1)
            call pio_write_darray(File, varid, iodesc3dc,&
                                  workr3, status, fillval=spval)
         endif
      enddo ! num_avail_hist_fields_3Dc
      deallocate(workr3)

    work_gr(:,:) = c0
    work_g1(:,:) = c0

      ! 3D (vertical ice)
      allocate(workr3(nx_block,ny_block,nblocks,nzilyr))
      do n = n3Dccum+1, n3Dzcum
         nn = n - n3Dccum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, nzilyr
               workr3(:,:,j,i) = a3Dz(:,:,i,nn,j)
            enddo
            enddo
            call pio_setframe(File, varid, FRAME_1)
            call pio_write_darray(File, varid, iodesc3di,&
                                  workr3, status, fillval=spval)
         endif
      enddo ! num_avail_hist_fields_3Dz
      deallocate(workr3)

      ! 3D (vertical ice biology)
      allocate(workr3(nx_block,ny_block,nblocks,nzlyrb))
      do n = n3Dzcum+1, n3Dbcum
         nn = n - n3Dzcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, nzlyrb
               workr3(:,:,j,i) = a3Db(:,:,i,nn,j)
            enddo
            enddo
            call pio_setframe(File, varid, FRAME_1)
            call pio_write_darray(File, varid, iodesc3db,&
                                  workr3, status, fillval=spval)
         endif
      enddo ! num_avail_hist_fields_3Db
      deallocate(workr3)

      allocate(workr4(nx_block,ny_block,nblocks,ncat_hist,nzilyr))
      ! 4D (categories, vertical ice)
      do n = n3Dbcum+1, n4Dicum
         nn = n - n3Dbcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nzilyr
               workr4(:,:,j,i,k) = a4Di(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
            call pio_setframe(File, varid, FRAME_1)
            call pio_write_darray(File, varid, iodesc4di,&
                                  workr4, status, fillval=spval)
         endif
      enddo ! num_avail_hist_fields_4Di
      deallocate(workr4)

      allocate(workr4(nx_block,ny_block,nblocks,ncat_hist,nzslyr))
      ! 4D (categories, vertical ice)
      do n = n4Dicum+1, n4Dscum
         nn = n - n4Dicum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nzslyr
               workr4(:,:,j,i,k) = a4Ds(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
            call pio_setframe(File, varid, FRAME_1)
            call pio_write_darray(File, varid, iodesc4ds,&
                                  workr4, status, fillval=spval)
         endif
      enddo ! num_avail_hist_fields_4Ds
      deallocate(workr4)

end subroutine write_3d_and_4d_variables


!=======================================================================

      end module ice_history_write

!=======================================================================
