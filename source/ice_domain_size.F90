!  SVN:$Id: ice_domain_size.F90 700 2013-08-15 19:17:39Z eclare $
!=======================================================================

! Defines the global domain size and number of categories and layers.
! Code originally based on domain_size.F in POP
!
! author Elizabeth C. Hunke, LANL
! 2004: Block structure and snow parameters added by William Lipscomb
!       Renamed (used to be ice_model_size)
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!       Removed hardwired sizes (NX...can now be set in compile scripts)

      module ice_domain_size

      use ice_kinds_mod

!=======================================================================

      implicit none
      private
      save

      integer (kind=int_kind), parameter, public :: &
        ncat      = NICECAT   , & ! number of categories
        nilyr     = NICELYR   , & ! number of ice layers per category
        nslyr     = NSNWLYR   , & ! number of snow layers per category

        max_aero  =   6       , & ! maximum number of aerosols 
        n_aero    = NTRAERO   , & ! number of aerosols in use

        nblyr     = NBGCLYR   , & ! number of bio/brine layers per category
        max_nbtrcr=   9       , & ! maximum number of biology tracers
!        nltrcr    = max_nbtrcr*TRBRI, & ! maximum layer bgc tracers (for zbgc)

        max_ntrcr =   1         & ! 1 = surface temperature              
                  + nilyr       & ! ice salinity
                  + nilyr       & ! ice enthalpy
                  + nslyr       & ! snow enthalpy
                              !!!!! optional tracers:
                  + TRAGE       & ! age
                  + TRFY        & ! first-year area
                  + TRLVL*2     & ! level/deformed ice
                  + TRPND*3     & ! ponds
                  + n_aero*4    & ! number of aerosols * 4 aero layers
                  + TRBRI       & ! brine height
                  + TRBGCS    , & ! skeletal layer BGC
!                  + TRBGCZ*nltrcr*nblyr ! for zbgc (off if TRBRI=0)
        max_nstrm =   5           ! max number of history output streams

   !*** The global grid size and the block decomposition are both set at run
   !*** time through the domain_nml namelist group; see init_domain_blocks
   !*** in ice_domain.F90.  The NXGLOB/NYGLOB and BLCKX/BLCKY/MXBLCKS CPP
   !*** macros are retained only as the defaults applied when the
   !*** corresponding namelist entry is absent, so that existing builds and
   !*** namelists keep their previous behaviour.
   !***
   !*** A max_blocks higher than necessary will not cause the code to
   !*** fail, but will allocate more memory than is necessary.  A value
   !*** that is too low will cause the code to exit.  Setting
   !*** max_blocks = -1 asks the model to derive it as
   !*** max_blocks = (nx_global/block_size_x)*(ny_global/block_size_y)/
   !***               num_procs

#ifndef NXGLOB
#define NXGLOB -1
#endif
#ifndef NYGLOB
#define NYGLOB -1
#endif
#ifndef BLCKX
#define BLCKX -1
#endif
#ifndef BLCKY
#define BLCKY -1
#endif
#ifndef MXBLCKS
#define MXBLCKS -1
#endif

      integer (kind=int_kind), parameter, public :: &
        default_nx_global    = NXGLOB , & ! compile-time default, nx_global
        default_ny_global    = NYGLOB , & ! compile-time default, ny_global
        default_block_size_x = BLCKX  , & ! compile-time default, block_size_x
        default_block_size_y = BLCKY  , & ! compile-time default, block_size_y
        default_max_blocks   = MXBLCKS    ! compile-time default, max_blocks

   !*** Set from domain_nml at run time, defaulting to the values above.
   !*** These are NOT parameters: every array dimensioned by them must be
   !*** allocatable and allocated by the relevant alloc_* routine after
   !*** init_grid1.

      integer (kind=int_kind), public :: &
        nx_global    = default_nx_global    , & ! i-axis size
        ny_global    = default_ny_global    , & ! j-axis size
        block_size_x = default_block_size_x , & ! block size, first horiz dimension
        block_size_y = default_block_size_y , & ! block size, second horiz dimension
        max_blocks   = default_max_blocks       ! max number of blocks per processor

!=======================================================================

      end module ice_domain_size

!=======================================================================
