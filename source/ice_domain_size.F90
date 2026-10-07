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

   !*** The global grid size and the block decomposition are set at run time
   !*** from the domain_nml namelist group; see init_domain_blocks in
   !*** ice_domain.F90.  nx_global, ny_global, block_size_x and block_size_y
   !*** must all be given there.
   !***
   !*** max_blocks may be given, or left at -1 to have the model derive it
   !*** from the block distribution it actually builds.  A max_blocks larger
   !*** than necessary is not fatal, but wastes memory; one that is too small
   !*** aborts the run.
   !***
   !*** These are NOT parameters: every array dimensioned by them must be
   !*** allocatable, and allocated by the relevant alloc_* routine after
   !*** init_grid1.

      integer (kind=int_kind), public :: &
        nx_global    = -1 , & ! i-axis size
        ny_global    = -1 , & ! j-axis size
        block_size_x = -1 , & ! block size, first horizontal dimension
        block_size_y = -1 , & ! block size, second horizontal dimension
        max_blocks   = -1     ! max number of blocks per processor

!=======================================================================

      end module ice_domain_size

!=======================================================================
