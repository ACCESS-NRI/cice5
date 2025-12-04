# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

cmake_minimum_required(VERSION 3.18)

#[==============================================================================[
#                         Basic project definition                              #
#]==============================================================================]

project(CICE 
        DESCRIPTION "CICE5 Sea Ice Model"
        HOMEPAGE_URL https://github.com/ACCESS-NRI/cice5
        LANGUAGES Fortran)

#[==============================================================================[
#                                   Options                                     #
#]==============================================================================]

message(STATUS "Build options")

# provide some default options from access "1degree" configurations
# its intended the user will set these

set(CICE_IO "PIO" CACHE STRING "CICE netcdf io method")

if(NOT CICE_IO MATCHES "^(NetCDF|PIO)$")
  message(FATAL_ERROR " CICE_IO   ${CICE_IO} not valid, choose NetCDF|PIO")
else()
  message(STATUS "  - CICE_IO   ${CICE_IO}")
endif()

set(CICE_DRIVER "auscom" CACHE STRING "CICE driver code to use")
message(STATUS "  - CICE_DRIVER   ${CICE_DRIVER}")

set(CICE_NXGLOB "360" CACHE STRING "Size of grid in X-direction")
set(CICE_NYGLOB "300" CACHE STRING "Size of grid in Y-direction")

# defaults from https://github.com/ACCESS-NRI/access-spack-packages/blob/d8b48b18a1400cf2a6850af0bf038af55804d92f/packages/cice5/spack-build.sh#L89-L91
# BLCKX = 360/24 = 15 , BLCKY = 300/1  = 300
set(CICE_BLCKX   "15"  CACHE STRING "Size of computational blocks in X-direction")
set(CICE_BLCKY   "300" CACHE STRING "Size of computational blocks in Y-direction")
set(CICE_MXBLCKS "1"   CACHE STRING "Maximum number of blocks per task")

# these are also configuration settings, but we leave them as internal only at this point
# could be made options if users need them

if (CICE_DRIVER MATCHES "access")
  set( NSNWLYR  1 )
  set( NICECAT  5 )
  set( NICELYR  1 )
else()
  set( NSNWLYR  1 )       # number of vertical layers in the snow
  set( NICECAT  5 )       # number of ice thickness categories
  set( NICELYR  4 )       # number of vertical layers in the ice
endif()

### Tracers               # match ice_in tracer_nml to conserve memory
set( TRBGCS  0 )       # number of skeletal layer bgc tracers 
set( TRAGE   1 )       # set to 1 for ice age tracer
set( TRFY    1 )       # set to 1 for first-year ice area tracer
set( TRLVL   1 )       # set to 1 for level and deformed ice tracers
set( TRPND   1 )       # set to 1 for melt pond tracers
set( NTRAERO 0 )       # number of aerosol tracers 
                          # (up to max_aero in ice_domain_size.F90) 
                          # CESM uses 3 aerosol tracers
set( TRBRI   1 )       # set to 1 for brine height tracer
set( NBGCLYR 0 )       # number of zbgc layers

### File unit numbers
set( NUMIN 11  )      # minimum file unit number
set( NUMAX 99  )      # maximum file unit number

#[==============================================================================[
#                           Project configuration                               #
#]==============================================================================]

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(CMakePackageConfigHelpers)

# https://github.com/ACCESS-NRI/access-spack-packages/blob/d8b48b18a1400cf2a6850af0bf038af55804d92f/packages/cice5/package.py#L197C1-L198C65
# note: previous deterministic and optimisation report options not included
set (NCI_INTEL_FLAGS "-r8 -i4 -w -fpe0 -ftz -convert big_endian -assume byterecl -check noarg_temp_created")
set (NCI_REPRO_FLAGS "-fp-model precise -fp-model source -align all")
set (NCI_OPTIM_FLAGS "-g3 -O2 -axCORE-AVX2 -debug all -check none -traceback -assume buffered_io")
set (NCI_DEBUG_FLAGS "-g3 -O0 -debug all -check all -no-vec -traceback -assume nobuffered_io")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Wall -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "${NCI_INTEL_FLAGS} ${NCI_REPRO_FLAGS} ${NCI_OPTIM_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${NCI_INTEL_FLAGS} ${NCI_REPRO_FLAGS} ${NCI_DEBUG_FLAGS}")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

add_compile_definitions(
  NXGLOB=${CICE_NXGLOB}
  NYGLOB=${CICE_NYGLOB}
  BLCKX=${CICE_BLCKX}
  BLCKY=${CICE_BLCKY}
  MXBLCKS=${CICE_MXBLCKS}
  NUMIN=${NUMIN}
  NUMAX=${NUMAX}
  TRAGE=${TRAGE}
  TRFY=${TRFY}
  TRLVL=${TRLVL}
  TRPND=${TRPND}
  NTRAERO=${NTRAERO}
  TRBRI=${TRBRI}
  NBGCLYR=${NBGCLYR}
  TRBGCS=${TRBGCS}
  NICECAT=${NICECAT}
  NICELYR=${NICELYR}
  NSNWLYR=${NSNWLYR}
)

add_compile_definitions( gather_scatter_barrier ncdf OASIS3_MCT )

if (CICE_IO MATCHES PIO)
  add_compile_definitions(PIO)
endif()

if(CICE_DRIVER MATCHES "auscom")
  add_compile_definitions(AusCOM coupled)
elseif(CICE_DRIVER MATCHES "access")
  add_compile_defitions(ACCESS AusCOM coupled)
endif()

#[==============================================================================[
#                              External packages                                #
#]==============================================================================]

find_package(PkgConfig REQUIRED)
pkg_check_modules(OASIS3PSMILE REQUIRED IMPORTED_TARGET "oasis3-psmile.MPI1")
pkg_check_modules(OASIS3MCT REQUIRED IMPORTED_TARGET "oasis3-mct")
pkg_check_modules(OASIS3MPEU REQUIRED IMPORTED_TARGET "oasis3-mpeu")
pkg_check_modules(OASIS3SCRIP REQUIRED IMPORTED_TARGET "oasis3-scrip")

find_package(MPI REQUIRED COMPONENTS Fortran)
if(CICE_DRIVER MATCHES "auscom")
  find_package(LIBACCESSOM2)
endif()

if(CICE_IO MATCHES "PIO" AND NOT TARGET PIO::PIO_Fortran)
  # Code has not been tested with versions older than 2.5.3, but probably still works fine
  # For multiple timesteps per history file, needs latest from main (>=2.6.7)
  find_package(PIO 2.5.3 REQUIRED COMPONENTS Fortran)
endif()
if(CICE_IO MATCHES "^(NetCDF|PIO)$" AND NOT TARGET NetCDF::NetCDF_Fortran)
  # Code has not been tested with versions older than 4.7.3, but probably still works fine 
  find_package(NetCDF 4.7.3 REQUIRED Fortran)
endif()

#[==============================================================================[
#                             Main definitions                                  #
#]==============================================================================]

set(CICE_SOURCE "${CMAKE_SOURCE_DIR}/source")
set(CICE_MPI "${CMAKE_SOURCE_DIR}/mpi")
set(CICE_DRIVER_SOURCE "${CMAKE_SOURCE_DIR}/drivers/${CICE_DRIVER}")

add_executable(cice    ${CICE_DRIVER_SOURCE}/CICE.F90)

target_link_libraries(cice 
  PUBLIC libaccessom2::accessom2
  PRIVATE 
  MPI::MPI_Fortran 
  PkgConfig::OASIS3PSMILE PkgConfig::OASIS3MCT PkgConfig::OASIS3MPEU PkgConfig::OASIS3SCRIP 
  # OASIS3MCT is a depedency of OASIS3PSMILE, so order matters here
  )

if(CICE_IO MATCHES "^(NetCDF|PIO)$")
  target_link_libraries(cice PRIVATE NetCDF::NetCDF_Fortran NetCDF::NetCDF_C)
  if(CICE_IO MATCHES "PIO")
    target_link_libraries(cice PRIVATE PIO::PIO_Fortran)
  endif()
endif()

target_sources(cice PRIVATE
  ${CICE_SOURCE}/ice_aerosol.F90
  ${CICE_SOURCE}/ice_dyn_evp.F90
  ${CICE_SOURCE}/ice_meltpond_lvl.F90
  ${CICE_SOURCE}/ice_step_mod.F90
  ${CICE_SOURCE}/ice_age.F90
  ${CICE_SOURCE}/ice_dyn_shared.F90
  ${CICE_SOURCE}/ice_history_mechred.F90
  ${CICE_SOURCE}/ice_meltpond_topo.F90
  ${CICE_SOURCE}/ice_therm_0layer.F90
  ${CICE_SOURCE}/ice_algae.F90
  ${CICE_SOURCE}/ice_fileunits.F90
  ${CICE_SOURCE}/ice_history_pond.F90
  ${CICE_SOURCE}/ice_ocean.F90
  ${CICE_SOURCE}/ice_therm_bl99.F90
  ${CICE_SOURCE}/ice_atmo.F90
  ${CICE_SOURCE}/ice_history_shared.F90
  ${CICE_SOURCE}/ice_orbital.F90
  ${CICE_SOURCE}/ice_therm_itd.F90
  ${CICE_SOURCE}/ice_blocks.F90
  ${CICE_SOURCE}/ice_read_write.F90
  ${CICE_SOURCE}/ice_therm_mushy.F90
  ${CICE_SOURCE}/ice_brine.F90
  ${CICE_SOURCE}/ice_firstyear.F90
  ${CICE_SOURCE}/ice_init.F90
  ${CICE_SOURCE}/ice_restart_driver.F90
  ${CICE_SOURCE}/ice_therm_shared.F90
  ${CICE_SOURCE}/ice_calendar.F90
  ${CICE_SOURCE}/ice_flux.F90
  ${CICE_SOURCE}/ice_itd.F90
  ${CICE_SOURCE}/ice_restart_shared.F90
  ${CICE_SOURCE}/ice_therm_vertical.F90
  ${CICE_SOURCE}/ice_diagnostics.F90
  ${CICE_SOURCE}/ice_forcing.F90
  ${CICE_SOURCE}/ice_kinds_mod.F90
  ${CICE_SOURCE}/ice_restart_shared.F90
  ${CICE_SOURCE}/ice_transport_driver.F90
  ${CICE_SOURCE}/ice_distribution.F90
  ${CICE_SOURCE}/ice_grid.F90
  ${CICE_SOURCE}/ice_restoring.F90
  ${CICE_SOURCE}/ice_transport_remap.F90
  ${CICE_SOURCE}/ice_domain.F90
  ${CICE_SOURCE}/ice_history_bgc.F90
  ${CICE_SOURCE}/ice_lvl.F90
  ${CICE_SOURCE}/ice_shortwave.F90
  ${CICE_SOURCE}/ice_zbgc.F90
  ${CICE_SOURCE}/ice_domain_size.F90
  ${CICE_SOURCE}/ice_history_drag.F90
  ${CICE_SOURCE}/ice_mechred.F90
  ${CICE_SOURCE}/ice_spacecurve.F90
  ${CICE_SOURCE}/ice_zbgc_shared.F90
  ${CICE_SOURCE}/ice_dyn_eap.F90
  ${CICE_SOURCE}/ice_history.F90
  ${CICE_SOURCE}/ice_meltpond_cesm.F90
  ${CICE_SOURCE}/ice_state.F90
  ${CMAKE_SOURCE_DIR}/csm_share/shr_orb_mod.F90
)

target_sources(cice PRIVATE
  ${CICE_MPI}/ice_boundary.F90
  ${CICE_MPI}/ice_communicate.F90
  ${CICE_MPI}/ice_gather_scatter.F90
  ${CICE_MPI}/ice_global_reductions.F90
  ${CICE_MPI}/ice_broadcast.F90
  ${CICE_MPI}/ice_exit.F90
  ${CICE_MPI}/ice_timers.F90
)

if(CICE_DRIVER MATCHES "access")
  target_sources(cice PRIVATE
    ${CICE_DRIVER_SOURCE}/CICE_RunMod.F90
    ${CICE_DRIVER_SOURCE}/cpl_forcing_handler.F90
    ${CICE_DRIVER_SOURCE}/cpl_parameters.F90
    ${CICE_DRIVER_SOURCE}/CICE_FinalMod.F90
    ${CICE_DRIVER_SOURCE}/cpl_interface.F90
    ${CICE_DRIVER_SOURCE}/ice_constants.F90
    ${CICE_DRIVER_SOURCE}/CICE_InitMod.F90
    ${CICE_DRIVER_SOURCE}/cpl_arrays_setup.F90
    ${CICE_DRIVER_SOURCE}/cpl_netcdf_setup.F90
  )
elseif(CICE_DRIVER MATCHES "auscom")
  target_sources(cice PRIVATE
    ${CICE_DRIVER_SOURCE}/cpl_forcing_handler.F90
    ${CICE_DRIVER_SOURCE}/sat_vapor_pres_k_mod.F90
    ${CICE_DRIVER_SOURCE}/CICE_FinalMod.F90
    ${CICE_DRIVER_SOURCE}/cpl_interface.F90
    ${CICE_DRIVER_SOURCE}/monin_obukhov_kernel.F90
    ${CICE_DRIVER_SOURCE}/sat_vapor_pres_mod.F90
    ${CICE_DRIVER_SOURCE}/CICE_InitMod.F90
    ${CICE_DRIVER_SOURCE}/cpl_netcdf_setup.F90
    ${CICE_DRIVER_SOURCE}/monin_obukhov_mod.F90
    ${CICE_DRIVER_SOURCE}/surface_flux_mod.F90
    ${CICE_DRIVER_SOURCE}/CICE_RunMod.F90
    ${CICE_DRIVER_SOURCE}/cpl_parameters.F90
    ${CICE_DRIVER_SOURCE}/ocean_rough_mod.F90
    ${CICE_DRIVER_SOURCE}/cpl_arrays_setup.F90
    ${CICE_DRIVER_SOURCE}/ice_constants.F90
  )
else()
  message(FATAL_ERROR "CICE_DRIVER: ${CICE_DRIVER} is currently not supported by the CMake build system")
endif()

# Select IO source files based on CICE_IO
if(CICE_IO MATCHES "NetCDF")
  target_sources(cice PRIVATE
    ${CMAKE_SOURCE_DIR}/io_netcdf/ice_history_write.F90
    ${CMAKE_SOURCE_DIR}/io_netcdf/ice_restart.F90
  )
elseif(CICE_IO MATCHES "PIO")
  target_sources(cice PRIVATE
    ${CMAKE_SOURCE_DIR}/io_pio/ice_history_write.F90
    ${CMAKE_SOURCE_DIR}/io_pio/ice_pio.F90
    ${CMAKE_SOURCE_DIR}/io_pio/ice_restart.F90
  )
endif()

#[==============================================================================[
#                             Install                                           #
#]==============================================================================]

set_target_properties(cice PROPERTIES
  LINKER_LANGUAGE Fortran
  OUTPUT_NAME "cice_${CICE_DRIVER}_${CICE_NXGLOB}x${CICE_NYGLOB}_${CICE_BLCKX}x${CICE_BLCKY}.exe"
)

install(TARGETS cice
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT IO_${CICE_IO}
)