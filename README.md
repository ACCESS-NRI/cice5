## Overview
This repository contains the access/auscom fork of cice5 used in the ACCESS-ESM1.6 and ACCESS-OM2. It was forked from https://github.com/CICE-Consortium/CICE-svn-trunk/ which in turn captured the trunk from the subversion (svn) repository of the Los Alamos Sea Ice Model, CICE, including release tags through version 5.1.2.

ACCESS-ESM1.6 related code came from https://code.metoffice.gov.uk/trac/cice/browser/main/branches/pkg/Config/vn5.1.2_GSI8.1_package_branch/cice?order=name&rev=334#source (MOSRS account required)

More recent versions are found in the [CICE](https://github.com/CICE-Consortium/CICE) and [Icepack](https://github.com/CICE-Consortium/Icepack) repositories, which are maintained by the CICE Consortium.

There is [PDF documentation](https://github.com/ACCESS-NRI/cice5/blob/master/doc/cicedoc.pdf)([DOI](https://doi.org/10.5281/zenodo.19207490)) available for CICE 5.1.2, however some changes were made to this fork to support coupling with ACCESS-OM2 and ACCESS-ESM1.6, Parallel IO, ERA5 Forcing, BGC modelling and for other updates. Some of these changes are described in the [ACCESS-OM2 Technical Report](https://github.com/COSIMA/ACCESS-OM2-1-025-010deg-report).


### Build Systems

We recommend and support using CMake to build cice5, however the Makefile build is kept in this repository for any legacy models still using these files.

Three drivers are supported by the CMake build, selected with `CICE_DRIVER`:

| `CICE_DRIVER` | Used by | Coupling |
| --- | --- | --- |
| `auscom` (default) | ACCESS-OM2 | OASIS3-MCT + libaccessom2 |
| `access` | ACCESS-ESM1.6 | OASIS3-MCT |
| `cice` | standalone / testing | none |


### Grid size, block decomposition and task count

The global grid size and the block decomposition are both chosen at **run
time**, in the `domain_nml` group of the namelist (`cice_in.nml` for the
`auscom` and `access` drivers, `ice_in` otherwise). A single executable can
therefore be run at any resolution and any number of MPI tasks, without
rebuilding:

```fortran
&domain_nml
    nprocs       = 24        ! must equal the MPI task count
  , nx_global    = 360       ! global grid size, i
  , ny_global    = 300       ! global grid size, j
  , block_size_x = 15        ! block size in i, excluding ghost cells
  , block_size_y = 300       ! block size in j, excluding ghost cells
  , max_blocks   = -1        ! -1 => derive from the actual distribution
/
```

Notes:

This runtime specification for domain grids, `nprocs`,
block sizes, and max blocks, was backported from CICE6 in this
[CICE5 PR](https://github.com/ACCESS-NRI/cice5/pull/113); please refer to
[CICE6 docs](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_case_settings.html)
for details on the namelist parameters.


* `block_size_x` and `block_size_y` need not divide `nx_global`/`ny_global`
  evenly — the decomposition is padded — but choosing sizes that do avoids
  wasted work.
* `max_blocks = -1` is the recommended setting. Each task sets it to the
  number of blocks that task actually owns, once the distribution has been
  built, so it is exact and needs no guessing. It is a **per-task** value: an
  uneven distribution gives different tasks different `max_blocks`, and each
  allocates only what it needs. The minimum and maximum across tasks are
  reported to the diagnostic log.
  Setting it explicitly still works, and still aborts if a task turns out to
  need more blocks than that; a value larger than necessary only wastes
  memory.
* There is nothing to tune for the `roundrobin`, `sectrobin` and `sectcart`
  distributions. They build an internal `blockIndex` work array before the
  per-task block count is known, so it is sized from an estimate and grown
  on demand if that estimate is too small; no run can fail because of it.
  `cartesian`, `rake` and `spacecurve` do not use that array at all.
* The values in effect are echoed to the ice diagnostic log under
  `Domain Information`.
* **The ESM (`access`) driver requires `max_blocks == 1`**, i.e. exactly one
  block per task, because its atmosphere coupling unpacks directly into block
  index 1. It aborts at startup otherwise. The OM2 (`auscom`) driver has no
  such restriction.

`nx_global`, `ny_global`, `block_size_x` and `block_size_y` are **required**
in `domain_nml`; `max_blocks` may be given or left at `-1` to be derived. The
build has no grid or block settings at all — the old `-DNXGLOB`, `-DNYGLOB`,
`-DBLCKX`, `-DBLCKY` and `-DMXBLCKS` macros and their `CICE_*` CMake variables
are gone. Since nothing about the domain is baked into the executable, it is
simply `cice_<driver>.exe`.

**Upgrading an existing namelist:** a `domain_nml` that relied on the
compiled-in values must now spell out `nx_global`, `ny_global`,
`block_size_x` and `block_size_y`. Use the values the build previously passed
via `-D`; the run aborts at startup with a message naming the missing setting
if any is absent.

Note that the vertical and tracer dimensions (`NICECAT`, `NICELYR`,
`NSNWLYR`, the tracer counts) are still compile-time.

## Useful links

* **Wiki**: https://github.com/CICE-Consortium/CICE-svn-trunk/wiki

   Information about the CICE model prior to version 6 including how to obtain the code

* **Version Index**: https://github.com/CICE-Consortium/CICE-svn-trunk/wiki/CICE-Versions-Index-(older)

   Numbered CICE releases prior to version 6. 

* **Resource Index**: https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index

   List of resources for information about the Consortium and its repositories as well as model documentation, testing, and development.
