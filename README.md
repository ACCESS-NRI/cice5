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

The `cice` driver builds an uncoupled executable for exercising the sea-ice
model on its own. It requires neither OASIS nor libaccessom2, and defines
neither `AusCOM` nor `coupled`:

```bash
cmake -S . -B build \
      -DCICE_DRIVER=cice -DCICE_IO=NetCDF -DCMAKE_BUILD_TYPE=Release \
      -DCICE_NXGLOB=100 -DCICE_NYGLOB=116
cmake --build build -j
```

See `input_templates/run_ice.gadi.nci.org.au` for a Gadi PBS script that runs
it. Note that the standalone driver is not part of any ACCESS configuration and
is not covered by the CI, so treat it as a development and testing aid.

For a smoke test that needs **no input data at all**, set `grid_type =
'rectangular'` in `grid_nml` and `atm_data_type = 'default'` in `forcing_nml`.
The grid is then generated analytically and the forcing is synthetic, so the
model runs from the namelist alone. This is a convenient way to check that a
build works and that results are independent of the block decomposition: run
the same executable at several `block_size_x`/`block_size_y`/`nprocs`
combinations and compare the restart files, which should be bitwise identical.
(The history field `blkmask` is expected to differ — it records which task and
block each cell belongs to.)

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
  , max_blocks   = -1        ! -1 => derive from the above and nprocs
/
```

Notes:

* `block_size_x` and `block_size_y` need not divide `nx_global`/`ny_global`
  evenly — the decomposition is padded — but choosing sizes that do avoids
  wasted work.
* `max_blocks = -1` derives
  `((nx_global-1)/block_size_x + 1) * ((ny_global-1)/block_size_y + 1) / nprocs`.
  If the chosen distribution then needs more blocks on some task than that, the
  model aborts and reports the value required; set `max_blocks` explicitly in
  that case. A value larger than necessary only wastes memory.
* The values in effect are echoed to the ice diagnostic log under
  `Domain Information`.
* **The ESM (`access`) driver requires `max_blocks == 1`**, i.e. exactly one
  block per task, because its atmosphere coupling unpacks directly into block
  index 1. It aborts at startup otherwise. The OM2 (`auscom`) driver has no
  such restriction.

`CICE_NXGLOB`, `CICE_NYGLOB`, `CICE_BLCKX`, `CICE_BLCKY` and `CICE_MXBLCKS`
are retained in the CMake build, but only as the **defaults** applied when the
corresponding `domain_nml` entry is absent, so existing namelists keep their
previous behaviour. Set any of them to `-1` to compile in no default and make
that namelist entry mandatory. Since nothing about the domain is baked into
the executable any more, it is simply `cice_<driver>.exe`.

Note that the vertical and tracer dimensions (`NICECAT`, `NICELYR`,
`NSNWLYR`, the tracer counts) are still compile-time.

## Useful links

* **Wiki**: https://github.com/CICE-Consortium/CICE-svn-trunk/wiki

   Information about the CICE model prior to version 6 including how to obtain the code

* **Version Index**: https://github.com/CICE-Consortium/CICE-svn-trunk/wiki/CICE-Versions-Index-(older)

   Numbered CICE releases prior to version 6. 

* **Resource Index**: https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index

   List of resources for information about the Consortium and its repositories as well as model documentation, testing, and development.
