# Standalone CICE driver

The `cice` driver builds an uncoupled executable for exercising the sea-ice
model on its own. It requires neither OASIS nor libaccessom2, and defines
neither `AusCOM` nor `coupled`:

```bash
cmake -S . -B build \
      -DCICE_DRIVER=cice -DCICE_IO=NetCDF -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

See `<CICE_ROOT>/input_templates/run_ice.gadi.nci.org.au` for a Gadi PBS script that runs
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
