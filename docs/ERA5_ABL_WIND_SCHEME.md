# ERA5 winds above the tallest active ABL: diagnostic implementation

This branch is an experimental change to the Earth_ABL_model. It is intended
for controlled comparisons with Heather's fixed-ERA5 wind runs. It is not a
validated production momentum-forcing scheme.

## Behaviour

For every geographical cell, the model computes a **common nudging BL top**
from the previous timestep's largest valid BL height among **active**
surface categories (SIC > 1e-8). A missing BL-height diagnosis uses a
conservative no-exact-forcing fallback. Upward changes in the diagnosed BL
top are followed immediately; downward changes are smoothed with a
300-second decay time. Inactive ice/open-water tiles cannot determine
the nudging transition. The model still integrates zero-area tiles but
excludes their diagnostics from the shared BL maximum.

Within this shared BL top, the **U and V components** are independently,
weakly nudged to ERA5 with a configurable e-folding time, default 5 hours.
Above the BL top a cubic smoothstep transition changes the per-step wind
adjustment factor continuously from the weak-ABL value to 1 (exact ERA5
matching). Transition width is at least 200 m **and** spans at least
three model levels above the BL top wherever there are enough levels.
These are height-based, not discrete nudging-index jumps.

At heights where ERA5 pressure-level forcing is available (700--1000 hPa),
the winds above the transition are set to the interpolated ERA5 U/V profile
on *every timestep*, after the dynamical solver. The solver's **top wind
boundary** is also changed to an ERA5 Dirichlet condition, but only
when the model top is above the transition and has valid ERA5 coverage.
Otherwise its original zero-shear top condition is retained.

Below 1000 hPa coverage the closest 1000 hPa ERA5 wind is used only
for weak nudging, **never** exact wind prescription. Above the 700 hPa
forcing range winds are not nudged at all. The surface boundary (model
level 1) is not nudged. Thus exact matching outside the supplied ERA5
pressure range is **not** claimed.

Temperature continues to use the **existing** `nudge_above_bl`,
`nudge_in_bl` and `nudge_grad` settings and is no longer tied to
the new wind blending coefficients. Pressure interpolation for the ERA5
reference is moved after hydrostatic p(z) recomputation.

## Code

- `wind_nudging.f90`: pure, separately testable per-level wind
  blending function; same factor for U and V.
- `Command_module_for_NeXtSIM.for_20261001`: updated interpolation
  sequence; separate post-solver wind correction; top-boundary flag;
  refresh current column's ERA5 850 hPa geostrophic-wind **proxy**
  and surface roughness.
- `Physics.for`: conditional ERA5 Dirichlet U/V at the model top;
  otherwise original zero-gradient top.
- `ABL.f90`: shared active-tile BL top and conservative filtering;
  area-weighted ERA5 reference diagnostics; per-thread target scratch
  arrays; temporary OpenMP CRITICAL around the column solver to protect
  shared `physics_vars`. This serializes the expensive part of the model.
- `Makefile_olivia_gfortran`: wind_nudging module compiled before
  the integration routine; original `Makefile` symlink is preserved.
- `setup_info.nml`: optional `wind_tau_bl_hours=5.0` and
  `wind_transition_m=200.0` settings.
- `tests/test_wind_nudging.f90`: small regression test for blending
  inside, across, and above the BL, and outside ERA5 coverage.

Both original repository symlinks, `Makefile` and
`Command_module_for_NeXtSIM.for`, are preserved, and their targets
are modified instead. `master` is unchanged.

## Tests to run before scientific interpretation

Compile the isolated blending regression test:

```sh
gfortran -Wall -Wextra -fcheck=all \
  wind_nudging.f90 tests/test_wind_nudging.f90 -o /tmp/test_wind
/tmp/test_wind
```

Then rebuild the full project with `make clean && make`, using the
existing cluster libraries. Start with `OMP_NUM_THREADS=1`, fixed ERA5
forcing (`repeat_forcing=-1`), and `nudge_level=0`.

1. Run a single column and 100%-sea-ice category for a short test.
   `do_tiling=1` still integrates zero-area tiles; use `do_tiling=0`
   if exactly one column is desired. `do_si_ics=3` can change SIC
   over time; verify actual concentrations.
2. Plot U and V separately from `Met.nc` against `u_from_ERA`
   and `v_from_ERA`. ERA5 comparison fields are now weighted using
   the same active-tile concentrations as the model mean winds.
3. Plot `blht_nudge` from `SRFV-ALL.nc`, and the height where
   the transition ends. The diagnosed BL height is from the previous
   timestep, updated at the end of the current one.
4. Above that transition **and inside 700--1000 hPa**, U and V should
   match their interpolated ERA5 target to floating-point roundoff.
   Failure indicates a mask, target, solver, or output issue.
5. Below the BL top, compare wind turning, vertical shear, friction
   velocity and surface stress, not just ERA5 mismatch. The ABL is
   supposed to respond independently to its surface.
6. Test convergence and avoid oscillations near the BL top when the
   diagnosed height varies. Compare 200 m and 400 m transition widths,
   and 3 h vs 5 h in-BL weak relaxation.

## Scientific and numerical limitations

- The existing model still treats ERA5's 850 hPa U/V as a **proxy**
  for geostrophic wind. Computing pressure-gradient forcing from
  geopotential remains separate future work.
- An invalid BL diagnosis suppresses exact forcing until the
  model produces a valid BL top; inspect early-time output.
- At a high BL top or a shallow model domain, there may be no levels
  above the transition, so the exact matching region does not exist.
- There is no ERA5 reference above 700 hPa; extrapolation is avoided.
- Strictly prescribed wind at each step acts as an external momentum
  source. This is by design, and can affect stress and TKE through
  the next timestep.
- OpenMP serialization is a temporary correctness safeguard;
  shared module state needs a future per-column refactor for efficient
  multi-threaded production.
- The local stand-alone coefficient test is not a full dynamical model
  validation; cluster compilation and numerical runs are still needed.
