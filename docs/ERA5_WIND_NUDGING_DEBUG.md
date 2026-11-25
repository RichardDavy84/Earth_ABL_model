# ERA5 wind nudging: diagnostic branch

This branch addresses problems in the fixed-ERA5 wind test; it is **not** a
validated production ABL physics change.

## Patch

1. **Exclude zero-area surface tiles from the nudging boundary-layer height.**
   After the category integrations, `ABL.f90` computes the shared previous-step
   `blht_max` from tiles whose `sic > 1e-8` and diagnosed `blht > 0`.
   If none qualifies, a sentinel `-1` triggers the existing fallback.
   This prevents a zero-area open-water tile from controlling ice nudging.
   It does **not** yet implement independent per-tile nudging heights.
2. **Refresh `ug`, `vg`, and `z0` for every integration call.**
   `Command_module_for_NeXtSIM.for` now assigns `ug=u_hPa(7)`,
   `vg=v_hPa(7)`, and `z0=z0_in`, rather than carrying module values
   from an earlier initialized column. Index 7 is the 850 hPa level.
   Caution: ERA5 wind at 850 hPa is *not* necessarily geostrophic wind.
3. **Prevent cross-thread interference in the column solver.**
   The call to `Integrate_NeXtSIM_ABL` is surrounded by an OpenMP
   `CRITICAL(ABL_PHYSICS_STATE)` region because `physics_vars` is shared.
   This is a temporary diagnostic protection, not a performance fix;
   column integrations are now serialized.
4. **Compare the correct ERA5 target with area-weighted model winds.**
   Each tile's pressure-interpolated `u_Pa_to_z` / `v_Pa_to_z` is
   retained in a per-thread scratch array during integration. The
   saved `u_from_ERA` and `v_from_ERA` are area-weighted using the same
   surface fractions as the output model `u` and `v`. Before this
   change the reference arrays held only the last processed tile.
5. **Use intended nudging timescales in `setup_info.nml`.**
   `nudge_above_bl=1.5` hours (90 min), `nudge_in_bl=5.0` hours.
   The implementation expects timescales in hours, not coefficients.
   At a 10-second step, `nudge_fac` is about 0.001852 above the BL,
   0.000556 inside, with a three-point transition.
6. **Initialize BLH if the stress threshold search has no crossing.**
   The two BL-height calculations use `blht=-1.` as an invalid sentinel
   rather than leaving BLH undefined.

## Suggested experiments

Keep the ERA5 forcing constant (`repeat_forcing=-1`), `nudge_level=0`,
and the grid and integration step unchanged.

1. Run `OMP_NUM_THREADS=1`, `do_tiling=1` with 100% sea ice in
   the first category, and uniform 30-min nudging (`0.5,0.5` hours).
2. Repeat with uniform 90-min nudging (`1.5,1.5` hours).
3. Repeat with 90 min above / 5 h below (`1.5,5.0` hours).
4. Repeat with nonzero fractions in multiple ice/ocean categories;
   compare area-averaged U and V separately to saved
   `u_from_ERA` and `v_from_ERA`, only at pressure levels where
   ERA5 1000–700 hPa forcing is available.
5. Compare `OMP_NUM_THREADS=1` and `4`. The integrator should now
   avoid cross-thread module-state races; it will not scale efficiently.

For the `do_si_ics=3` test mode, verify actual SIC at every output time:
its special experiments may change SIC, independently of this patch.
With `do_tiling=1`, **zero-area categories are still integrated**,
even though they no longer affect grid-mean output or the active-tile BLH
maximum. `do_tiling=0` runs a single surface column.

## Known limitations and next work

- There is no guarantee that the ERA5 850 hPa wind is the correct
  *geostrophic* wind for the model momentum equations.
- The model's upper boundary imposes zero vertical shear and the ERA5
  profile reaches only 700 hPa. Compare only nudged pressure levels.
- The shared maximum active-tile BLH is an interim algorithm;
  per-active-tile previous-step BL heights would be more consistent.
- OpenMP protection around the integrator is a correctness diagnostic.
  Other shared state in the main driver may still require an audit.
- A failed BL-height diagnosis uses the existing initial-threshold
  fallback; this needs explicit validation in model runs.
- This branch has **not been compiled or tested on the cluster**.
  Compile with runtime bounds/invalid checks where possible before
  interpreting simulated results.
