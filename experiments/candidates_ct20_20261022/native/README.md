# Candidate-search native ray kernel

This is a new isolated float64 numerical kernel, not the unavailable earlier
speed experiment and not a replacement of the repository's original engine.
The original `src/` files remain unchanged.

Build explicitly before importing:

```bash
python experiments/candidates_20261022/native/build.py
```

With `experiments/candidates_20261022` on `sys.path`:

```python
from native import trace
result = trace(surfaces, origins, directions, indices)
```

- `surfaces`: `(S,9)`, columns `z,c,k,A4,A6,A8,A10,A12,clear_radius`.
- `origins`, `directions`: `(N,3)`. Directions are normalized by the kernel.
- `indices`: `(N,S+1)` or `(S+1,)`, indices in successive media.
- All lengths are mm. A4 through A12 are dimensional coefficients multiplying
  `r**4` through `r**12`, not normalized sag amplitudes.
- `pos`, `dir`, `opl`: last successfully crossed boundary, retaining the initial
  ray state when the first boundary fails.
- `history`, `direction_history`: `(N,S,3)` after each successful refraction.
  `opl_history`: `(N,S)`. The failed boundary and later entries are NaN.
- `failed_surface`: zero-based index, `-1` for success.
- Status `0`: success; `1`: no hit on the finite clear surface;
  `2`: total internal reflection; `3`: unresolved/invalid intersection or
  numerical degeneracy; `4`: transmitted ray does not progress in positive z.

Every launched ray retains its position in the output. Optimization and
evaluation consumers must not discard failed rays and then pretend that the
remaining sample represents full throughput.

The kernel starts accumulated OPL at zero. A tilted plane wave launched from a
fixed plane needs its incident phase offset `initial_direction @ initial_origin`
added to the ray's OPL before coherent wavefront analysis. STOP aiming, image
projection, Fresnel losses, global geometry/nonintersection validation and
material dispersion are caller responsibilities.

The surface equation uses the standard stable conic sag plus even polynomial
orders4–12. A finite-cylinder/sag bracket is searched in increasing ray
distance. Numerical interval bounds of the surface function and its derivative
discard empty intervals and certify simple monotonic brackets; safeguarded
Newton/bisection solves those brackets. An unresolved earlier or tangent
interval fails explicitly instead of skipping to a later root. This is not a
formal interval-arithmetic proof for arbitrary extreme coefficients.

Global lens overlap and mechanical zones outside clear radii are not inferred.
Candidate acceptance must validate them separately. The implementation targets
six distinct physical refracting boundaries; coincident curved dummy STOP
chains are not a promised interface. Spherical glass restrictions are also
caller-side prescription constraints.

`validation.json` records408 comparisons to the original engine's independently
implemented interval intersection and vector Snell:353 successful rays,
51 escaped rays and4TIR, with zero status/surface-index disagreements. Maximum
position, direction and OPL differences are respectively2.613e−11mm,
1.079e−13 and4.855e−11mm. These include polynomial orders through12 and a
six-boundary triplet. `first_root_controls.json` additionally verifies first of
two positive paraboloid roots, tangent rejection and zero-direction rejection.

The small `benchmark_smoke.json` records399 rays in about1.515ms median through
the new model wrapper versus1.365s for one scalar-reference call. It contains
240 successful rays and159 optical failures, so it is neither an all-valid
throughput benchmark nor an end-to-end speedup claim.

Compilation is GCC `-O3 -fno-fast-math -ffp-contract=off`, with compiler/source
and binary hashes in `build.json`. The failed development library-name attempt
is retained as `development_raytrace.so`; it is not loaded. Do not confuse a
ctypes library with a CPython extension named `raytrace.so`.
