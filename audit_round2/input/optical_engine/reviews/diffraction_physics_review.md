# Independent diffraction physics code review

Reviewer role: real-ray/geometry implementer independently reviewing the other author's `wave.py`.

Scope: actual ray eikonal → physical boundary → referenced OPD → scalar FFT/RS → complex OTF/MTF and angular power. This review evaluates the declared independent profile in `docs/calculation_spec_ko.md`, especially W.1–W.14. It does not establish Zemax numerical equivalence or a Maxwell solution.

## Findings and disposition

1. **Focal snapshot reused for afocal analysis silently changed illumination.** The focal snapshot's original uniform Fourier-pupil convention implies a nonuniform physical boundary amplitude A=γ₀²/L₀. The afocal path initially reset A=1 when reusing that snapshot. This violated frozen-boundary illumination. The author added explicit rejection of cross-mode reuse for that case. The existing fixed-boundary test now includes this rejection and passes.
2. **Insufficient FFT padding must affect the high-level validity state.** A finite support sampled without enough zero padding produces circular rather than linear pupil autocorrelation. Low-level `fft_psf` already reports this. The high-level status originally depended only on phase sampling; the author changed it so inadequate linear-autocorrelation padding also produces `underresolved`, with a regression test. Even the remaining `conditional` status is not an unconditional accuracy certificate.
3. **T/S axis names require explicit provenance, not a formula change.** Initial review suspected that global projected X/Y slices should be rotated to meridional/radial axes. Re-reading W.5 establishes that the document deliberately defines this independent FFT profile as S=f_x, T=f_y in its stored chief frame. Therefore this is not a defect against the specification. The author added explicit axis-convention metadata. These labels are not interchangeable with the geometric metrics module's field-aligned meridional/sagittal axes.

## Physics checks performed

- Infinite-field launch OPL uses a shared plane-wave eikonal; finite object rays share their source point. Transport to the chief-normal plane adds n times the actual positive distance.
- The boundary plane bounds the entire declared last clear surface, rather than relying on a maximum of sparse ray samples.
- Focal OPD uses S−S_chief+n(R−L). The reference sphere's position chart χ=u/R, η=v/R is distinct from the aberrated ray's actual direction cosines.
- The physical-area Jacobian L²/γ⁴ and Fourier amplitude B=A L/γ² are consistent with the declared local scalar RS expansion. Through-focus reuse retains A and changes B accordingly.
- FFT output coordinates scale with λ_medium=λ_vacuum/n; the OTF retains complex sign and phase. Same-pupil diffraction comparison removes phase while retaining amplitude and support.
- Far-field spectral power uses γ|F|² dl dm, while density per solid angle uses γ²|F|². Exact angular chart conversion is separate from treating direction cosine as angle.
- Physical aperture failures and unresolved ray failures remain distinct. Original STOP topology is retained in pupil interpolation; successful-only triangulation does not bridge known clipped holes.

## New independent tests

`tests/test_review_wave.py` adds four tests whose oracles are separate physical identities:

1. **Uniform circular aperture, exact on-axis complex RS field.** Direct integration must equal exp(ikL)−(L/R_edge)exp(ikR_edge). This checks complex phase, absolute amplitude and the obliquity kernel together, not just a normalized PSF.
2. **Tilted periodic plane wave, raw power.** A single exact spectral mode must carry n|E|² A cosθ/(2Z₀). This catches omitted or duplicated angular-cosine factors and incorrect mm²→m² conversion before normalization can hide them.
3. **Immersed focal frequency scaling.** Changing output index from 1 to 1.5 scales the physical PSF coordinates and matching OTF frequency axis by the corresponding medium-wavelength ratio.
4. **Actual refracting lens eikonal-gradient identity.** Independently re-aimed STOP perturbations, in both finite-object and oblique infinite-object modes, verify ∇⊥S=n d⊥ on the same physical boundary. This connects ray OPL, material index and frame projection, rather than constructing a perfect pupil in isolation.

At review execution, all 27 tests in `test_wave.py` and `test_review_wave.py` passed. The four new independent tests were included in that count. Later suite totals may change as other reviewers add cases.

## Remaining model limits

- Piecewise-linear unwrapped OPD and conservative sampled pupil polygons require independent ray-grid and FFT-grid refinement. A small phase increment is necessary, not sufficient, for interpolation accuracy.
- The local near-focus Fourier expansion must be assessed against direct scalar propagation over the region of interest. The provided full-window residual is a diagnostic, not a rigorous MTF error bound.
- Sampled Jacobian orientation and overlap checks do not prove global ray-map uniqueness for arbitrary high-order aspheres.
- Scalar, no-coating profile calculations do not become vector diffraction or vendor-compatible just because the analytic tests pass.
- The far-field result is normalized over the represented propagating angular window; omitted-window power still needs convergence checks.

The review supports implementing and exercising the documented independent scalar profile with explicit diagnostics. It does not support labeling the engine as validated against commercial software or measured lenses.

## Additional independent review: finite-source physical diffraction bridge

Reviewed `source_wave.py` after the physical source-conditioned path was added.

The source measure is calculated from the initial unit-direction spherical triangles, using exact great-circle solid angle and projected solid angle. The latter is the Lambertian integral of cosθ_source dΩ, not ordinary solid angle with a missing source cosine. Stable `atan2` edge angles retain small-angle precision. These geodesic triangles approximate a refracted nonlinear STOP-to-source map; map and transmission refinement are still required.

Each triangle transports already attenuated source power to its actual physical boundary area. The electric amplitude uses the actual mean crossing cosine, absolute output index and mm²→m² conversion. There is no extra n² multiplication of a conserved power packet. The trace transmission is included once; later normalized kernel synthesis does not apply Fresnel, absorption or fill again.

The raw ray/triangle, raster-boundary and angular-spectrum powers are retained before PSF normalization. The triangle closure is explicitly acknowledged as construction, while raster and spectral ratios are independent numerical diagnostics. Focal-reference OPD is not reused as the physical boundary phase for that raw spectral check. An inadequately sampled raw phase remains flagged rather than being promoted to a successful power check.

Two changes were needed during integration:

- **Physical STOP support was absent from Gauss quadrature nodes.** The outer quadrature ring lies inside the rim. Treating its hull as the wave aperture caused approximately 5% missing collection in a five-ring circular-aperture fixture, and adding azimuth samples could not fix that radial omission. `trace_bundle` now optionally includes actual rim samples and the chief with zero geometric quadrature weight, enabled by `SamplingSettings.include_pupil_boundary`. This fixes support reconstruction while leaving weighted geometric moments unchanged. Azimuthal polygon and clipped-edge refinement remain necessary.
- **Rim aiming and clipping tolerances initially disagreed.** Subsequent camera review found five falsely clipped rim samples for a finite-object triplet field: the default aiming tolerance was looser than the actual STOP aperture predicate. The aiming termination condition now requires both its xy residual tolerance and the identical physical STOP radial predicate used in the final trace. It performs another Newton step when needed; neither the aperture nor the requested pupil coordinate is enlarged or shifted. Both the geometry regression and the independently authored triplet rim regression pass.
- **Source focal propagation needed the same approximation diagnostics as the lens metric.** The author added full-window local-expansion phase/amplitude residuals and optional exact-RS comparison points to the source path. These are model diagnostics, not a bound on MTF error.

The detector bridge conserves each Fourier cell's power by polygon overlap on a common global grid. It neither recenters each source image nor renormalizes a cropped tail. At chief tilt, the bridge explicitly requires opt-in to a local chief-parallel affine projection; it does not claim exact wave propagation to an inclined detector. Afocal combination rotates full unit vectors to the global frame before binning.

`tests/test_review_source_wave.py` adds four independently designed tests; all four pass:

1. Changing a homogeneous medium's absolute index from 1 to 1.44 leaves collected Lambertian power unchanged and scales the physical electric amplitude by 1/√1.44. Both raw boundary and raw wave powers are compared against the analytic circular-aperture collection, before normalization, with a declared coarse-grid tolerance.
2. A nearly paraxial, index-matched absorbing plate with μ=ln(2)/mm halves collected and synthesized detector power. This detects applying transmission twice in the Q–kernel interface.
3. An independently constructed tilted-frame single-cell response lands at the global detector location obtained by an explicit 3-D line-plane intersection, while preserving total cell power.
4. A same-field bundle traced from an old object-z is rejected when reused with a different physical source distance; matching only field metadata and wavelength is insufficient.

The author's seven source-wave integration tests also pass after the support fix. They cover the projected spherical measure against an independent surface integral, analytic Lambertian aperture collection, radiance/fill/source-area scaling, cropped image power, angular registration, explicit off-axis projection, and invalid source/absorbing-output inputs.

Remaining limits include piecewise-constant triangle amplitude and vertex-mean transmission, sampled rather than globally certified map regularity, source macro-cell center quadrature requiring area refinement, scalar polarization reduction, and a nonabsorbing output medium for wave propagation. A normalized Q×h image is not by itself evidence that the raw wave power or physical approximation is accurate.
