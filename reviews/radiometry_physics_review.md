# Independent physics review: radiometry and source image synthesis

Reviewer role: optical physicist, distinct from module author. Date: 2026-10-18.
Reviewed: `radiometry.py`, `imaging.py`, shared geometry contracts and the source-to-wave interface.
Evidence: `tests/test_review_radiometry.py`, 5 independent regression tests pass. These are mathematical identities and deliberately simple synthetic fixtures, not comparison with Zemax or measurements.

## Findings and dispositions

1. **Reverse RI media, Snell directions, and radiance measure are consistent.** Reversing surface order and reflecting z requires negating curvature and every even coefficient, retaining conic, and changing each `medium_after` to its original incident medium. The implementation does this. On a flat image detector, `cos(theta) dOmega = dl dm`; therefore no extra detector cosine belongs in the direction-cosine quadrature. The factor `(n_image/n_source)^2` belongs in radiance transport, not packet power. The new single planar interface test uses n=1 to n=1.5: the reverse TIR escape cone and index-squared factor cancel to give E=pi L with disabled losses. This simultaneously challenges index orientation, reverse Snell, and the integration measure.

2. **Unpolarized Jones-power transport retains polarization history.** Each transverse input basis vector is propagated coherently through the changing S/P bases; their two output powers are averaged incoherently. A nonmeridional curved two-surface test retraces the actual ray backwards and verifies equal direct transmission. The scalar-interface-average path deliberately resets polarization and is correctly flagged as an approximation. Neither path is vector diffraction: no coating phase or coherent interference among pupil-dependent polarization vectors is supplied to the scalar FFT.

3. **Equal-index grazing Fresnel branch corrected.** The initial function tested `sin(theta_t)^2 >= 1` before the equal-index branch, incorrectly labelling n1=n2, cos(theta)=0 as TIR. The author moved the equal-index interface case first. Independent regression now passes. At grazing incidence this is mainly a physically correct boundary condition and API edge case; emitted projected power at exactly grazing angle remains zero.

4. **Absorbing output medium in a literal far-field calculation corrected.** A ray's last-surface power cannot be called power at infinity if a positive homogeneous absorption coefficient continues downstream. The author added an explicit rejection for supplied nonzero final-medium absorption, unless the caller intentionally selects `lossless_scalar`. Missing absorption remains a declared transparent assumption, not fabricated data. Finite image kernels separately propagate bulk absorption from the final surface to the detector.

5. **Power, radiance, source area and fill are separated correctly.** Source hemisphere power is pi times panel area times panel-average Lambertian radiance. Cosine hemisphere sampling assigns constant packet weight; entrance-plane importance sampling uses both projected cosines and r^-2. The full finite-source path uses a joint 4D sample and its corresponding measure. `panel_average` avoids applying fill a second time. Image synthesis multiplies a matched Q by source-area weight and unit/cropped kernel mass; it does not reapply Fresnel, a geometric warp Jacobian, or fill.

6. **Reported-bin solid angle and cropping are physically explicit.** Spherical bins use `(cos(theta0)-cos(theta1))*Delta_phi`; the slope rectangle uses an exact corner formula. Finite image kernels retain mass below one when the reported ROI crops the response. Independent unequal-area pixel tests verify raw detector power, density conversion, and cropped power separately. Normalizing the cropped image would hide loss and is not done.

7. **Kernel provenance and wavelength registration are explicit contracts, not proof.** `SourceKernel` rejects a declared lens-metric pupil and requires global physical coordinates. The registered spectral sum preserves different image positions and adds intensities. A user can still mislabel arbitrary supplied arrays; type metadata cannot certify their physics. A production source-conditioned diffraction builder should construct Q and amplitude from the same source/ray data rather than requiring the user to assert provenance.

## Remaining model and numerical limits

- RI midpoint direction grids can miss a small transmitted cone completely; zero accepted rays is not proof of physical zero illumination. The result reports sample counts and assumptions, and refinement/coverage checks are still required. A truncated cone must be justified externally.
- Geometric image and far-field paths retain aberrated ray directions and positions but do not include diffraction. They should not be described as the requested FFT image kernel by themselves.
- Finite-source macro-cell midpoint approximation requires both source-grid and angular-ray refinement; adding only rays per cell cannot resolve missing source/bitmap structure.
- The direct-path engine excludes scattered light, lens ghosts, coherent etalons, surface roughness and polarization-vector diffraction. Its transparent missing-absorption option is an explicit modelling choice.
- A Jones power ledger is valid only for successful direct paths. Aggregate nondirect power includes clipping, misses, Fresnel reflection and absorption; it is not an absorption measurement.

## Source-conditioned wave bridge recommendation

A feasible bounded first implementation can reuse the boundary state's original STOP triangulation:

1. Require a finite physical source point and its actual first-segment directions. The STOP coordinates parameterize the map, not a solid-angle measure.
2. Compute each launch direction triangle's solid angle and the Lambertian projected solid-angle integral. Use the actual transmission along that same triangle's rays; retain its numerical quadrature approximation.
3. Compute its transmitted power per source area `Q_triangle = L_panel integral(tau cos(theta_source) dOmega)` and multiply by an explicit source cell area when constructing electric fields.
4. Map the same triangle to the canonical physical boundary and divide by its physical area to obtain boundary-plane power density. Do not divide by a direction-cosine/Fourier-chart area instead.
5. Convert with `A_E=sqrt(2 Z0 p_boundary /(n cos(theta_actual)))`, after converting W/mm2 to W/m2. This actual ray cosine is distinct from the reference-sphere gamma.
6. Reconstruct the amplitude field, retain raw boundary/ray power, and pass the same phase and fixed boundary into the wave engine. Use its source amplitude path, not the uniform Fourier pupil.
7. Compare raw angular-spectrum wave power with boundary and ray power before normalizing a PSF. Also refine STOP triangles and boundary amplitude reconstruction. A normalized `Q*h` integral alone is an identity and cannot validate this bridge.
8. Mark the scalar and sampled reconstruction limitations. A fold, multiple phase branches or numerical aiming failure is an unsupported state, never a dark pupil.

## Follow-up implementation and independent review handoff

The proposed bridge was subsequently implemented as `source_wave.py` with seven author tests in `test_source_wave.py`. It uses exact spherical projected vector area, vertex-mean transmission, and constant physical electric-field amplitude per mapped triangle. It returns matched Q, a focal/afocal wave response, unnormalized ray/boundary/wave power diagnostics, optional conservatively remapped `SourceKernel`, and global-frame angular binning. Off-axis detector projection requires explicit opt-in and is labelled as a local parallel-ray approximation. The focal path also exposes local-expansion residuals and optional exact scalar RS checks.

The analytic aperture-power test found that the original Gauss STOP quadrature contained no actual pupil rim. Increasing azimuth alone therefore retained about 5% missing source power at five radial rings. The geometry author added configured zero-weight rim and chief rays; quadrature weights are unchanged and the physical support is now represented. The aperture projected-solid-angle result converges to its independent analytic value with boundary refinement.

Because this reviewer also authored the follow-up bridge, its independent review was assigned to the separate geometry/diffraction reviewer (`test_review_source_wave.py` and the diffraction review), with an additional camera integration review. The seven author tests do not substitute for that independent review. Coarse triangular reconstruction, scalar boundary physics, clipped-edge refinement and local detector projection remain explicit limits rather than guaranteed accuracy for arbitrary lenses.
