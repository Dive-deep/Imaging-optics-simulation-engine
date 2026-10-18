# Project guidance

Read `HANDOFF.md`, `README.md`, and `FIX9_REPORT.md` before continuing work. Read the historical `audit_round2/CODE_REVIEW.md` only as the pre-fix evidence. Consult `docs/calculation_spec_ko.md` for equations, references, conventions, and unresolved model questions.

## Current state

- Python sequential optical engine v0.2.0; source to develop is `src/sequential_optics/` at repository root.
- Historical baseline: 183 tests passed. A subsequent independent review reproduced 9 additional defects; v0.2.0 fixes all nine with 52 new regressions. The complete 235-test suite passes. Migration is postponed.
- `audit_round2/input/optical_engine/` is the unchanged review snapshot. Preserve historical evidence and implement future changes in root `src/` and `tests/`.
- Reproduction assertions under `audit_round2/` demonstrate existing bugs; their passing is not a correctness acceptance test. Some scripts deliberately import the frozen snapshot. To verify a fix, port the expected behavior into root tests or explicitly target the working source.

## Optical and software contracts

- Real-ray intersections, vector Snell and actual OPL underpin analyses. Paraxial definitions remain appropriate for EFL/BFL and initial aiming.
- Keep mm, vacuum nm, absolute refractive index, coordinate frames, pupil measure, and optical power versus field amplitude explicit.
- G is the configurable default wavelength (packaged example 550 nm). Nominal monochromatic work must not require RGB, CTE, dn/dT, or absorption data it does not use.
- Preserve modular APIs and explicit settings. Do not silently infer physical catalog properties from a material name.
- Lens-only FFT MTF and physical Lambertian source illumination are separate models.
- Uniform temperature only; stress, GRIN, birefringence, mount FEA and manufacturing tolerances are outside the accepted scope.
- Physical clipping, TIR, missing data, numerical failure, and insufficient sampling must remain distinguishable. Do not turn unresolved rays into zero power to make an analysis succeed.
- Do not claim Zemax numerical equivalence or measured accuracy from analytic invariants or passing tests.

## Validation when fixes are requested

Install: `python -m pip install -e .`.
Run the root suite: `python run_tests.py`.
Examples: `python examples/run_demo.py` and `python examples/run_source_wave.py`.

For each confirmed defect, add an expected-behavior regression, make the targeted correction, and preserve independent physical controls. After substantive engine changes, perform independent optics-physics, lens-design, and software reviews; user explicitly requires this workflow. Record actual results and remaining model limitations. Do not rewrite historical review outcomes as though they were obtained with corrected code.

Respond to the user in Korean. Migration is postponed. The 9 reviewed fixes are completed; do not initiate migration without a new request.
