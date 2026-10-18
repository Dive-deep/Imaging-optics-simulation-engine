# Project guidance

Read `CODEX_START.md`, `HANDOFF.md`, `README.md`, `FIX_C_REPORT.md`, and `FIX9_REPORT.md` before continuing work. Read the historical `audit_round2/CODE_REVIEW.md` only as the pre-fix evidence. Consult `docs/calculation_spec_ko.md` for equations, references, conventions, and unresolved model questions.

## Current state

- Python sequential optical engine v0.2.1; source to develop is `src/sequential_optics/` at repository root.
- Historical baseline: v0.1.0 had 183 tests; v0.2.0 fixed R1–R9 with 52 regressions (235 total). v0.2.1 fixes C1–C4 and adds explicit per-field FFT MTF axes, with 42 more regressions (277 total). Three independent post-fix role reviews completed. D1 code is unchanged; its component-wise F-theta definition is documented. GitHub migration was authorized on 2026-10-18. See `CODEX_START.md` for environment setup.
- `audit_round2/input/optical_engine/` retains the historical source and review files unchanged. Large generated outputs are omitted from Git as listed in `MIGRATION_MANIFEST.json`; the original ZIP preserves them. Implement future changes in root `src/` and `tests/`.
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
Examples: `python examples/run_demo.py`, `python examples/run_source_wave.py`, and `python examples/run_field_mtf.py`.

For each confirmed defect, add an expected-behavior regression, make the targeted correction, and preserve independent physical controls. After substantive engine changes, perform independent optics-physics, lens-design, and software reviews; user explicitly requires this workflow. Record actual results and remaining model limitations. Do not rewrite historical review outcomes as though they were obtained with corrected code.

Respond to the user in Korean. GitHub migration was authorized on 2026-10-18. See `CODEX_START.md` for environment setup. The R1–R9 and C1–C4 reviewed fixes are completed; continue from this repository; do not repeat the import.
