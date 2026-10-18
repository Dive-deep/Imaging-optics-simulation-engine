# Independent software and reproducibility review

Reviewer role: implementation author of material/thermal modules, reviewing **other authors'** `models.py`, `settings.py`, `io.py`, `api.py`, `__main__.py`, and geometry/paraxial integration. Material/thermal formula correctness is not self-certified by this review.

Review date: 2026-10-18. Scope: declared units, selective data requirements, run snapshots, user-input failure semantics, portable JSON, reusable APIs, and deterministic conventions. No Zemax output or physical measurement is available.

## Findings and concrete reproductions

| ID | Severity | Finding | Reproduction / effect | Status |
|---|---|---|---|---|
| SW01 | Medium | Count settings accept noninteger values | `EngineSettings(ray=RaySettings(aim_max_iterations=1.5))` succeeded, but later `range(1.5)` raised `TypeError`; fractional frequency counts similarly reached NumPy invalid count calls. | Fixed: integral count validation; regression passes |
| SW02 | Medium | STOP index did not require an integer | `replace(system, stop_index=.5)` succeeded when no radius override existed, then STOP indexing raised `TypeError`. | Fixed: integral STOP index validation; regression passes |
| SW03 | High | Ray direction norm overflow silently produced zero direction | `Ray([0,0,0], [1e308,0,1e308])` passed finite component checks, overflow in norm yielded infinity, division yielded zero direction. This violated the central unit-direction contract rather than reporting invalid data. | Fixed: scale by largest component before normalization; regression passes |
| SW04 | Medium | Settings and material registries aliased caller dictionaries | Changing the originally supplied `{'G':550}` mapping to 532 changed the frozen settings object's behavior. Replacing caller-owned material map entries changed a frozen system's material without constructing a new system. | Fixed: immutable copied mappings; regression passes |
| SW05 | Medium | Portable thermal lens-group JSON was not adapted | `thermal_options={'lens_groups':[{'front_index':0,'back_index':1}]}` reached `g.back_index` on a dict and raised `AttributeError` outside the structured error path. | Fixed: RunRequest converts dict groups to LensGroup; regression passes |
| SW06 | Medium / convention | Report EFL reference was chosen by incidental material registry key | Adding an unused key `'air'` changed `efl_report_mm`, despite unchanged physical prescription. A synthetic lens changed 50.8474576271 to 50.8627118644 mm when unused air n=1.0003 was inserted. Physical image-space EFL stayed unchanged. Formula was not the issue: the reference role needed to be an explicit setting rather than a name heuristic. | Fixed: `EngineSettings.report_reference_material`; default is ambient; both default-invariance and explicit-reference regressions pass |
| SW07 | Medium / API robustness | Some invalid RunRequest structures failed with raw Python exceptions | Fields containing arbitrary objects, unknown thermal kwargs, and nested JSON structures needed boundary validation. Do not fix by catching every exception as a physics failure; validate the actual request schema. | Fixed for reviewed cases: Field validation, allowed thermal keys, nested group conversion; regressions pass |

## Passing behavior before corrections

The new independent suite `tests/test_review_software.py` already verified these positive properties on the initial code:

- A glass with only a measured G index computes first-order results without CTE, dn/dT, absorption, source, or RGB data.
- A missing B index affects B's requested result, retaining the completed G result.
- Missing thermal data at a hot state retains the baseline result.
- Changing the configurable G wavelength from 550 to 532 nm propagates to both lookup and reported wavelength without changing module constants.
- Prescription JSON round-trip preserves a single-wavelength material and the deterministic input hash.
- Complex OTF arrays retain real/imaginary components; nonfinite result values are explicitly tagged rather than emitted as invalid JSON NaN.

Initial run: 12 tests, 8 failed assertions (including four count-setting subcases), 1 error. Each failure is a concrete boundary defect above, not a physical agreement claim. The final rerun status is recorded below after fixes.

## Required interpretation of successful execution

`evaluate` returning a result container does not establish usable physical output for every requested metric. Consumers must inspect nested `status` and wave applicability/sampling diagnostics. A successful input hash is an input identity, not a code revision certificate; retain the distributed source/package version alongside generated data. A conventional report EFL should not be used in place of the explicitly named physical image-space EFL.

## Remaining scope limits

- No vendor numerical equivalence or measurement correlation is established.
- This review does not certify arbitrary pathological geometry, Fourier pupil reconstruction, or radiometric amplitude mapping; other domain reviews cover those calculations.
- Performance results are meaningful only with the stored requested field, wavelength, thermal placement, STOP, focus, sampling and model-applicability metadata.

## Final correction verification

Root corrected the reported boundaries. Independent rerun:

```text
PYTHONPATH=src python -m unittest discover -s tests -p test_review_software.py -v
Ran 14 tests in 0.029s
OK
```

The source changes were read after the rerun: overflow-safe ray normalization, read-only ray arrays, copied immutable configuration mappings, explicit report reference, count/STOP validation, and typed thermal-group conversion are present. No remaining blocking defect was found within this reviewed scope. This is an engineering review and regression result, not a certification of the full physical solver.
