# Candidate report rendering review

The renderer in `render_report.py` generates one physical candidate page per stored record, an offline HTML viewer, a multipage PDF, and a layout contact sheet. It does not generate missing MTF data or infer a successful wave result.

- Layouts use the exact stored surface sag and fresh 588 nm meridional traces at 0, 0.5, 0.8, and 1F. Chief and both marginal STOP rays are included. Failed traces end at their last recorded point and are marked with a cross. No failed segment is extended to the image.
- Separate clear radii on front/back surfaces are honored. Connecting clear-aperture endpoints is a schematic edge, not a certified mechanical blank shape.
- Spot panels use an explicitly labeled visual 8×24 sample and a common RGB centroid at the same fixed image plane. Summary spot numbers are read from the independent full-sampling radiometry result, not recomputed from this sparse display.
- MTF curves are loaded from calculated wave data. Primary scalar pupil autocorrelation, optional actual FFT samples, and geometric diagnostics have distinct labels. Not-applicable and missing results are distinguished. RGB complex-OTF combination is supplied by the wave evaluator and not redefined by the renderer.
- Relative illumination is labeled as the uniform-radiance, uncoated Fresnel, chief-map approximation. The display is not an OpticStudio irradiance certification. EPD uses the radiometry result and remains null where nominal STOP clipping prevents declaring a complete circular entrance pupil.
- Pilot 0083 was rendered and visually inspected. The surface table, chief/marginal traces, RI/spot/distortion values and absent-wave annotations were checked. Output paths are under `results/report_pilot_v2/` and are development artifacts, not final selected candidates.
- HTML data is embedded with `<` escaped and accessed using text content. Figures use embedded data URLs. No external script/image request is needed. Node syntax validation and an embedded-data round trip passed. Summary NA matched the independent radiometry field exactly.

Run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/candidates_20261022/render_report.py FINAL_SELECTED_JSON --output experiments/candidates_20261022/reports
```

Input accepts a list or `{"candidates": [...]}`. Each item has `candidate_id`, `design`, `metrics` (radiometry output), and `wave` (wave evaluator output). Final record count and IDs are preserved in `render_manifest.json`.

## Final selection display update

The final candidate schema can have unavailable scalar wave curves at 0.8F/1F because the declared angular-pupil model folds. An entire field panel switches to explicitly labeled **geometric MTF diagnostic** when any of its RGB scalar curves is missing; different model types are not silently mixed in a single panel. All three available RGB geometric T/S curves remain visible. The black RGB geometric curve comes from `wave.geometric_polychromatic_diagnostic`, or the identically defined stored `geometric_mtf_proxy` if that field has not yet been assembled.

Read-only verification of `wave_evaluation.geometric_otf` confirmed it evaluates `exp(-2πif x)` at absolute shared sensor coordinates. Its reported centroid is not subtracted independently per wavelength. The root `summary_proxy` subtracts one common RGB centroid, which adds only a common complex phase and preserves polychromatic MTF magnitude. Neither implementation averages monochromatic MTF magnitudes.

C01's final page was rendered from frozen design + completed wave evidence. All four field panels have calculated curves: scalar/autocorrelation diagnostics at 0/0.5F, and labeled geometric RGB diagnostics at 0.8/1F. Underresolved scalar values remain labeled as unconverged. Its schematic clear-aperture edge joins, optical surface prescription table, chief/marginal rays, independently evaluated NA/RI/distortion/spot labels were visually inspected. Candidate IDs and final selection are not altered by the renderer.

## Final artifact and interface QA

Final `candidates.json` generated 24 PDF pages, 24 per-candidate PNG figures, one 24-layout contact sheet, and one standalone HTML viewer. The 96 field panels contain 288 RGB tangential/sagittal curve pairs, with unavailable scalar fields explicitly replaced by geometric diagnostic panels. Summary RI and spot values exactly match the stored full-sampling metrics. The header and HTML metric cards specify 588 nm and RMS radius, avoiding confusion with the common-centroid RGB spot plots. Machine-readable checks are in `results/report_content_qa.json`.

Final C01 and C24 pages and the complete contact sheet were visually inspected. Surface order, different front/back clear radii, chief and marginal traces, displayed values, and conditional/underresolved model labels are visible without page clipping. Desktop and mobile HTML screenshots were also inspected. The full-page scientific figure is intentionally scaled on a narrow phone viewport; the accompanying metric cards and controls reflow independently.

Headless Chromium interface checks passed: initial C01 selection, next/previous navigation, C24 search, selected-candidate JSON download, empty-search disabled navigation, zero strict-EPD50 matches, and 24 first-lens-diameter50 matches. All 24 unique records, summaries and embedded images were present, with zero JavaScript errors and zero external network requests. The exact generated HTML was loaded with Playwright `page.setContent` because the managed browser administrator policy blocks `file://` navigation; no browser policy was changed. This verifies the standalone document's rendering and JavaScript, while direct local-file opening remains dependent on the viewer's browser policy. Evidence is in `results/report_browser_qa.json` and the two browser screenshots.

Final principal artifacts are `reports/candidate_report.html` (25.3 MB), `reports/candidate_layouts_and_mtf.pdf` (1.50 MB), `reports/candidate_contact_sheet.png`, and `reports/report_summary.json`. All figures use calculated data; the report does not assert that high-NA scalar wave results have achieved 0.5% convergence or Zemax equivalence.
