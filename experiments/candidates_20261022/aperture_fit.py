"""Fit unused clear apertures to a dense, freshly traced footprint envelope.

This changes the explicit design aperture, not the trace mask after calculation.
By default each element keeps equal front/back radii. Optional independent-face
mode specifies only clear apertures, not a proven manufacturable edge profile.
No existing transmitted sampled ray may be removed. Both the before and after
traces remain subject to comparison.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from model import trace_design, native_surfaces, sag, WAVELENGTHS, geometry_checks


def _digest(design):
    return hashlib.sha256(json.dumps(design, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _compare(before, after, mask, key, atol):
    a = np.asarray(before[key])[mask]
    b = np.asarray(after[key])[mask]
    same_finite = np.array_equal(np.isfinite(a), np.isfinite(b))
    finite = np.isfinite(a) & np.isfinite(b)
    delta = np.abs(a[finite] - b[finite])
    return {
        "array": key,
        "compared_scalar_count": int(a.size),
        "bitwise_equal": bool(a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()),
        "nonfinite_pattern_equal": bool(same_finite),
        "max_absolute_difference": float(delta.max()) if len(delta) else 0.0,
        "absolute_tolerance": atol,
        "passed": bool(same_finite and (not len(delta) or delta.max() <= atol)),
    }


def fit_apertures(design, fieldfractions=None, rings=16, arms=48, wavelengths=WAVELENGTHS, *, coupled=True):
    """Return (new_design, audit); reject a changed prescription on any lost ray.

    The 0.30 mm /2% margin is a sampled engineering envelope, not a continuum
    extremum proof, manufacturing margin, or permission to omit vignetting.
    """
    fields = np.linspace(0.0, 1.0, 13) if fieldfractions is None else np.asarray(fieldfractions, float)
    wavelengths = np.asarray(wavelengths, float)
    original = copy.deepcopy(design)
    before = trace_design(original, fields=fields, wavelengths=wavelengths, rings=rings, arms=arms)
    success = before["status"] == 0
    audit = {
        "schema": "dense-footprint-aperture-fit/1",
        "original_design_sha256": _digest(original),
        "field_fractions": fields.tolist(),
        "field_angles_deg": (12.0 * fields).tolist(),
        "wavelengths_nm": wavelengths.tolist(),
        "rings": rings,
        "arms": arms,
        "chief_and_rim_included": True,
        "sampled_ray_count": int(success.size),
        "original_transmitted_ray_count": int(np.count_nonzero(success)),
        "original_status_counts": {str(int(k)): int(v) for k, v in zip(*np.unique(before["status"], return_counts=True))},
        "coupled_element_radii": bool(coupled),
        "policy": "Coupled mode: equal front/back radius, capped by smaller original face. Independent mode: each face capped by its own original radius. Reserve max(0.30mm,2% of maximum sampled footprint). Front-diameter50 first front face remains exactly25mm radius; coupled mode keeps both first-element faces25mm.",
        "limitation": "Sampled nonnegative Y fields in a rotationally symmetric system. No rigorous extremum proof between fields/rings/azimuths and no manufacturing tolerance allocation. Dense downstream qualification remains necessary.",
        "mechanical_profile": "Only clear semi-diameters are specified. Mechanical outer diameter, edge/chamfer or tapered solid profile is unspecified; independent clear faces are not manufacturing certification. The wave model treats its declared last clear surface extent as the boundary constraint; moving that boundary with identical rays is model sensitivity, not a physical improvement.",
    }
    if not success.any():
        audit.update(accepted=False, applied=False, reason="No transmitted rays; no aperture inference", result_design_sha256=_digest(original))
        return original, audit
    radii = np.sqrt(np.sum(before["history"][success, :, :2] ** 2, axis=-1))
    surface_maximum = np.max(radii, axis=0)
    new = copy.deepcopy(original)
    old_radii = np.asarray(original["clear_radii_mm"], float)
    fitted = old_radii.copy()
    elements = []
    groups = [(2*i, 2*i+2) for i in range(3)] if coupled else [(i, i+1) for i in range(6)]
    for start, end in groups:
        foot = float(np.max(surface_maximum[start:end]))
        margin = max(0.30, 0.02 * foot)
        cap = float(np.min(old_radii[start:end]))
        requested = max(0.5, foot + margin)
        if start == 0 and original["aperture_mode"] == "front_diameter50":
            proposal = 25.0
            if cap < 25.0:
                audit.update(accepted=False, applied=False, reason="Original first aperture is not50mm; fitting cannot enlarge it", result_design_sha256=_digest(original))
                return original, audit
        else:
            proposal = min(cap, requested)
        fitted[start:end] = proposal
        elements.append({
            "lens": start // 2 + 1,
            "surfaces": list(range(start+1, end+1)),
            "original_face_radii_mm": old_radii[start:end].tolist(),
            "maximum_sampled_footprint_radius_mm": foot,
            "requested_margin_mm": margin,
            "proposed_equal_face_radius_mm": proposal,
            "actual_sampled_clearance_mm": proposal - foot,
            "full_requested_margin_available": bool(proposal >= requested),
        })
    new["clear_radii_mm"] = fitted.tolist()
    audit["surface_footprint_maxima_mm"] = surface_maximum.tolist()
    audit["elements"] = elements
    audit["proposed_design_sha256"] = _digest(new)
    audit["original_geometry"] = geometry_checks(original, points=513)
    audit["proposed_geometry"] = geometry_checks(new, points=513)
    # Shrinking a spherical clear radius cannot violate a previously valid
    # spherical domain, but explicitly check all six domains rather than assume.
    domains = audit["proposed_geometry"]["minimum_conic_domains"]
    if min(domains) < 0 or min(fitted) < 0.5:
        audit.update(accepted=False, applied=False, reason="Proposed surface domain/radius invalid", result_design_sha256=_digest(original))
        return original, audit
    after = trace_design(new, fields=fields, wavelengths=wavelengths, rings=rings, arms=arms)
    after_success = after["status"] == 0
    lost = success & ~after_success
    gained = ~success & after_success
    changed_status = before["status"] != after["status"]
    comparison = [
        _compare(before, after, success, "history", 1e-10),
        _compare(before, after, success, "direction_history", 1e-12),
        _compare(before, after, success, "opl_history", 1e-10),
        _compare(before, after, success, "pos", 1e-10),
        _compare(before, after, success, "dir", 1e-12),
        _compare(before, after, success, "image_xy", 1e-10),
    ]
    flags_equal = bool(np.array_equal(before["image_forward"], after["image_forward"]))
    passed = not np.any(changed_status) and flags_equal and all(x["passed"] for x in comparison)
    audit.update({
        "lost_transmitted_ray_count": int(np.count_nonzero(lost)),
        "newly_transmitted_ray_count": int(np.count_nonzero(gained)),
        "changed_ray_status_count": int(np.count_nonzero(changed_status)),
        "all_ray_statuses_equal": bool(np.array_equal(before["status"], after["status"])),
        "all_image_forward_flags_equal": flags_equal,
        "transmitted_ray_comparison": comparison,
        "transmitted_optics_bitwise_equal": bool(all(x["bitwise_equal"] for x in comparison)),
        "accepted": bool(passed),
        "applied": bool(passed and not np.array_equal(old_radii, fitted)),
        "reason": "All sampled statuses and transmitted optical histories preserved" if passed else "Rejected: a status or transmitted optical history changed",
    })
    result = new if passed else original
    audit["result_design_sha256"] = _digest(result)
    return result, audit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--independent-faces", action="store_true")
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    design = source.get("design", source)
    result, audit = fit_apertures(design, coupled=not args.independent_faces)
    output = {**source, "design": result, "aperture_fit": audit} if "design" in source else {"design": result, "aperture_fit": audit}
    from optimize import quick_metrics
    output["quick_metrics"] = quick_metrics(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: audit.get(k) for k in ("accepted", "applied", "sampled_ray_count", "lost_transmitted_ray_count", "changed_ray_status_count", "transmitted_optics_bitwise_equal", "reason")}, indent=2))
    print("clear radii:", result["clear_radii_mm"])
