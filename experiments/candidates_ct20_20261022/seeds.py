"""Deterministic, unqualified analytical starts for the 2026-10-22 search.

Lengths are mm; curvatures are inverse mm. Asphere numbers are sag amplitudes
at r=25 mm, not conventional dimensional A4...A12 coefficients. This file
does not claim that a paraxially prepared start passes any real-ray requirement.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
N = {"PMMA": 1.492, "PC": 1.586, "FDS18W": 1.94594, "TAFD45": 1.95375}
V = {"PMMA": 57.47, "PC": 32.0, "FDS18W": 17.98, "TAFD45": 32.32}
TARGET_EFL = 6.4 / np.tan(np.deg2rad(12.0))


def sag(c, k, r):
    disc = 1.0 - (1.0 + k) * c * c * np.asarray(r) ** 2
    if np.min(disc) < 0:
        raise ValueError("conic domain")
    return c * np.asarray(r) ** 2 / (1.0 + np.sqrt(disc))


def first_order(materials, curvatures, distances):
    matrix = np.eye(2)
    n0 = 1.0
    for j, c in enumerate(curvatures):
        n1 = N[materials[j // 2]] if j % 2 == 0 else 1.0
        matrix = np.array([[1.0, 0.0], [-(n1 - n0) * c, 1.0]]) @ matrix
        if j != 5:
            matrix = np.array([[1.0, distances[j] / n1], [0.0, 1.0]]) @ matrix
        n0 = n1
    power = -float(matrix[1, 0])
    efl = 1.0 / power if abs(power) > 1e-14 else 1e14
    bfl = float(matrix[0, 0]) / power if abs(power) > 1e-14 else 1e14
    return efl, bfl, matrix


def prepare_geometry(materials, c, k, mode, offset, existing_distances=None):
    """Conservative sampled rim-clearance construction; no optical qualification."""
    element_radii = [25.0 if mode == "front_diameter50" else 26.0, 24.0, 22.0]
    for i, material in enumerate(materials):
        if material in ("FDS18W", "TAFD45"):
            if i == 1:
                element_radii[i] = 23.0
            elif i == 2:
                element_radii[i] = 19.5
    stop_d = float(offset)
    for _ in range(8):
        if mode == "epd50":
            element_radii[0] = 25.0 + stop_d * np.tan(np.deg2rad(12.0)) + 0.25
        r = np.linspace(0.0, element_radii[0], 257)
        stop_d = max(float(offset), -float(np.min(sag(c[0], k[0], r))) + 0.3)
    radii = np.repeat(element_radii, 2)
    distances = []
    min_clearance = []
    for j in range(5):
        r = np.linspace(0.0, min(radii[j], radii[j + 1]), 257)
        a = sag(c[j], k[j], r)
        b = sag(c[j + 1], k[j + 1], r)
        minimum = 2.4 if j % 2 == 0 else 1.2
        distance = max(minimum, float(np.max(a - b)) + minimum)
        if existing_distances is not None:
            distance = max(distance, float(existing_distances[j]))
        distances.append(distance)
        min_clearance.append(float(np.min(distance + b - a)))
    efl, bfl, matrix = first_order(materials, c, distances)
    rear_sag = sag(c[-1], k[-1], np.linspace(0.0, radii[-1], 257))
    image_gap = max(3.2, bfl, float(np.max(rear_sag)) + 0.5)
    distances.append(image_gap)
    return radii, distances, stop_d, efl, bfl, matrix, min_clearance


def analytic_start(materials, mode, family, variant):
    glass = next(i for i, material in enumerate(materials) if material not in ("PMMA", "PC"))
    plastic = [i for i in range(3) if i != glass]
    # Thin powers are only an initial direction. Thick-lens compensation below
    # changes them and records its achieved EFL, BFL and remaining discrepancy.
    powers = np.zeros(3)
    if family == "negative_glass_positive_plastics":
        vp = sum(V[materials[i]] for i in plastic) / 2.0
        vg = V[materials[glass]]
        color_ratio = min(0.70, max(0.25, vg / vp))
        total_positive = (1.0 / TARGET_EFL) / (1.0 - color_ratio)
        powers[glass] = -total_positive * color_ratio
        split = (0.43, 0.57)[variant % 2]
        powers[plastic[0]] = total_positive * split
        powers[plastic[1]] = total_positive * (1.0 - split)
    else:
        powers[glass] = 0.035 + 0.005 * (variant % 2)
        powers[plastic[0]] = 0.022
        powers[plastic[1]] = -0.022
    bend = np.array([0.18, -0.12, -0.15]) * (1 if variant % 2 == 0 else -1)
    k = np.repeat([-1.0 if material in ("PMMA", "PC") else 0.0 for material in materials], 2)
    # A modestly hyperbolic plastic can form another bending family without a
    # conic-domain singularity. Glass always remains exactly spherical.
    if variant % 2:
        k[np.repeat([m in ("PMMA", "PC") for m in materials], 2)] = -2.0

    def curvatures_from_power(p):
        c = []
        for i, material in enumerate(materials):
            delta = p[i] / (N[material] - 1.0)
            c.extend([delta * (0.5 + bend[i]), delta * (-0.5 + bend[i])])
        return np.asarray(c)

    # Domain-safe power bounds also keep unusually strong plastic starts finite.
    upper = np.array([0.11, 0.11, 0.11])
    upper[glass] = (0.9 / (28.0 if glass == 0 else 24.0 if glass == 1 else 20.0)) * (N[materials[glass]] - 1.0) / max(abs(0.5 + bend[glass]), abs(-0.5 + bend[glass]))
    lower = -upper
    powers = np.clip(powers, lower * 0.95, upper * 0.95)
    prior = powers.copy()

    def residual(p):
        try:
            c = curvatures_from_power(p)
            radii, distances, stop_d, efl, bfl, matrix, clearance = prepare_geometry(materials, c, k, mode, 0.75 + 1.5 * variant)
            power = -matrix[1, 0]
            return np.array([
                2500.0 * (power - 1.0 / TARGET_EFL),
                0.04 * min(0.0, bfl - 4.0),
                0.005 * max(0.0, sum(distances) - 130.0),
                *((p - prior) * 0.25),
            ])
        except (ValueError, FloatingPointError):
            return np.array([1e4, 1e4, 1e4, *(p - prior)])

    fit = least_squares(residual, powers, bounds=(lower, upper), max_nfev=80, ftol=1e-9, xtol=1e-9, gtol=1e-9)
    c = curvatures_from_power(fit.x)
    radii, distances, stop_d, efl, bfl, matrix, clearance = prepare_geometry(materials, c, k, mode, 0.75 + 1.5 * variant)
    return {
        "seed_id": "",
        "materials": list(materials),
        "curvatures": c.tolist(),
        "conics": k.tolist(),
        "asphere_sag_coefficients_mm": np.zeros((6, 5)).tolist(),
        "axial_distances_mm": distances,
        "stop_distance_mm": stop_d,
        "stop_radius_mm": 25.0 if mode == "epd50" else (21.7 if variant % 2 == 0 else 23.3),
        "clear_radii_mm": radii.tolist(),
        "aperture_mode": mode,
        "provenance": {
            "method": "analytic_thin_power_then_exact_reduced_ABCD_compensation",
            "family": family,
            "bending_variant": variant,
            "thin_initial_element_powers_per_mm": prior.tolist(),
            "fitted_thin_element_powers_per_mm": fit.x.tolist(),
            "focal_target_mm": float(TARGET_EFL),
            "paraxial_EFL_mm": float(efl),
            "paraxial_BFL_mm": float(bfl),
            "actual_image_gap_mm": distances[-1],
            "first_order_matrix_reduced_y_nu": matrix.tolist(),
            "minimum_sampled_material_and_air_axial_clearances_mm": clearance,
            "clearance_radius_samples": 257,
            "fit_nfev": fit.nfev,
            "fit_success": bool(fit.success),
            "real_ray_qualified": False,
            "limits": "Only paraxial power and sampled nonintersection prepared. No claim of NA, image-height, throughput, RGB correction or MTF feasibility.",
        },
    }


def patent_starts():
    source = HERE / "references/US3944337A_examples.json"
    if not source.exists():
        return []
    examples = json.loads(source.read_text())["examples"]
    seeds = []
    for example, glass, mode in itertools.product(examples, ("FDS18W", "TAFD45"), ("epd50", "front_diameter50")):
        scale = TARGET_EFL / float(example["focal_length_mm"])
        materials = ["PMMA", "PC", glass]
        c = np.asarray(example["curvatures_per_mm"], float) / scale
        k = np.asarray(example["conic_constants"], float)
        old = [example["center_thicknesses_mm"][0], example["air_gaps_mm"][0], example["center_thicknesses_mm"][1], example["air_gaps_mm"][1], example["center_thicknesses_mm"][2]]
        old = np.asarray(old) * scale
        # The enlarged opening can exceed an original positive-conic domain.
        # Hyperbolic replacement is an explicit new seed perturbation, not a
        # claim that the enlarged prescription remains the patented design.
        conic_changes = []
        for j in range(4):
            if 1.0 - (1.0 + k[j]) * c[j] ** 2 * 27.0 ** 2 <= 0.05:
                conic_changes.append({"surface": j + 1, "original_scaled_conic": float(k[j]), "seed_conic": -1.0})
                k[j] = -1.0
        scaled_patent_curvatures = c.copy()
        upper = np.array([1.5, 1.5, 0.97 / (19.5 * max(abs(c[4]), abs(c[5])))])
        for j in range(4):
            if k[j] > -1.0 and abs(c[j]) > 0:
                upper[j // 2] = min(upper[j // 2], 0.95 / (27.0 * abs(c[j]) * np.sqrt(1.0 + k[j])))

        def compensate(factors):
            trial = scaled_patent_curvatures * np.repeat(factors, 2)
            try:
                _, ds, _, _, bf, mat, _ = prepare_geometry(materials, trial, k, mode, 1.8 * scale, old)
                return np.array([2500.0 * (-mat[1, 0] - 1.0 / TARGET_EFL), 0.04 * min(0.0, bf - 4.0), 0.005 * max(0.0, sum(ds) - 130.0), *((factors - 1.0) * 0.002)])
            except (ValueError, FloatingPointError):
                return np.array([1e4, 1e4, 1e4, *(factors - 1.0)])

        fit = least_squares(compensate, np.minimum([0.9, 0.7, 0.9], upper * 0.95), bounds=(np.full(3, 0.1), upper), max_nfev=100, ftol=1e-9, xtol=1e-9, gtol=1e-9)
        c = scaled_patent_curvatures * np.repeat(fit.x, 2)
        radii, distances, stop_d, efl, bfl, matrix, clearance = prepare_geometry(materials, c, k, mode, 1.8 * scale, old)
        seeds.append({
            "seed_id": "",
            "materials": materials,
            "curvatures": c.tolist(),
            "conics": k.tolist(),
            "asphere_sag_coefficients_mm": np.zeros((6, 5)).tolist(),
            "axial_distances_mm": distances,
            "stop_distance_mm": stop_d,
            "stop_radius_mm": 25.0 if mode == "epd50" else 21.7,
            "clear_radii_mm": radii.tolist(),
            "aperture_mode": mode,
            "provenance": {
                "method": "public_patent_prescription_scaling_material_substitution_and_aperture_growth",
                "patent": "US3944337A",
                "source_url": "https://patents.google.com/patent/US3944337A/en",
                "example": example["example"],
                "source_data_file": "references/US3944337A_examples.json",
                "linear_scale": float(scale),
                "scaled_patent_curvatures_before_ABCD_compensation_per_mm": scaled_patent_curvatures.tolist(),
                "element_curvature_factors_after_material_aperture_change": fit.x.tolist(),
                "ABCD_compensation_nfev": fit.nfev,
                "material_changes": ["original PMMA-like index to user PMMA nd/Vd", "original polystyrene to user PC nd/Vd", "original n1.734 V51 glass to user " + glass],
                "conic_domain_perturbations": conic_changes,
                "focal_target_mm": float(TARGET_EFL),
                "paraxial_EFL_mm": float(efl),
                "paraxial_BFL_mm": float(bfl),
                "actual_image_gap_mm": distances[-1],
                "first_order_matrix_reduced_y_nu": matrix.tolist(),
                "minimum_sampled_material_and_air_axial_clearances_mm": clearance,
                "clearance_radius_samples": 257,
                "real_ray_qualified": False,
                "limits": "The original patent is about f/1.4, not this aperture/NA. Material, aperture and clearance changes break original correction; full real-ray optimization is required.",
            },
        })
    return seeds


def generate():
    seeds = []
    for glass, pair, position, mode in itertools.product(("FDS18W", "TAFD45"), itertools.product(("PMMA", "PC"), repeat=2), range(3), ("epd50", "front_diameter50")):
        materials = list(pair)
        materials.insert(position, glass)
        for variant, family in enumerate(("negative_glass_positive_plastics", "positive_glass_mixed_plastics")):
            seeds.append(analytic_start(materials, mode, family, variant))
    seeds.extend(patent_starts())
    for i, seed in enumerate(seeds, 1):
        seed["seed_id"] = f"seed_{i:04d}"
    return seeds


if __name__ == "__main__":
    seeds = generate()
    out = HERE / "seeds.json"
    payload = {"schema": "initial-prescription-seeds/1", "status": "unqualified_optimization_starts", "target_EFL_from_rectilinear_mapping_mm": float(TARGET_EFL), "asphere_convention": "z=c*r^2/(1+sqrt(1-(1+k)c^2*r^2))+sum(b_order*(r/25)^order), order4,6,8,10,12; b is mm", "seed_count": len(seeds), "seeds": seeds}
    out.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    focal = np.array([x["provenance"]["paraxial_EFL_mm"] for x in seeds])
    print(json.dumps({"file": str(out), "count": len(seeds), "analytic_count": sum(x["provenance"]["method"].startswith("analytic") for x in seeds), "EFL_min_max_mm": [float(focal.min()), float(focal.max())], "EFL_within1percent": int(np.sum(abs(focal / TARGET_EFL - 1) <= 0.01)), "qualification": "No real-ray qualification"}, indent=2))
