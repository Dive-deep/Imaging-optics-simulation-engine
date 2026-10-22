"""Independent geometric/radiometric candidate audit, not a detector simulator.

All geometric rays retain their original entrance-pupil quadrature weights.
Spot moments are conditional on transmitted rays; throughput and irradiance are
never renormalized to the surviving pupil. RI assumes uniform distant radiance,
a smooth rotational chief-ray image map and a fixed object-side planar STOP.
The optional bare-interface Fresnel model transports scalar unpolarized power,
not polarization/Jones vectors. Every field/wavelength shares one image plane.
"""
from __future__ import annotations

from pathlib import Path
import json
import math
import numpy as np

from model import (MATERIALS, ORDERS, WAVELENGTHS, first_order, geometry_checks,
                   index_normal_line, native_surfaces, trace_design)
from native.raytrace import STATUS

FIELDS = (0., .5, .8, 1.)
THETA_MAX = math.radians(12.)


def _finite(value):
    """Explicit JSON-safe conversion; undefined numeric values remain null."""
    if isinstance(value, dict):
        return {str(k): _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v) for v in value]
    if isinstance(value, np.ndarray):
        return _finite(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    return value


def _moments(xy, weights, valid):
    keep = np.asarray(valid, bool) & np.isfinite(xy).all(axis=-1) & (weights > 0)
    mass = float(np.sum(weights[keep]))
    if mass <= 0:
        return {'centroid_xy_mm': None, 'rms_radius_mm': None,
                'rms_diameter_mm': None, 'covariance_mm2': None,
                'transmitted_quadrature_weight': mass,
                'status': 'no_transmitted_positive_weight_rays'}
    centroid = np.sum(xy[keep] * weights[keep, None], axis=0) / mass
    delta = xy[keep] - centroid
    cov = (delta * weights[keep, None]).T @ delta / mass
    rms = float(np.sqrt(max(0., np.trace(cov))))
    return {'centroid_xy_mm': centroid, 'rms_radius_mm': rms,
            'rms_diameter_mm': 2 * rms, 'covariance_mm2': cov,
            'transmitted_quadrature_weight': mass,
            'status': 'conditional_transmitted_ray_moment'}


def _surface_normal(surface, positions):
    x, y = positions[..., 0], positions[..., 1]
    q = x*x + y*y
    c, k = surface[1:3]
    disc = 1 - (1+k)*c*c*q
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        radial_factor = c / np.sqrt(disc)
        for order, coefficient in zip(ORDERS, surface[3:8]):
            radial_factor = radial_factor + order * coefficient * q**((order-2)//2)
        normal = np.stack([-radial_factor*x, -radial_factor*y, np.ones_like(x)], axis=-1)
        normal /= np.linalg.norm(normal, axis=-1)[..., None]
    return normal


def fresnel_transmission(design, traced):
    """Six bare interfaces, scalar unpolarized power; failed rays have zero T."""
    status = traced['status']
    valid = (status == 0) & traced['image_forward']
    transmission = np.ones_like(status, dtype=float)
    maximum_snell_residual = 0.
    surfaces = native_surfaces(design)
    wavelengths = np.asarray(traced['wavelengths_nm'])
    fields = np.asarray(traced['fields'])
    n = np.ones((len(wavelengths), 7))
    for lens, material in enumerate(design['materials']):
        n[:, 2*lens+1] = index_normal_line(material, wavelengths)
    incoming = np.zeros(status.shape + (3,))
    incoming[..., 1] = np.sin(THETA_MAX*fields)[None, :, None]
    incoming[..., 2] = np.cos(THETA_MAX*fields)[None, :, None]
    for surface_index, surface in enumerate(surfaces):
        normal = _surface_normal(surface, traced['history'][..., surface_index, :])
        outgoing = traced['direction_history'][..., surface_index, :]
        ci = np.abs(np.sum(incoming * normal, axis=-1))
        ct = np.abs(np.sum(outgoing * normal, axis=-1))
        n1, n2 = n[:, surface_index, None, None], n[:, surface_index+1, None, None]
        with np.errstate(invalid='ignore', divide='ignore'):
            rs = ((n1*ci - n2*ct)/(n1*ci + n2*ct))**2
            rp = ((n1*ct - n2*ci)/(n1*ct + n2*ci))**2
            factor = 1 - .5*(rs+rp)
            snell = np.abs(n1*np.sqrt(np.maximum(0., 1-ci*ci)) -
                           n2*np.sqrt(np.maximum(0., 1-ct*ct)))
        factor = np.where(valid & np.isfinite(factor), np.clip(factor, 0., 1.), 0.)
        transmission *= factor
        if np.any(valid):
            maximum_snell_residual = max(maximum_snell_residual, float(np.max(snell[valid])))
        incoming = outgoing
    transmission[~valid] = 0.
    return transmission, {'maximum_surface_snell_sine_residual': maximum_snell_residual,
                          'model': 'product_of_unpolarized_lossless_bare_interface_power_transmissions',
                          'coatings_absorption_polarization_transport': 'not_included'}


def chief_mapping(design, fields=FIELDS, wavelengths=WAVELENGTHS, step_rad=2e-5):
    """Image-map Jacobian with h/h2 convergence evidence, fixed image plane."""
    angles = np.asarray(fields)*THETA_MAX
    sample_angles = np.r_[angles, angles-step_rad, angles+step_rad,
                          angles-step_rad/2, angles+step_rad/2,
                          -step_rad/2, step_rad/2]
    t = trace_design(design, sample_angles/THETA_MAX, wavelengths, xy=[[0., 0.]])
    good = (t['status'][..., 0] == 0) & t['image_forward'][..., 0]
    height = t['image_xy'][..., 0, 1]
    count = len(fields)
    h = (height[:, 2*count:3*count] - height[:, count:2*count])/(2*step_rad)
    h2 = (height[:, 4*count:5*count] - height[:, 3*count:4*count])/step_rad
    axis_slope = (height[:, -1] - height[:, -2])/step_rad
    g = good[:, :count] & good[:, count:2*count] & good[:, 2*count:3*count]
    g &= good[:, 3*count:4*count] & good[:, 4*count:5*count]
    g &= good[:, -1, None] & good[:, -2, None]
    return {'chief_height_mm': height[:, :count], 'derivative_mm_per_rad': h2,
            'derivative_h_mm_per_rad': h, 'derivative_h2_minus_h_mm_per_rad': h2-h,
            'axis_chief_slope_mm_per_rad': axis_slope, 'valid': g,
            'status_counts': {STATUS.get(int(k), str(k)): int(v)
                              for k, v in zip(*np.unique(t['status'], return_counts=True))},
            'step_rad': step_rad,
            'status': 'finite_difference_chief_map_not_full_aberrated_detector_radiometry'}


def _na_metrics(directions, valid, pupil_xy, stop_radius, chief_index=0):
    rim = np.isclose(np.linalg.norm(pupil_xy, axis=-1), stop_radius, rtol=1e-10, atol=1e-9)
    good = np.asarray(valid) & np.isfinite(directions).all(axis=-1)
    axial = np.sqrt(np.sum(directions[..., :2]**2, axis=-1))
    rim_values = axial[rim & good]
    rim_complete = bool(np.any(rim) and np.all(good[rim]))
    result = {'image_medium_index': 1., 'definition': 'n_image*sin(angle_to_optical_z_axis)',
              'rim_ray_count': int(np.sum(rim)), 'valid_rim_ray_count': int(np.sum(rim & good)),
              'complete_rim_transmitted': rim_complete,
              'valid_rim_minimum': float(np.min(rim_values)) if len(rim_values) else None,
              'valid_rim_median': float(np.median(rim_values)) if len(rim_values) else None,
              'valid_rim_maximum': float(np.max(rim_values)) if len(rim_values) else None,
              'conservative_complete_rim_NA': float(np.min(rim_values)) if rim_complete else None,
              'all_valid_sample_maximum_NOT_acceptance_value': float(np.max(axial[good])) if np.any(good) else None}
    if good[chief_index]:
        chief = directions[chief_index]
        angle = np.arccos(np.clip(np.sum(directions*chief, axis=-1), -1., 1.))
        relative = np.sin(angle[rim & good])
        result['chief_angle_to_z_deg'] = math.degrees(math.acos(np.clip(chief[2], -1., 1.)))
        result['rim_sine_angle_to_chief_minimum'] = float(np.min(relative)) if len(relative) else None
        result['rim_sine_angle_to_chief_maximum'] = float(np.max(relative)) if len(relative) else None
    else:
        result['chief_angle_to_z_deg'] = None
    for axis_name, axis in [('tangential', 1), ('sagittal', 0)]:
        lo, hi = int(np.argmin(pupil_xy[:, axis])), int(np.argmax(pupil_xy[:, axis]))
        pair_good = good[lo] and good[hi] and rim[lo] and rim[hi]
        half_angle = .5*math.acos(np.clip(float(np.dot(directions[lo], directions[hi])), -1., 1.)) if pair_good else None
        result[axis_name+'_opposite_rim_cone_NA'] = math.sin(half_angle) if half_angle is not None else None
    return result


def _constraint(passed, measured, requirement, **extra):
    return {'passed': bool(passed), 'measured': measured, 'requirement': requirement, **extra}


def evaluate_design(design, trace_result=None, *, rings=16, arms=48,
                    fields=FIELDS, wavelengths=WAVELENGTHS,
                    image_height_tolerance_mm=.032, geometry_points=1025):
    """Return JSON-safe metrics and explicit per-constraint gates.

    image_height_tolerance_mm=.032 is a declared numeric proximity rule (0.5%
    of6.4), not a user-approved change to the target or manufacturing tolerance.
    This evaluator never moves the image plane or refocuses any field/color.
    """
    t = trace_result if trace_result is not None else trace_design(
        design, fields, wavelengths, rings=rings, arms=arms)
    fields = np.asarray(t['fields'], dtype=float)
    wavelengths = np.asarray(t['wavelengths_nm'], dtype=float)
    weights = np.asarray(t['weights'], dtype=float)
    if np.any(weights < 0) or not np.isclose(weights.sum(), 1., rtol=0., atol=2e-14):
        raise ValueError('Expected original normalized nonnegative entrance pupil quadrature weights')
    if not np.isclose(fields[0], 0.) or not np.isclose(fields[-1], 1.):
        raise ValueError('Evaluation fields must include0 first and1 last')
    if np.count_nonzero(np.isclose(wavelengths, 588.)) != 1:
        raise ValueError('Evaluation wavelengths must include588nm exactly once')
    primary = int(np.flatnonzero(np.isclose(wavelengths, 588.))[0])
    valid = (t['status'] == 0) & t['image_forward'] & np.isfinite(t['image_xy']).all(axis=-1)
    bare, fresnel_evidence = fresnel_transmission(design, t)
    mapping = chief_mapping(design, fields, wavelengths)
    area = math.pi*float(design['stop_radius_mm'])**2
    throughput = np.sum(valid*weights, axis=-1)
    bare_throughput = np.sum(bare*weights, axis=-1)
    irradiance_g = np.full(throughput.shape, np.nan)
    irradiance_b = np.full(throughput.shape, np.nan)
    mapping_valid = mapping['valid'].copy()
    f0 = mapping['axis_chief_slope_mm_per_rad']
    for fi, fraction in enumerate(fields):
        theta = THETA_MAX*fraction
        if abs(theta) < 1e-12:
            factor = 1/f0**2
        else:
            denominator = mapping['chief_height_mm'][:, fi]*mapping['derivative_mm_per_rad'][:, fi]
            factor = np.cos(theta)*np.sin(theta)/denominator
            mapping_valid[:, fi] &= denominator > 0
        mapping_valid[:, fi] &= (f0 > 0) & np.isfinite(factor) & (factor > 0)
        irradiance_g[:, fi] = np.where(mapping_valid[:, fi], area*throughput[:, fi]*factor, np.nan)
        irradiance_b[:, fi] = np.where(mapping_valid[:, fi], area*bare_throughput[:, fi]*factor, np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        ri_g = irradiance_g/irradiance_g[:, [0]]
        ri_b = irradiance_b/irradiance_b[:, [0]]
    by_wavelength = {}
    all_na = []
    for wi, wavelength in enumerate(wavelengths):
        fo = first_order(design, float(wavelength))
        per_field = {}
        for fi, fraction in enumerate(fields):
            good = valid[wi, fi]
            moment = _moments(t['image_xy'][wi, fi], weights, good)
            moment_bare = _moments(t['image_xy'][wi, fi], weights*bare[wi, fi], good)
            na = _na_metrics(t['dir'][wi, fi], good, t['pupil_xy'], design['stop_radius_mm'])
            if fi == 0:
                all_na.append(na['conservative_complete_rim_NA'])
            height = float(t['image_xy'][wi, fi, 0, 1]) if good[0] else None
            theta = THETA_MAX*fraction
            reference = fo['efl_mm']*math.tan(theta)
            slope_reference = float(f0[wi])*math.tan(theta)
            distortion = 100*(height/reference-1) if height is not None and abs(reference)>1e-15 else (0. if good[0] and fi==0 else None)
            distortion_slope = 100*(height/slope_reference-1) if height is not None and abs(slope_reference)>1e-15 else (0. if good[0] and fi==0 else None)
            status_counts = {STATUS.get(int(k), str(k)): int(v)
                             for k, v in zip(*np.unique(t['status'][wi, fi], return_counts=True))}
            failure_surface = {str(int(k)+1): int(v) for k,v in
                zip(*np.unique(t['failed_surface'][wi, fi][t['status'][wi, fi]!=0], return_counts=True))}
            per_field[str(float(fraction))] = {
                'field_angle_deg': float(fraction*12), 'spot_geometric': moment,
                'spot_uncoated_fresnel': moment_bare, 'image_space_NA': na,
                'chief_height_mm': height, 'chief_image_xy_mm': t['image_xy'][wi, fi, 0] if good[0] else None,
                'rectilinear_reference_height_paraxial_EFL_mm': reference,
                'distortion_percent_paraxial_EFL': distortion,
                'distortion_percent_axis_chief_slope': distortion_slope,
                'entrance_pupil_area_mm2': area,
                'geometric_transmitted_pupil_fraction': float(throughput[wi, fi]),
                'uncoated_power_transmission_fraction': float(bare_throughput[wi, fi]),
                'RI_geometric': float(ri_g[wi, fi]), 'RI_uncoated_fresnel': float(ri_b[wi, fi]),
                'irradiance_per_unit_input_radiance_geometric_sr': float(irradiance_g[wi, fi]),
                'irradiance_per_unit_input_radiance_uncoated_sr': float(irradiance_b[wi, fi]),
                'chief_mapping_valid': bool(mapping_valid[wi, fi]),
                'ray_count_including_chief_rim': len(weights),
                'positive_weight_ray_count': int(np.sum(weights>0)),
                'invalid_ray_count_including_chief_rim': int(np.sum(~good)),
                'native_status_counts': status_counts,
                'failed_surface_counts_1_based': failure_surface,
                'backward_image_projection_count': int(np.sum((t['status'][wi, fi]==0)&~t['image_forward'][wi, fi])),
                'vignetting_status': 'none_on_positive_weight_quadrature' if throughput[wi, fi]>=1-1e-14 else 'vignetted_or_failed',
                'RI_status': 'conditional_local_chief_mapping_and_sampled_pupil' if mapping_valid[wi, fi] else 'undefined_chief_mapping'}
        by_wavelength[str(float(wavelength))] = {'first_order': fo, 'fields': per_field}
    # Equal incident spectral radiance weights. No RGB PSF recentering here.
    spectral = np.full(len(wavelengths), 1/len(wavelengths))
    poly_g = np.sum(spectral[:, None]*irradiance_g, axis=0)
    poly_b = np.sum(spectral[:, None]*irradiance_b, axis=0)
    polychromatic = {}
    for fi, fraction in enumerate(fields):
        poly_weights = (spectral[:, None]*np.broadcast_to(weights, (len(wavelengths), len(weights)))).ravel()
        xy = t['image_xy'][:, fi].reshape(-1, 2)
        good = valid[:, fi].ravel()
        polychromatic[str(float(fraction))] = {
            'spot_geometric_common_RGB_centroid': _moments(xy, poly_weights, good),
            'spot_uncoated_common_RGB_centroid': _moments(xy, poly_weights*bare[:, fi].ravel(), good),
            'RI_geometric': float(poly_g[fi]/poly_g[0]) if poly_g[0]>0 else None,
            'RI_uncoated_fresnel': float(poly_b[fi]/poly_b[0]) if poly_b[0]>0 else None}
    geom = geometry_checks(design, points=geometry_points)
    materials = design['materials']
    plastic = [i for i,m in enumerate(materials) if m in ('PMMA','PC')]
    glass = [i for i,m in enumerate(materials) if m in ('FDS18W','TAFD45')]
    glass_standard = all(float(design['conics'][j])==0 and
                         np.all(np.asarray(design['asphere_sag_coefficients_mm'][j])==0)
                         for lens in glass for j in (2*lens,2*lens+1))
    axis_na = all_na[primary]
    primary_full = by_wavelength['588.0']['fields']['1.0']['chief_height_mm']
    image_error = primary_full-6.4 if primary_full is not None else None
    epd = 2*float(design['stop_radius_mm'])
    first_diameter = 2*float(design['clear_radii_mm'][0])
    nominal_pupil_fully_transmitted_onaxis = bool(np.all(valid[:,0]))
    actual_epd_supported = epd if nominal_pupil_fully_transmitted_onaxis else None
    acceptable_aperture = (math.isclose(epd,50.,rel_tol=0.,abs_tol=1e-8) and nominal_pupil_fully_transmitted_onaxis) or math.isclose(first_diameter,50.,rel_tol=0.,abs_tol=1e-8)
    fo_primary = first_order(design, 588.)
    constraints = {
        'three_separated_elements_two_plastics_one_glass': _constraint(len(materials)==3 and len(plastic)==2 and len(glass)==1, materials, 'twoPMMA/PC plusoneFDS18W/TAFD45'),
        'glass_surfaces_standard': _constraint(glass_standard, glass_standard, 'glass k=0 andA4..A12=0'),
        'only_even_orders4_through12': _constraint(np.asarray(design['asphere_sag_coefficients_mm']).shape==(6,5), 'A4,A6,A8,A10,A12', 'no other aspheric orders'),
        'aperture_alternative50mm': _constraint(acceptable_aperture, {'nominal_STOP_diameter_mm':epd,'supported_unvignetted_EPD_mm':actual_epd_supported,'first_surface_clear_diameter_mm':first_diameter}, 'EPD50fullytransmittedonaxis or first_lens_object_side_clear_diameter50'),
        'center_thickness_strictly_over2mm': _constraint(all(x>2 for x in geom['center_thickness_mm']), geom['center_thickness_mm'], '>2mm each'),
        'center_thickness_at_most20mm': _constraint(all(x<=20 for x in geom['center_thickness_mm']), geom['center_thickness_mm'], '<=20mm each; exact floating-point comparison, no positive tolerance'),
        'axis_air_gaps_at_least1mm': _constraint(all(x>=1 for x in geom['axis_air_gaps_mm']), geom['axis_air_gaps_mm'], '>=1mm each'),
        'vertex_BFL_at_least3mm': _constraint(design['axial_distances_mm'][5]>=3, design['axial_distances_mm'][5], 'last_surface_vertex_to_fixed_image_plane>=3mm'),
        'positive_finite_paraxial_EFL': _constraint(math.isfinite(fo_primary['efl_mm']) and fo_primary['efl_mm']>0, fo_primary['efl_mm'], 'finitepositiveimageformingpower'),
        'paraxial_BFL_at_least3mm': _constraint(math.isfinite(fo_primary['paraxial_bfl_mm']) and fo_primary['paraxial_bfl_mm']>=3, fo_primary['paraxial_bfl_mm'], 'finiteparaxialbackfocaldistance>=3mm', interpretation='additionalconservativeBFLdefinitioncheck;distinctfromactualfixedimageplaneBFL'),
        'real_surface_sag_domain': _constraint(all(x>0 for x in geom['minimum_conic_domains']), geom['minimum_conic_domains'], 'strictlypositiveconicradicand acrossclearaperture'),
        'no_adjacent_surface_crossing_sampled': _constraint(all(x>0 for x in geom['minimum_adjacent_surface_separation_mm']), geom['minimum_adjacent_surface_separation_mm'], '>0mm overcommonclearaperture', radial_samples_per_interval=geometry_points, continuous_proof=False),
        'stop_plane_before_entire_first_clear_surface': _constraint(geom['stop_front_min_clearance_mm']>0, geom['stop_front_min_clearance_mm'], 'positiveclearance'),
        'image_plane_after_entire_last_clear_surface': _constraint(geom['image_back_min_clearance_mm']>0, geom['image_back_min_clearance_mm'], 'positiveclearance'),
        'onaxis_588nm_NA_at_least0p7': _constraint(axis_na is not None and axis_na>=.7, axis_na, 'allnominalrimraysvalidandminimumrimNA>=0.7'),
        'fullfield_588nm_height_numerically_close_to6p4mm': _constraint(image_error is not None and abs(image_error)<=image_height_tolerance_mm, primary_full, 'target6.4mm', residual_mm=image_error, numerical_proximity_tolerance_mm=image_height_tolerance_mm, user_approved_relaxation=False),
        'all_requested_sampled_rays_transmitted': _constraint(bool(np.all(valid)), {'invalid_count':int(np.sum(~valid))}, 'no deletedfailedorchiefrimrays;candidatequalificationpolicy_notanextrauserthroughputspec')}
    hard_names = [k for k in constraints if k!='all_requested_sampled_rays_transmitted']
    numerical_requested = all(constraints[k]['passed'] for k in hard_names)
    fully_sampled = numerical_requested and constraints['all_requested_sampled_rays_transmitted']['passed']
    return _finite({
        'schema': 'candidate-geometric-radiometric-audit/1',
        'candidate_id': design.get('candidate_id',design.get('seed_id')),
        'status': 'numerically_qualified_initial_candidate' if fully_sampled else 'unqualified_initial_candidate',
        'not_global_optimum_or_Zemax_validation': True,
        'all_declared_numeric_constraints_passed': numerical_requested,
        'all_declared_numeric_constraints_and_sampled_rays_passed': fully_sampled,
        'constraint_results': constraints,
        'failed_constraints': [name for name,c in constraints.items() if not c['passed']],
        'aperture': {'EPD_mm': actual_epd_supported, 'nominal_STOP_diameter_mm': epd, 'first_surface_clear_diameter_mm': first_diameter,
                     'first_lens_maximum_surface_clear_diameter_mm':2*max(design['clear_radii_mm'][:2]),
                     'EPD_definition': 'physical_stop_diameter;no_powered_surface_before_stop',
                     'onaxis_nominal_pupil_fully_transmitted_allRGB':nominal_pupil_fully_transmitted_onaxis,
                     'onaxis_effective_area_equivalent_diameter_RGB_mm':epd*np.sqrt(throughput[:,0]),
                     'effective_area_equivalent_diameter_is_not_a_pupil_boundary_measurement':True,
                     'stop_distance_before_first_vertex_mm':design['stop_distance_mm'],
                     'stop_distance_less_than_positive_EFL_preference_met': bool(fo_primary['efl_mm']>0 and design['stop_distance_mm']<fo_primary['efl_mm'])},
        'fixed_image_plane_z_mm': float(sum(design['axial_distances_mm'])),
        'vertex_BFL_mm': float(design['axial_distances_mm'][5]),
        'first_order_primary588nm': fo_primary,
        'geometry_dense_sample': geom,
        'by_wavelength_nm': by_wavelength,
        'polychromatic_equal_incident_radiance_weights': spectral,
        'polychromatic_fields': polychromatic,
        'chief_mapping_evidence': mapping,
        'fresnel_evidence': fresnel_evidence,
        'NA_rgb_onaxis_conservative_complete_rim': all_na,
        'RI_model': {'formula':'E/L=A_eff*cos(theta)*sin(theta)/(r(theta)*dr/dtheta);axislimit=A_eff/f0^2',
                     'pupil_measure':'originalentranceSTOPplanearea;positiveweightquadrature;chief/rimzeroarea',
                     'object':'uniformdistantincoherentradiance;equalRGBspectralradiance',
                     'image':'onefixedflatplane;localchiefmap;aberratedPSFredistributionnotintegrated',
                     'status':'conditionalgeometricradiometry;notfullOpticStudioexitpupilsolidangleintegration'},
        'sampling': {'rings':rings,'arms':arms,'ray_count_per_field_wavelength':len(weights),
                     'geometry_radial_samples':geometry_points},
        'all_requested_native_status_counts':{STATUS.get(int(k),str(k)):int(v) for k,v in zip(*np.unique(t['status'],return_counts=True))},
        'all_invalid_ray_count_including_chief_rim':int(np.sum(~valid))})


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    design=json.loads(args.input.read_text())
    design=design.get('design',design)
    result=evaluate_design(design)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'failed_constraints':result['failed_constraints']},indent=2))
