"""Geometrical and wavefront diagnostics (calculation_spec_ko.md G8--G12).

All lengths and OPD are mm; directions use radians.  Real rays, not an EFL
substitution, determine intercepts.  T is the radial field meridian and S its
perpendicular; on axis T=Y,S=X.  These labels are not a vendor FFT axis claim.

Analyses preserve the chosen IMA, quadrature measure, and reference.  Only
``common_geometric_refocus`` moves one common plane, and that optimizes geometric
centroid RMS, not MTF.  Zernike fitting never mutates the input OPD/pupil.
"""
from __future__ import annotations

import math
from typing import Callable, Sequence
import numpy as np
from scipy.special import eval_jacobi

from .errors import InputError, MissingDataError, NumericalError, ModelApplicabilityError
from .models import Field
from .provenance import require_bundle_match, validate_bundle_integrity
from .settings import EngineSettings


_PHYSICAL_LOSSES = frozenset({'aperture_blocked', 'escaped', 'tir', 'TIR'})


def tangential_sagittal_basis(field: Field) -> tuple[np.ndarray, np.ndarray]:
    """Global XY radial field basis, returned as ``(T,S)`` unit vectors."""
    v = np.array([field.x, field.y], dtype=float)
    if field.kind == 'angle':
        v = np.tan(v)
    r = np.linalg.norm(v)
    if r == 0:
        return np.array([0., 1.]), np.array([1., 0.])
    t = v / r
    return t, np.array([t[1], -t[0]])


def _weights(weights, n):
    w = np.asarray(weights, dtype=float)
    if w.shape != (n,) or not np.all(np.isfinite(w)) or np.any(w < 0) or w.sum() <= 0:
        raise InputError('Weights must be a finite nonnegative N-vector with positive sum')
    return w


def _valid_bundle(bundle):
    """Return geometric membership, separately from the quadrature weights.

    Successfully traced chief/rim diagnostics with zero weight belong in a
    spot plot and sampled maximum. They do not contribute to pupil moments.
    """
    n = len(bundle.traces)
    w = _weights(bundle.weights, n)
    if np.shape(bundle.pupil_xy) != (n, 2):
        raise InputError('Bundle pupil_xy must have shape (N,2)')
    for tr in bundle.traces:
        if not tr.success and tr.status not in _PHYSICAL_LOSSES:
            raise NumericalError(f'Unresolved ray cannot be counted as vignetting: {tr.status}: {tr.message}')
    mask = np.array([t.success for t in bundle.traces], dtype=bool)
    if not np.any(mask & (w > 0)):
        raise ModelApplicabilityError('No successful positive-weight rays for this metric')
    return mask, w


def _last(tr):
    if not tr.success:
        if tr.status not in _PHYSICAL_LOSSES:
            raise NumericalError(f'Unresolved ray: {tr.status}: {tr.message}')
        raise ModelApplicabilityError(f'Physical ray unavailable: {tr.status}')
    if len(tr.positions) == 0 or len(tr.directions) == 0:
        raise NumericalError('Successful trace contains no final ray')
    p, d = np.asarray(tr.last_position, float), np.asarray(tr.last_direction, float)
    if p.shape != (3,) or d.shape != (3,) or not np.all(np.isfinite([p, d])):
        raise NumericalError('Nonfinite or malformed final ray')
    if not np.isclose(np.linalg.norm(d), 1.):
        raise NumericalError('Final ray direction is not unit length')
    return p, d


def _image_z(system, image_z):
    z = system.image_z if image_z is None else image_z
    if z is None:
        raise MissingDataError('A finite image_z is required; use angular_spot for afocal output')
    if not np.isfinite(z):
        raise InputError('image_z must be finite')
    return float(z)


def ray_plane_intercept(trace, image_z, *, direction_tolerance=None):
    """Intersect final ray with a global-Z plane; report signed virtual extension.

    A negative distance is valid geometric extrapolation but is not a physical
    screen hit. This function does not add screen clipping or ray power.
    """
    tol = EngineSettings().ray.direction_tol if direction_tolerance is None else direction_tolerance
    p, d = _last(trace)
    if abs(d[2]) <= tol:
        raise ModelApplicabilityError('Ray is parallel to the image plane')
    distance = (float(image_z) - p[2]) / d[2]
    return (p + distance * d)[:2], float(distance), bool(distance < 0)


def _moments(x, w, reference):
    # Keep the integration subset unchanged by geometric-only samples. In
    # particular, adding a zero-weight rim must not perturb the centroid/RMS.
    positive = w > 0
    quadrature_x, quadrature_w = x[positive], w[positive]
    wn = quadrature_w / quadrature_w.sum()
    centroid = wn @ quadrature_x
    dx = quadrature_x - centroid
    covariance = (dx * wn[:, None]).T @ dx
    residual = x - reference
    return {
        'centroid': centroid,
        'covariance': covariance,
        'rms_radius': float(np.sqrt(np.sum(wn * np.sum(residual[positive]**2, axis=1)))),
        'geo_radius_sampled': float(np.max(np.linalg.norm(residual, axis=1))),
    }


def spot_diagram(system, bundle, image_z=None, *, reference='centroid', allow_virtual=True,
                 settings=None):
    """Weighted geometric spot on one explicit global-Z image plane.

    Weight is the bundle's pupil quadrature, not LED radiance or traced power.
    All successful coordinates, including zero-weight chief/rim samples,
    retain their source ray indices. Only positive quadrature weights enter
    the centroid/RMS; the geometric radius is a sampled, not continuous-pupil,
    maximum. Physical clipping is allowed; unresolved rays raise rather than
    silently improve spot RMS.
    """
    cfg = settings or EngineSettings()
    require_bundle_match(system, bundle, cfg, allow_unverified=True)
    z = _image_z(system, image_z)
    mask, allw = _valid_bundle(bundle)
    ids = np.flatnonzero(mask)
    intercepts = [ray_plane_intercept(bundle.traces[i], z, direction_tolerance=cfg.ray.direction_tol) for i in ids]
    x = np.array([v[0] for v in intercepts])
    virtual = np.array([v[2] for v in intercepts])
    if np.any(virtual) and not allow_virtual:
        raise ModelApplicabilityError('Requested image plane requires virtual ray extension')
    w = allw[ids]
    positive = w > 0
    centroid = (w[positive] / w[positive].sum()) @ x[positive]
    chief = None
    if bundle.chief.success:
        chief = ray_plane_intercept(bundle.chief, z, direction_tolerance=cfg.ray.direction_tol)[0]
    elif bundle.chief.status not in _PHYSICAL_LOSSES:
        raise NumericalError(f'Chief ray unresolved: {bundle.chief.status}')
    if reference == 'centroid':
        ref = centroid
    elif reference == 'chief':
        if chief is None:
            raise ModelApplicabilityError('Chief-referenced spot requires a transmitted chief')
        ref = chief
    else:
        raise InputError('reference must be centroid or chief')
    t, s = tangential_sagittal_basis(bundle.field)
    out = _moments(x, w, ref)
    out.update(status='ok', coordinates_mm=x, ray_indices=ids, weights=w.copy(),
               weight_measure='bundle_pupil_quadrature', quadrature_mask=positive,
               geometric_radius_scope='all_successful_traced_samples', image_z_mm=z,
               reference=reference, reference_xy_mm=ref.copy(), chief_xy_mm=chief,
               pupil_xy=np.asarray(bundle.pupil_xy)[ids].copy(), virtual_intercepts=virtual,
               transmitted_quadrature_fraction=float(w.sum() / allw.sum()),
               basis_T=t, basis_S=s, coordinates_TS_mm=(x-ref) @ np.column_stack([t, s]),
               sigma_T_mm=float(np.sqrt(max(0., t @ out['covariance'] @ t))),
               sigma_S_mm=float(np.sqrt(max(0., s @ out['covariance'] @ s))),
               wavelength_nm=bundle.wavelength_nm)
    return out


def _chief_frame(chief):
    _, dc = _last(chief)
    e1 = np.array([1., 0., 0.]) - dc[0] * dc
    if np.linalg.norm(e1) < np.sqrt(np.finfo(float).eps):
        e1 = np.array([0., 1., 0.]) - dc[1] * dc
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(dc, e1)
    return dc, e1, e2


def _angular_coordinates(directions, dc, e1, e2):
    forward = directions @ dc
    if np.any(forward <= 0):
        raise ModelApplicabilityError('Angular metric requires a common forward chief hemisphere')
    return np.column_stack([np.arctan2(directions @ e1, forward), np.arctan2(directions @ e2, forward)])


def angular_spot(bundle, *, reference='centroid'):
    """Chief-local atan2 spot [rad], with all successful geometric samples.

    Zero-weight chief/rim samples affect the plot and sampled maximum only.
    The chart radius is not a great-circle angle on the direction sphere.
    """
    validate_bundle_integrity(bundle, allow_unverified=True)
    mask, weights = _valid_bundle(bundle)
    dc, e1, e2 = _chief_frame(bundle.chief)
    ids = np.flatnonzero(mask)
    directions = np.array([_last(bundle.traces[i])[1] for i in ids])
    angles = _angular_coordinates(directions, dc, e1, e2)
    w = weights[ids]
    positive = w > 0
    centroid = (w[positive] / w[positive].sum()) @ angles[positive]
    if reference not in ('centroid', 'chief'):
        raise InputError('reference must be centroid or chief')
    ref = centroid if reference == 'centroid' else np.zeros(2)
    out = _moments(angles, w, ref)
    out.update(status='ok', coordinates_rad=angles, units='rad', reference=reference,
               ray_indices=ids, weights=w.copy(), quadrature_mask=positive,
               geometric_radius_scope='all_successful_traced_samples',
               basis_x=e1, basis_y=e2, chief_direction=dc.copy(),
               chart='atan2(d.e1,d.dc), atan2(d.e2,d.dc)', wavelength_nm=bundle.wavelength_nm)
    return out


def _trace_at(system, field, pupil_xy, wavelength_nm, settings):
    from .geometry import aim_ray, trace_ray
    ray = aim_ray(system, field, np.asarray(pupil_xy), wavelength_nm, settings)
    return trace_ray(system, ray, wavelength_nm, settings)


def chief_mapping(system, fields, wavelength_nm=None, *, image_z=None, angular=False, settings=None):
    """Real chief map. Afocal coordinates use fixed global atan2(dx,dz), atan2(dy,dz)."""
    cfg = settings or EngineSettings()
    wavelength = cfg.wavelength(wavelength_nm)
    fields = tuple(fields)
    z = None if angular else _image_z(system, image_z)
    xy = np.full((len(fields), 2), np.nan)
    directions = np.full((len(fields), 3), np.nan)
    statuses, virtual = [], []
    for i, field in enumerate(fields):
        tr = _trace_at(system, field, (0., 0.), wavelength, cfg)
        statuses.append(tr.status)
        if not tr.success:
            if tr.status not in _PHYSICAL_LOSSES:
                raise NumericalError(f'Chief map unresolved: {tr.status}')
            virtual.append(False)
            continue
        _, d = _last(tr)
        directions[i] = d
        if angular:
            if d[2] <= cfg.ray.direction_tol:
                raise ModelApplicabilityError('Global angular map requires a forward hemisphere')
            xy[i] = np.arctan2(d[:2], d[2])
            virtual.append(False)
        else:
            xy[i], _, vi = ray_plane_intercept(tr, z, direction_tolerance=cfg.ray.direction_tol)
            virtual.append(vi)
    return dict(coordinates=xy, directions=directions, fields=fields, statuses=tuple(statuses),
                units='rad' if angular else 'mm', image_z_mm=z, wavelength_nm=wavelength,
                virtual_intercepts=np.array(virtual), reference='chief',
                chart='global_atan2' if angular else 'global_XY')


def map_distortion(coordinates, actual_mapping, reference_jacobian, *, origin_coordinate=(0., 0.),
                   origin_mapping=(0., 0.), nominal_mapping=None):
    """Distortion against an explicitly supplied local linear map.

    Signed radial percent projects map residual onto the reference-map radius.
    At the expansion origin its defined local-derivative limit is zero. At any
    other zero reference radius it is undefined (NaN), not forced to zero.
    Raw maps and optional nominal differences are always retained.
    """
    q, g = np.asarray(coordinates, float), np.asarray(actual_mapping, float)
    j = np.asarray(reference_jacobian, float)
    q0, g0 = np.asarray(origin_coordinate, float), np.asarray(origin_mapping, float)
    if q.ndim != 2 or q.shape[1] != 2 or g.shape != q.shape or j.shape != (2, 2) or q0.shape != (2,) or g0.shape != (2,):
        raise InputError('Map coordinates and mappings must be Nx2; Jacobian 2x2')
    if not all(np.all(np.isfinite(v)) for v in (q, g, j, q0, g0)):
        raise InputError('Distortion needs finite transmitted chief mappings')
    reference = g0 + (q-q0) @ j.T
    error = g-reference
    radial = reference-g0
    norm2 = np.einsum('ij,ij->i', radial, radial)
    percent = np.full(len(q), np.nan)
    nonzero = norm2 > 0
    percent[nonzero] = 100 * np.einsum('ij,ij->i', error[nonzero], radial[nonzero]) / norm2[nonzero]
    origin = np.all(q == q0, axis=1)
    percent[origin & (norm2 == 0)] = 0.
    out = dict(actual_mapping=g.copy(), reference_mapping=reference, residual=error,
               radial_percent=percent, origin_limit_used=origin & (norm2 == 0),
               reference_jacobian=j.copy(), origin_coordinate=q0.copy(), origin_mapping=g0.copy())
    if nominal_mapping is not None:
        nominal = np.asarray(nominal_mapping, float)
        if nominal.shape != g.shape or not np.all(np.isfinite(nominal)):
            raise InputError('Nominal mapping must match actual finite mapping shape')
        out['nominal_mapping'] = nominal.copy()
        out['nominal_change'] = g-nominal
    return out


def distortion(system, fields, wavelength_nm=None, *, image_z=None, angular=False,
               angle_reference='rectilinear', reference_field=None, derivative_step=None,
               convergence_rtol=1e-5, nominal_mapping=None, settings=None):
    """Real-ray distortion with a converged small-field chief Jacobian.

    Infinite fields use tan(angle) for rectilinear or angle for f-theta. Finite
    fields use object height. This never substitutes EFL for a defocused map.
    """
    cfg = settings or EngineSettings()
    fields = tuple(fields)
    if not fields or len({f.kind for f in fields}) != 1:
        raise InputError('Provide nonempty fields of one kind')
    if angle_reference not in ('rectilinear', 'f_theta'):
        raise InputError('angle_reference must be rectilinear or f_theta')
    kind = fields[0].kind
    ref = reference_field or Field(kind=kind)
    if ref.kind != kind:
        raise InputError('Reference field kind must match fields')
    def chart(f):
        v = np.array([f.x, f.y])
        return np.tan(v) if kind == 'angle' and angle_reference == 'rectilinear' else v
    def from_chart(v):
        a = np.arctan(v) if kind == 'angle' and angle_reference == 'rectilinear' else v
        return Field(float(a[0]), float(a[1]), kind)
    q0 = chart(ref)
    step = cfg.ray.derivative_step if derivative_step is None else float(derivative_step)
    if not np.isfinite(step) or step <= 0 or convergence_rtol <= 0:
        raise InputError('Derivative step and convergence tolerance must be positive')
    q = np.array([chart(f) for f in fields])
    mapped = chief_mapping(system, fields, wavelength_nm, image_z=image_z, angular=angular, settings=cfg)
    g0 = chief_mapping(system, [ref], wavelength_nm, image_z=image_z, angular=angular, settings=cfg)['coordinates'][0]
    jacobians = []
    for h in (step, step/2, step/4):
        perturbed = [from_chart(q0 + sign*h*np.eye(2)[axis]) for axis in range(2) for sign in (-1, 1)]
        gm = chief_mapping(system, perturbed, wavelength_nm, image_z=image_z, angular=angular, settings=cfg)['coordinates']
        jacobians.append(np.column_stack([(gm[1]-gm[0])/(2*h), (gm[3]-gm[2])/(2*h)]))
    j = jacobians[-1]
    if not np.all(np.isfinite(j)):
        raise ModelApplicabilityError('Reference chief neighborhood is clipped')
    err = np.linalg.norm(jacobians[-1]-jacobians[-2]) / max(np.linalg.norm(j), np.finfo(float).tiny)
    if err > convergence_rtol:
        raise NumericalError(f'Distortion reference derivative did not converge: relative change {err:g}')
    out = map_distortion(q, mapped['coordinates'], j, origin_coordinate=q0, origin_mapping=g0,
                         nominal_mapping=nominal_mapping)
    out.update(status='ok', chief_map=mapped, coordinate_chart='object_height' if kind == 'height' else angle_reference,
               derivative_step=step, derivative_relative_change=float(err))
    return out


def grid_mapping(system, x_coordinates, y_coordinates, *, field_kind='height', wavelength_nm=None,
                 image_z=None, angular=False, include_spots=False, settings=None):
    """Map every source grid node to a common image plane/directional chart.

    Include-spots traces full bundles at those SAME nodes and plane. Rows/columns
    can be connected by callers; grid drawing and display scaling are not physics.
    """
    from .geometry import trace_bundle
    x, y = np.asarray(x_coordinates, float), np.asarray(y_coordinates, float)
    if x.ndim != 1 or y.ndim != 1 or not len(x) or not len(y) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise InputError('Grid axes must be nonempty finite 1-D arrays')
    xx, yy = np.meshgrid(x, y)
    fields = [Field(a, b, field_kind) for a, b in zip(xx.ravel(), yy.ravel())]
    out = chief_mapping(system, fields, wavelength_nm, image_z=image_z, angular=angular, settings=settings)
    out['grid_shape'] = xx.shape
    out['source_grid'] = np.stack([xx, yy], axis=-1)
    out['mapped_grid'] = out['coordinates'].reshape(*xx.shape, 2)
    if include_spots:
        spots = []
        for f in fields:
            bundle = trace_bundle(system, f, wavelength_nm, settings)
            spots.append(angular_spot(bundle) if angular else spot_diagram(system, bundle, image_z, settings=settings))
        out['spots'] = spots
    return out


def ray_fans(system, field, wavelength_nm=None, *, image_z=None, pupil_coordinates=None,
             points=21, angular=False, settings=None):
    """Both components of real T/S STOP diameter fans, chief-referenced.

    Angular mode records the same chief-local atan2 XY chart as ``angular_spot``.
    Its T/S components use exact direction projections onto the 3-D meridional
    tangent frame; finite atan2 angles are never linearly rotated. It requires
    no IMA and does not claim a vendor angular convention.
    """
    cfg = settings or EngineSettings()
    wavelength = cfg.wavelength(wavelength_nm)
    z = None if angular else _image_z(system, image_z)
    q = np.linspace(-1., 1., points) if pupil_coordinates is None else np.asarray(pupil_coordinates, float)
    if q.ndim != 1 or not len(q) or not np.all(np.isfinite(q)) or np.any(abs(q) > 1):
        raise InputError('Fan coordinates must be finite normalized STOP coordinates in [-1,1]')
    t, s = tangential_sagittal_basis(field)
    basis = np.column_stack([t, s])
    chief = _trace_at(system, field, (0., 0.), wavelength, cfg)
    if angular:
        dc, e1, e2 = _chief_frame(chief)
        t3 = np.array([*t, 0.]); t3 -= (t3 @ dc)*dc
        if np.linalg.norm(t3) <= cfg.ray.direction_tol:
            raise ModelApplicabilityError('Angular fan meridional tangent frame is degenerate')
        t3 /= np.linalg.norm(t3)
        s3 = np.cross(t3, dc)
        xc = np.zeros(2)
    else:
        xc = ray_plane_intercept(chief, z, direction_tolerance=cfg.ray.direction_tol)[0]
    result = dict(pupil_q=q.copy(), basis_T=t, basis_S=s, chief_xy_mm=xc,
                  image_z_mm=z, wavelength_nm=wavelength, reference='chief', units='rad' if angular else 'mm')
    if angular:
        result.pop('chief_xy_mm')
        result.update(chief_direction=dc.copy(), angular_basis_x=e1, angular_basis_y=e2,
                      angular_basis_T=t3, angular_basis_S=s3,
                      angular_chart='atan2(d.e1,d.dc), atan2(d.e2,d.dc)',
                      TS_projection='exact_atan2_direction_projection_in_meridional_tangent_frame')
    for name, vector in (('T', t), ('S', s)):
        values = np.full((len(q), 2), np.nan)
        chart_xy = np.full((len(q), 2), np.nan) if angular else None
        statuses = []
        for i, qi in enumerate(q):
            tr = _trace_at(system, field, qi*vector, wavelength, cfg)
            statuses.append(tr.status)
            if tr.success:
                if angular:
                    direction = _last(tr)[1][None, :]
                    chart_xy[i] = _angular_coordinates(direction, dc, e1, e2)[0]
                    values[i] = _angular_coordinates(direction, dc, t3, s3)[0]
                else:
                    xy = ray_plane_intercept(tr, z, direction_tolerance=cfg.ray.direction_tol)[0]
                    values[i] = (xy-xc) @ basis
            elif tr.status not in _PHYSICAL_LOSSES:
                raise NumericalError(f'Fan ray unresolved: {tr.status}')
        result[f'{name}_fan_TS_rad' if angular else f'{name}_fan_TS_mm'] = values
        if angular:
            result[f'{name}_fan_XY_rad'] = chart_xy
        result[f'{name}_statuses'] = tuple(statuses)
    return result


def coma_from_fan(pupil_q, transverse_values, *, pairing_atol=None):
    """Even/odd paired fan descriptors, NOT Seidel coma coefficients.

    Clipped values stay NaN. For a fan Δ(q), even=(Δ(q)+Δ(-q))/2,
    odd=(Δ(q)-Δ(-q))/2. Supports both transverse components.
    """
    q, v = np.asarray(pupil_q, float), np.asarray(transverse_values, float)
    if q.ndim != 1 or v.shape[0] != len(q) or not np.all(np.isfinite(q)):
        raise InputError('Fan coordinate and value lengths must match')
    tol = np.finfo(float).eps*32 if pairing_atol is None else float(pairing_atol)
    indices = []
    for x in q:
        matches = np.flatnonzero(np.abs(q+x) <= tol)
        if len(matches) != 1:
            raise InputError('Fan must have one unambiguous negative-coordinate pair for every sample')
        indices.append(matches[0])
    return dict(pupil_q=q.copy(), even=(v+v[indices])/2, odd=(v-v[indices])/2,
                label='paired_real_ray_fan_descriptor_not_Seidel')


def longitudinal_spherical_aberration(system, wavelength_nm=None, *, pupil_coordinates=None,
                                     points=21, reference_focus_z=None, settings=None):
    """Axial meridional ray-axis crossings relative to an explicit reference.

    The axial ray's crossing is the finite paraxial conjugate limit. A custom
    reference translates LSA values; it never changes those physical crossings.
    """
    from .paraxial import paraxial_properties
    cfg = settings or EngineSettings()
    wavelength = cfg.wavelength(wavelength_nm)
    kind = 'height' if system.object_z is not None else 'angle'
    field = Field(kind=kind)
    paraxial_focus_z = paraxial_properties(system, wavelength, cfg)['paraxial_image_z']
    if paraxial_focus_z is None or not np.isfinite(paraxial_focus_z):
        raise ModelApplicabilityError('LSA requires a finite paraxial conjugate focus')
    if reference_focus_z is None:
        reference_focus_z = paraxial_focus_z
    if reference_focus_z is None or not np.isfinite(reference_focus_z):
        raise ModelApplicabilityError('LSA requires a finite paraxial conjugate reference')
    q = np.linspace(0., 1., points) if pupil_coordinates is None else np.asarray(pupil_coordinates, float)
    if q.ndim != 1 or not len(q) or not np.all(np.isfinite(q)) or np.any(abs(q) > 1):
        raise InputError('LSA pupil coordinate must lie in [-1,1]')
    crossing = np.full(len(q), np.nan)
    statuses = []
    for i, qi in enumerate(q):
        if qi == 0:
            crossing[i] = paraxial_focus_z
            statuses.append('paraxial_limit')
            continue
        tr = _trace_at(system, field, (0., qi), wavelength, cfg)
        if not tr.success:
            if tr.status not in _PHYSICAL_LOSSES:
                raise NumericalError(f'LSA ray unresolved: {tr.status}')
            statuses.append(tr.status)
            continue
        p, d = _last(tr)
        if abs(d[1]) <= cfg.ray.direction_tol:
            statuses.append('no_finite_axis_crossing')
            continue
        crossing[i] = p[2]-p[1]*d[2]/d[1]
        statuses.append('ok')
    return dict(pupil_q=q.copy(), axis_crossing_z_mm=crossing,
                lsa_mm=crossing-reference_focus_z, reference_focus_z_mm=float(reference_focus_z),
                paraxial_focus_z_mm=float(paraxial_focus_z),
                statuses=tuple(statuses), wavelength_nm=wavelength, field=field)


def differential_field_curvature(system, field, wavelength_nm=None, *, reference_z=None,
                                 image_z=None, derivative_step=None, convergence_rtol=1e-4,
                                 coupling_rtol=1e-5, settings=None):
    """Parabasal T/S focus from actually reaimed chief-neighbor real rays.

    Three step sizes test derivative convergence. Significant T/S coupling is
    reported unsupported instead of inventing diagonal focal values. The
    returned converging-positive local vergence is an independent convention.
    """
    cfg = settings or EngineSettings()
    wavelength = cfg.wavelength(wavelength_nm)
    zima = system.image_z if image_z is None else image_z
    z0 = reference_z if reference_z is not None else zima
    if z0 is None:
        raise MissingDataError('Field curvature needs a finite reference_z or image_z')
    if not np.isfinite(z0):
        raise InputError('reference_z must be finite')
    h0 = cfg.ray.derivative_step if derivative_step is None else float(derivative_step)
    if not np.isfinite(h0) or h0 <= 0 or h0 > 1 or convergence_rtol <= 0 or coupling_rtol <= 0:
        raise InputError('Step must be in (0,1]; convergence tolerances positive')
    t, s = tangential_sagittal_basis(field)
    basis = np.column_stack([t, s])
    values = []
    for h in (h0, h0/2, h0/4):
        jx, jb = np.empty((2, 2)), np.empty((2, 2))
        for axis in range(2):
            xs, slopes = [], []
            for sign in (-1, 1):
                tr = _trace_at(system, field, sign*h*basis[:, axis], wavelength, cfg)
                p, d = _last(tr)
                if abs(d[2]) <= cfg.ray.direction_tol:
                    raise ModelApplicabilityError('Parabasal ray grazing reference plane')
                xs.append(p[:2]+(z0-p[2])*d[:2]/d[2])
                slopes.append(d[:2]/d[2])
            jx[:, axis] = basis.T @ (xs[1]-xs[0])/(2*h)
            jb[:, axis] = basis.T @ (slopes[1]-slopes[0])/(2*h)
        values.append((jx, jb))
    jx, jb = values[-1]
    # J_X vanishes at a perfect local focus. Relative error divided by J_X
    # alone would reject the most nearly focused case. Use a geometric scale
    # (STOP radius transported by the local slope) as an absolute floor.
    pupil_scale = float(system.stop_semi_diameter)
    last_z = float(system.surfaces[-1].z)
    flight_scale = max(abs(z0-last_z), pupil_scale)
    scales = (max(np.linalg.norm(jx), pupil_scale, flight_scale*np.linalg.norm(jb)),
              max(np.linalg.norm(jb), pupil_scale/flight_scale))
    absolute_changes = [float(np.linalg.norm(a-b)) for a, b in zip(values[-1], values[-2])]
    changes = [difference/scale for difference, scale in zip(absolute_changes, scales)]
    if max(changes) > convergence_rtol:
        raise NumericalError(f'Parabasal derivatives did not converge: {changes}')
    offdiag = lambda a: np.linalg.norm(a-np.diag(np.diag(a)))
    coupled = any(offdiag(a) > coupling_rtol*scale for a, scale in zip((jx, jb), scales))
    out = dict(status='unsupported_coupling' if coupled else 'ok', jacobian_position_TS=jx,
               jacobian_slope_TS=jb, reference_z_mm=float(z0), image_z_mm=zima,
               derivative_relative_change=changes, derivative_absolute_change=absolute_changes,
               derivative_comparison_scales=scales, derivative_step=h0, basis_T=t, basis_S=s,
               wavelength_nm=wavelength, vergence_convention='converging_positive_independent')
    if coupled:
        return out
    delta = np.full(2, np.inf)
    stable = np.abs(np.diag(jb)) > cfg.ray.direction_tol
    delta[stable] = -np.diag(jx)[stable]/np.diag(jb)[stable]
    degenerate = (~stable) & (np.abs(np.diag(jx)) <= cfg.ray.intersection_tol_mm)
    delta[degenerate] = np.nan
    ni = system.index(system.surfaces[-1].medium_after, wavelength)
    with np.errstate(divide='ignore', invalid='ignore'):
        vergence = 1000*ni/delta
    out.update(focus_z_TS_mm=z0+delta, offset_TS_mm=delta,
               field_curvature_TS_mm=None if zima is None else z0+delta-zima,
               astigmatic_focus_separation_mm=float(delta[0]-delta[1]),
               vergence_TS_per_m=vergence, finite_focus_mask=stable & ~degenerate)
    if np.any(degenerate):
        out['status'] = 'degenerate_ray_family'
    elif not np.all(stable):
        out['status'] = 'no_finite_focus'
    return out


def zernike_terms(n_max: int):
    """Explicit (n,m,0/c/s) ordering, full-unit-disk RMS normalization."""
    if isinstance(n_max, bool) or not isinstance(n_max, (int, np.integer)) or n_max < 0:
        raise InputError('n_max must be a nonnegative integer')
    return tuple((n, m, kind) for n in range(n_max+1) for m in range(n+1)
                 if (n-m) % 2 == 0 for kind in (('0',) if m == 0 else ('c', 's')))


def zernike_basis(pupil_xy, terms):
    """Stable radial evaluation using Jacobi polynomials, equivalent to G9.1."""
    xy = np.asarray(pupil_xy, float)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy)):
        raise InputError('Zernike coordinates must be finite Nx2 normalized STOP points')
    rho = np.linalg.norm(xy, axis=1)
    if np.any(rho > 1+32*np.finfo(float).eps):
        raise InputError('Zernike STOP coordinates must lie inside the unit disk')
    rho = np.minimum(rho, 1.)
    phi = np.arctan2(xy[:, 1], xy[:, 0])
    columns = []
    for n, m, kind in terms:
        if n < 0 or m < 0 or m > n or (n-m) % 2 or kind not in (('0',) if m == 0 else ('c', 's')):
            raise InputError('Invalid Zernike tuple')
        radial = rho**m * eval_jacobi((n-m)//2, 0, m, 2*rho*rho-1)
        if m == 0:
            columns.append(np.sqrt(n+1)*radial)
        else:
            angular = np.cos(m*phi) if kind == 'c' else np.sin(m*phi)
            columns.append(np.sqrt(2*(n+1))*radial*angular)
    return np.column_stack(columns) if columns else np.empty((len(xy), 0))


def zernike_fit(pupil_xy, opd_mm, *, n_max=6, weights=None, mask=None,
                removed_terms=(), rcond=None, wavelength_nm=None, settings=None,
                reference_metadata=None):
    """Weighted SVD fit on declared normalized STOP chart; immutable diagnostics.

    OPD must ALREADY include the wave module's launch phase and reference
    geometry. Raw ray OPL differences are not a valid input substitute. Rank
    deficiency is explicit; minimum-norm diagnostic coefficients are retained
    but never called uniquely identified aberration coefficients.
    """
    if rcond is not None and (not np.isfinite(rcond) or rcond < 0):
        raise InputError('rcond must be None or a finite nonnegative SVD threshold')
    xy, wraw = np.asarray(pupil_xy, float), np.asarray(opd_mm, float)
    if xy.ndim != 2 or xy.shape[1] != 2 or wraw.shape != (len(xy),):
        raise InputError('OPD and pupil coordinate shape mismatch')
    w = np.ones(len(xy)) if weights is None else _weights(weights, len(xy))
    valid = np.ones(len(xy), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if valid.shape != (len(xy),):
        raise InputError('Mask shape mismatch')
    valid = valid & (w > 0)
    if not np.any(valid) or not np.all(np.isfinite(wraw[valid])):
        raise InputError('At least one finite positive-weight OPD sample is required')
    terms = zernike_terms(n_max)
    basis = zernike_basis(xy[valid], terms)
    sw = np.sqrt(w[valid])
    matrix = basis*sw[:, None]
    coeff, _, rank, singular = np.linalg.lstsq(matrix, wraw[valid]*sw, rcond=rcond)
    fit = basis @ coeff
    residual = wraw[valid]-fit
    remove = tuple(tuple(t) for t in removed_terms)
    if any(t not in terms for t in remove):
        raise InputError('Removed terms must occur in the fitted basis')
    selection = [terms.index(t) for t in remove]
    diagnostic = wraw[valid].copy()
    if selection:
        diagnostic -= basis[:, selection] @ coeff[selection]
    wn = w[valid]/w[valid].sum()
    threshold = (np.finfo(float).eps*max(matrix.shape) if rcond is None else rcond) * singular[0]
    condition = np.inf if rank < len(terms) else float(singular[0]/singular[-1])
    cfg = settings or EngineSettings()
    wavelength = cfg.wavelength(wavelength_nm)
    return dict(status='ok' if rank == len(terms) else 'fit_rank_deficient',
                terms=terms, coefficients_mm=coeff, coefficients_waves=coeff/(wavelength*1e-6),
                wavelength_nm=wavelength, rank=int(rank), singular_values=singular,
                rank_threshold=float(threshold), condition_number=condition,
                valid_mask=valid, raw_opd_mm=wraw.copy(), fitted_opd_mm=fit,
                residual_mm=residual, residual_rms_mm=float(np.sqrt(wn @ residual**2)),
                residual_pv_sampled_mm=float(np.ptp(residual)), diagnostic_opd_mm=diagnostic,
                diagnostic_rms_about_mean_mm=float(np.sqrt(wn @ (diagnostic-wn@diagnostic)**2)),
                removed_terms=remove, removal_policy='joint_full_basis_subtraction',
                reference_metadata=dict(reference_metadata or {}), chart='normalized_physical_STOP_XY',
                normalization='full_disk_uniform_area_RMS_1')


def focus_from_ray_lines(positions, slopes, *, pupil_weights=None, field_weights=None,
                         reference_z=0., bounds=None):
    """One closed-form geometric centroid-RMS focus for multiple field bundles.

    Each field's moments are normalized independently before field weighting.
    Input ray positions refer to the same ``reference_z`` plane; slopes are
    dx/dz,dy/dz. Membership, pupil illumination and sensor clipping must remain
    fixed. This does not optimize a per-field best-focus envelope.
    """
    aa = [np.asarray(a, float) for a in positions]
    bb = [np.asarray(b, float) for b in slopes]
    if not aa or len(aa) != len(bb) or not np.isfinite(reference_z):
        raise InputError('Provide matching nonempty field ray-line groups and finite reference_z')
    if pupil_weights is None:
        pupil_weights = [np.ones(len(a)) for a in aa]
    if len(pupil_weights) != len(aa):
        raise InputError('Pupil weight group count differs from fields')
    fw = _weights(np.ones(len(aa)) if field_weights is None else field_weights, len(aa))
    fw = fw/fw.sum()
    moments = []
    for a, b, weights in zip(aa, bb, pupil_weights):
        if a.ndim != 2 or a.shape[1] != 2 or b.shape != a.shape or not np.all(np.isfinite([a, b])):
            raise InputError('Ray position and slope groups must be finite matching Nx2 arrays')
        w = _weights(weights, len(a)); w = w/w.sum()
        da, db = a-w@a, b-w@b
        moments.append([np.sum(w*np.sum(da*db, axis=1)), np.sum(w*np.sum(db*db, axis=1)),
                        np.sum(w*np.sum(da*da, axis=1))])
    moments = np.array(moments)
    total = fw @ moments
    if total[1] == 0:
        return dict(status='no_unique_finite_focus', common_image_z_mm=None,
                    offset_mm=None, field_moments=moments, field_weights=fw)
    unconstrained = -total[0]/total[1]
    delta = unconstrained
    if bounds is not None:
        if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or bounds[0] > bounds[1]:
            raise InputError('Bounds must be increasing finite absolute image_z limits')
        delta = float(np.clip(reference_z+delta, *bounds)-reference_z)
    with np.errstate(divide='ignore', invalid='ignore'):
        local = np.where(moments[:, 1] > 0, -moments[:, 0]/moments[:, 1], np.nan)
    squared = np.maximum(0., moments[:, 2]+2*delta*moments[:, 0]+delta**2*moments[:, 1])
    return dict(status='ok', common_image_z_mm=float(reference_z+delta), offset_mm=float(delta),
                unconstrained_image_z_mm=float(reference_z+unconstrained), reference_z_mm=float(reference_z),
                local_best_focus_diagnostic_mm=reference_z+local,
                per_field_rms_mm=np.sqrt(squared), merit_rms_mm=float(np.sqrt(fw @ squared)),
                field_moments=moments, field_weights=fw, bounds=bounds,
                criterion='field_weighted_centroid_geometric_RMS_squared',
                illumination_policy='frozen_ray_membership_and_pupil_weights')


def common_geometric_refocus(system, bundles, *, image_z=None, field_weights=None,
                             bounds=None, settings=None):
    """Common geometric plane optimization from traced bundles, then all-field spots."""
    cfg = settings or EngineSettings()
    z = _image_z(system, image_z)
    bundles = tuple(bundles)
    positions, slopes, weights = [], [], []
    for bundle in bundles:
        require_bundle_match(system, bundle, cfg, allow_unverified=True)
        mask, w = _valid_bundle(bundle)
        # Refocus minimizes the quadrature RMS, not the geometric envelope.
        ids = np.flatnonzero(mask & (w > 0))
        p, b = [], []
        for i in ids:
            tr = bundle.traces[i]
            p.append(ray_plane_intercept(tr, z, direction_tolerance=cfg.ray.direction_tol)[0])
            _, d = _last(tr)
            b.append(d[:2]/d[2])
        positions.append(p); slopes.append(b); weights.append(w[ids])
    result = focus_from_ray_lines(positions, slopes, pupil_weights=weights,
                                  field_weights=field_weights, reference_z=z, bounds=bounds)
    if result['status'] == 'ok':
        result['fixed_plane_spots'] = [spot_diagram(system, b, z, settings=cfg) for b in bundles]
        result['common_focus_spots'] = [spot_diagram(system, b, result['common_image_z_mm'], settings=cfg) for b in bundles]
    return result
