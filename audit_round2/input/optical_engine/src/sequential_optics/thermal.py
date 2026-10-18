"""Uniform free expansion with explicit idealized axial placement constraints.

This module builds geometry, not a thermo-mechanical FEM solution. Every state
is derived from the reference prescription, never from the preceding hot state.
Index and absorption data are requested lazily by the analyses that need them.
"""
from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.optimize import minimize_scalar

from .errors import InputError, GeometryError, ModelApplicabilityError


@dataclass(frozen=True)
class LensGroup:
    """An uncemented homogeneous lens bounded by consecutive surface vertices.

    anchor_fraction=0/.5/1 fixes front/vertex-midpoint/back, respectively.
    anchor_shift_per_k_mm is an explicitly prescribed idealized mount motion;
    no such displacement is inferred from the glass CTE.
    """
    front_index: int
    back_index: int
    material_id: str | None = None
    anchor_fraction: float = .5
    anchor_shift_per_k_mm: float = 0.

    def __post_init__(self):
        if self.back_index != self.front_index + 1 or self.front_index < 0:
            raise InputError('A LensGroup requires consecutive front/back surfaces')
        if not math.isfinite(self.anchor_fraction) or not 0 <= self.anchor_fraction <= 1:
            raise InputError('Anchor fraction must be in [0,1]')
        if not math.isfinite(self.anchor_shift_per_k_mm):
            raise InputError('Anchor thermal displacement must be finite')

    def to_dict(self):
        from dataclasses import asdict
        return asdict(self)

    @classmethod
    def from_dict(cls, value):
        return cls(**value)


def infer_lens_groups(system, anchor_fraction=.5):
    """Infer only unambiguous isolated two-surface lenses; reject cemented data."""
    groups = []
    surfaces = system.surfaces
    i = 0
    while i < len(surfaces):
        s = surfaces[i]
        if s.medium_after != system.ambient_material:
            if i + 1 >= len(surfaces) or surfaces[i + 1].medium_after != system.ambient_material:
                raise ModelApplicabilityError('Thermal auto grouping supports only isolated homogeneous two-surface lenses')
            groups.append(LensGroup(i, i + 1, s.medium_after, anchor_fraction))
            i += 2
        else:
            i += 1
    return tuple(groups)


def _validated_groups(system, lens_groups, anchor_fraction):
    groups = infer_lens_groups(system, anchor_fraction) if lens_groups is None else tuple(lens_groups)
    occupied = set()
    previous = -1
    for g in groups:
        if g.back_index >= len(system.surfaces) or g.front_index <= previous:
            raise InputError('Lens groups must be ordered, disjoint, and inside the surface list')
        f, b = system.surfaces[g.front_index], system.surfaces[g.back_index]
        before = system.ambient_material if g.front_index == 0 else system.surfaces[g.front_index - 1].medium_after
        if f.medium_after == system.ambient_material or before != system.ambient_material or b.medium_after != system.ambient_material:
            raise ModelApplicabilityError('Each thermal lens must be isolated in the ambient medium')
        if g.material_id is not None and g.material_id != f.medium_after:
            raise InputError('LensGroup material must match its internal propagation material')
        if b.z <= f.z: raise GeometryError('Lens center thickness must be positive')
        occupied.update((g.front_index, g.back_index)); previous = g.back_index
    for i, s in enumerate(system.surfaces):
        if s.medium_after != system.ambient_material and i not in occupied:
            raise InputError('Every solid-material interval must belong to a thermal LensGroup')
    return groups


def scale_surface(surface, scale, z=None):
    """z_T(r)=s z_0(r/s); handles c=0 with nonzero A terms correctly."""
    if not math.isfinite(scale) or scale <= 0: raise InputError('Positive finite expansion scale required')
    return replace(surface, z=surface.z if z is None else z,
                   curvature=surface.curvature / scale,
                   coefficients=tuple(a * scale ** (1 - 2 * (i + 1)) for i, a in enumerate(surface.coefficients)),
                   clear_radius=surface.clear_radius * scale,
                   inactive_height_mm=None if surface.inactive_height_mm is None else surface.inactive_height_mm*scale,
                   mechanical_radius=None if surface.mechanical_radius is None else surface.mechanical_radius * scale)


def validate_radial_separation(system, samples=257, tolerance_mm=1e-10):
    """Numerical full-aperture surface separation check (not a proof certificate).

    Samples the shared optical clear aperture and refines every sampled local
    minimum of the adjacent gap. Mechanical rims have no specified shape, so
    this cannot certify their contact or mount clearances. Strongly oscillatory
    near-singular prescriptions should independently increase samples.
    """
    from .geometry import sag
    if samples < 5: raise InputError('Separation validation needs at least 5 radial samples')
    for i, (a, b) in enumerate(zip(system.surfaces, system.surfaces[1:])):
        radius = min(a.clear_radius, b.clear_radius)
        rs = np.linspace(0., radius, int(samples))
        def gap(r): return float(b.z + sag(b, r, 0.) - a.z - sag(a, r, 0.))
        gs = np.array([gap(r) for r in rs])
        if not np.all(np.isfinite(gs)):
            raise GeometryError(f'Nonfinite sag during separation validation at surfaces {i},{i+1}')
        min_gap = float(gs.min())
        for k in range(1, len(rs) - 1):
            if gs[k] <= gs[k - 1] and gs[k] <= gs[k + 1]:
                opt = minimize_scalar(gap, bounds=(rs[k-1], rs[k+1]), method='bounded',
                                      options={'xatol': tolerance_mm})
                if opt.success: min_gap = min(min_gap, float(opt.fun))
        if min_gap < -tolerance_mm:
            raise GeometryError(f'Surfaces {i},{i+1} overlap inside clear aperture (gap {min_gap:g} mm)')
        # Coincident dummy surfaces are allowed, but a solid lens cannot have
        # zero/negative material thickness even when its vertices are ordered.
        if a.medium_after != system.ambient_material and min_gap <= tolerance_mm:
            raise GeometryError(f'Lens at surfaces {i},{i+1} has zero/nonpositive optical thickness')


def thermal_system(system, temperature_c, *, lens_groups=None,
                   anchor_fraction=.5, layout='fixed_anchors',
                   stop_policy='follow_surface', image_policy='fixed',
                   object_policy='fixed', validate_geometry=True,
                   validation_samples=257, tolerance_mm=1e-10):
    """Create an independent T state from a baseline OpticalSystem.

    ``layout='fixed_anchors'`` holds each lens's declared anchor (plus supplied
    shift). ``fixed_gaps`` holds all lens-to-lens vertex air gaps and only uses
    the first lens anchor; non-first anchor motions are rejected as conflicting.
    Non-lens dummy surfaces retain their global vertices.

    STOP follows its assigned surface's position. ``follow_surface`` scales
    its radius with that lens, while ``fixed_radius`` holds a separate physical
    aperture radius. An independent dummy STOP is fixed under either policy.
    A held radius exceeding the expanded lens clear aperture is rejected, not
    silently resized. EPD/F-number must be resolved at baseline beforehand.

    Object/image policies are fixed global coordinate or first/last surface
    offset. No refocus is performed here. Source dimensions remain unchanged.
    This geometry operation does not request n(T), loss, or spectral data;
    index() during subsequent analysis validates the needed thermo-optic data.
    At T==T0, no CTE is required. Nonbaseline input geometry is rejected rather
    than cumulatively re-expanded. Call this function afresh on the baseline.
    """
    t, t0 = float(temperature_c), float(system.reference_temperature_c)
    if not math.isfinite(t): raise InputError('Temperature must be finite')
    if system.temperature_c is not None and system.temperature_c != t0:
        raise InputError('thermal_system requires original reference geometry; do not re-expand a nonreference state')
    if layout not in ('fixed_anchors', 'fixed_gaps'):
        raise InputError('Unknown lens thermal layout')
    if stop_policy not in ('follow_surface', 'fixed_radius'):
        raise InputError('Unknown STOP thermal policy')
    if image_policy not in ('fixed', 'last_surface_offset') or object_policy not in ('fixed', 'first_surface_offset'):
        raise InputError('Unknown object/image thermal position policy')
    groups = _validated_groups(system, lens_groups, anchor_fraction)
    if layout == 'fixed_gaps' and any(g.anchor_shift_per_k_mm != 0 for g in groups[1:]):
        raise InputError('Fixed gaps and independent downstream anchor motion overconstrain layout')
    if t == t0:
        # Preserve exactly the supplied coordinates (including large-offset
        # prescriptions), rather than subtracting/readding rounded anchors.
        result = replace(system, temperature_c=t,
                         stop_radius=system.stop_semi_diameter if stop_policy == 'fixed_radius' else system.stop_radius)
        if validate_geometry: validate_radial_separation(result, validation_samples, tolerance_mm)
        return result
    surfaces = list(system.surfaces)
    scales = {}
    previous_group = None
    for g in groups:
        f0, b0 = system.surfaces[g.front_index], system.surfaces[g.back_index]
        material = system.materials[f0.medium_after]
        scale = material.expansion_scale(t, t0)
        ct = b0.z - f0.z
        if layout == 'fixed_gaps' and previous_group is not None:
            old_gap = f0.z - system.surfaces[previous_group.back_index].z
            front_z = surfaces[previous_group.back_index].z + old_gap
        else:
            anchor = f0.z + g.anchor_fraction * ct + g.anchor_shift_per_k_mm * (t - t0)
            front_z = anchor - g.anchor_fraction * scale * ct
        back_z = front_z + scale * ct
        surfaces[g.front_index] = scale_surface(f0, scale, front_z)
        surfaces[g.back_index] = scale_surface(b0, scale, back_z)
        scales[g.front_index] = scale; scales[g.back_index] = scale
        previous_group = g
    stop_radius = system.stop_radius
    if stop_policy == 'fixed_radius': stop_radius = system.stop_semi_diameter
    elif stop_radius is not None: stop_radius *= scales.get(system.stop_index, 1.)
    image_z, object_z = system.image_z, system.object_z
    if image_z is not None and image_policy == 'last_surface_offset':
        image_z += surfaces[-1].z - system.surfaces[-1].z
    if object_z is not None and object_policy == 'first_surface_offset':
        object_z += surfaces[0].z - system.surfaces[0].z
    result = replace(system, surfaces=tuple(surfaces), temperature_c=t,
                     stop_radius=stop_radius, image_z=image_z, object_z=object_z)
    if validate_geometry: validate_radial_separation(result, validation_samples, tolerance_mm)
    return result
