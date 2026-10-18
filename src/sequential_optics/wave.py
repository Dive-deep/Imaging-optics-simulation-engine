"""Independent scalar diffraction profiles, with explicit physical measures.

Arrays use [y, x]. Lengths are mm; ``wavelength_mm`` is the vacuum wavelength.
No pupil fit removes tilt or defocus. This is not a Zemax compatibility adapter.
Equations and applicability: docs/calculation_spec_ko.md, chapter 4 (W-1..W-14).
"""
from dataclasses import dataclass, replace
import hashlib
import numpy as np
from scipy.spatial import Delaunay
from .errors import InputError, ModelApplicabilityError, NumericalError
from .settings import EngineSettings
from .provenance import (require_bundle_match, trace_configuration_fingerprint,
                         bundle_content_fingerprint)

PROFILE_ID = 'independent_scalar_fft_v1'
VACUUM_IMPEDANCE_OHM = 376.730313668


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise InputError(f'{name} must be finite and positive')
    return value


def _spacing(spacing):
    if np.isscalar(spacing):
        return (_positive(spacing, 'spacing'),) * 2
    if len(spacing) != 2:
        raise InputError('spacing must be dx or (dx,dy)')
    return _positive(spacing[0], 'dx'), _positive(spacing[1], 'dy')


def _grid_field(field):
    a = np.asarray(field, dtype=complex)
    if a.ndim != 2 or min(a.shape) < 4 or any(n % 2 for n in a.shape):
        raise InputError('The scalar profile requires an even 2-D grid, at least 4 per axis')
    if not np.all(np.isfinite(a)):
        raise InputError('Pupil field must be finite')
    if not np.any(np.abs(a) > 0):
        raise InputError('Zero field has no normalized PSF')
    return a


def _fft(field, spacing, padding):
    a = _grid_field(field)
    if int(padding) != padding or padding < 1:
        raise InputError('padding must be a positive integer')
    dx, dy = _spacing(spacing)
    ny, nx = a.shape
    py, px = ny * int(padding), nx * int(padding)
    padded = np.zeros((py, px), dtype=complex)
    sy, sx = (py - ny) // 2, (px - nx) // 2
    padded[sy:sy + ny, sx:sx + nx] = a
    f = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(padded))) * dx * dy
    fx = np.fft.fftshift(np.fft.fftfreq(px, dx))
    fy = np.fft.fftshift(np.fft.fftfreq(py, dy))
    return f, fx, fy, padded


def _edge_fraction(power):
    edge = np.zeros(power.shape, dtype=bool)
    edge[[0, -1], :] = True
    edge[:, [0, -1]] = True
    return float(power[edge].sum() / power.sum())


def otf_from_psf(psf, dx, dy=None):
    """Complex OTF using the forward-negative Fourier sign, with OTF(0)=1.

    ``psf`` is a nonnegative density or intensity on an origin-centered grid.
    No centroid registration and no MTF monotonicization are performed.
    """
    a = np.asarray(psf, dtype=float)
    if a.ndim != 2 or min(a.shape) < 2 or np.any(~np.isfinite(a)) or np.any(a < 0):
        raise InputError('PSF must be a finite nonnegative 2-D array')
    dx = _positive(dx, 'dx'); dy = _positive(dx if dy is None else dy, 'dy')
    total = float(a.sum() * dx * dy)
    if total <= 0:
        raise InputError('Zero PSF has no normalized OTF')
    otf = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(a))) * dx * dy / total
    return dict(otf=otf, mtf=np.abs(otf),
                frequency_x=np.fft.fftshift(np.fft.fftfreq(a.shape[1], dx)),
                frequency_y=np.fft.fftshift(np.fft.fftfreq(a.shape[0], dy)),
                normalization=total)


def sample_axis_otf(otf_result, frequencies, axis='x', support_cutoff=None):
    """Interpolate complex OTF, never its magnitude. Reject unresolved bandwidth.

    A proven pupil autocorrelation support bound can explicitly justify zero
    beyond ``support_cutoff``; it is not inferred from the measured MTF.
    """
    f = np.asarray(frequencies, dtype=float)
    if not np.all(np.isfinite(f)) or axis not in ('x', 'y'):
        raise InputError('Finite frequencies and x/y axis are required')
    data = otf_result['otf']
    if axis == 'x':
        grid = otf_result['frequency_x']; line = data[data.shape[0] // 2, :]
    else:
        grid = otf_result['frequency_y']; line = data[:, data.shape[1] // 2]
    outside = (f < grid[0]) | (f > grid[-1])
    known_zero = np.zeros(f.shape, dtype=bool)
    if support_cutoff is not None:
        known_zero = np.abs(f) >= float(support_cutoff)
    if np.any(outside & ~known_zero):
        raise ModelApplicabilityError('Requested OTF bandwidth exceeds PSF Nyquist; increase padding')
    result = np.interp(f, grid, line.real) + 1j * np.interp(f, grid, line.imag)
    return np.where(known_zero, 0j, result)


def _mtf_axis_metadata(basis):
    """Declare frequency directions separately from perpendicular target lines.

    The independent FFT profile uses a fixed chief-normal frame: S=e1 and
    T=e2, not a field-radial rotation. This follows the fixed object-X/Y naming
    intent of OpticStudio FFT MTF, without claiming identical vendor mapping.
    For afocal curves these are angular-chart axes about the chief direction.
    """
    s, t = np.asarray(basis, float)[:2]
    return dict(axis_convention_id='fixed_chief_normal_xy',
                mtf_frequency_axes_global={'s':s.copy(), 't':t.copy()},
                mtf_target_line_axes_global={'s':t.copy(), 't':s.copy()},
                mtf_axis_meaning='frequency direction; target lines are perpendicular within the chief-normal plane',
                field_radial_axes=False)


def fft_psf(pupil, spacing, wavelength_mm, *, index=1., padding=4):
    """Focal Fourier profile from uniform (chi,eta) complex pupil.

    ``spacing=(dchi,deta)`` is dimensionless. The returned output is the
    chief-normal virtual plane in mm. The same-pupil zero-phase diffraction
    reference can be computed with ``fft_psf(abs(pupil),...)``.
    """
    wavelength = _positive(wavelength_mm, 'vacuum wavelength')
    n = _positive(index, 'index'); lm = wavelength / n
    a = _grid_field(pupil); dc, de = _spacing(spacing)
    amplitude, fx, fy, padded = _fft(a, (dc, de), padding)
    x = lm * fx; y = lm * fy
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0])
    intensity = np.abs(amplitude) ** 2
    raw = float(intensity.sum() * dx * dy)
    if not np.isfinite(raw) or raw <= 0:
        raise NumericalError('Invalid Fourier normalization')
    psf = intensity / raw
    otf = otf_from_psf(psf, dx, dy)
    yy, xx = np.where(np.abs(a) > 0)
    widths = (int(xx.max() - xx.min() + 1), int(yy.max() - yy.min() + 1))
    linear = padded.shape[1] >= 2 * widths[0] - 1 and padded.shape[0] >= 2 * widths[1] - 1
    return dict(psf=psf, coherent_field=amplitude, x_mm=x, y_mm=y,
                dx_mm=dx, dy_mm=dy, raw_integral=raw,
                pupil_squared_integral=float(np.sum(np.abs(a) ** 2) * dc * de),
                edge_energy_fraction=_edge_fraction(intensity),
                physical_period_mm=(lm / dc, lm / de),
                input_grid_shape=a.shape, padded_shape=padded.shape,
                linear_autocorrelation_supported=bool(linear),
                wavelength_vacuum_mm=wavelength, refractive_index=n,
                frequency_units='cycles/mm', output_plane='chief_normal_virtual', **otf)


def circular_mtf(frequency, cutoff):
    """Analytic incoherent MTF of a uniform unobscured circular Fourier pupil."""
    fc = _positive(cutoff, 'cutoff')
    q = np.abs(np.asarray(frequency, dtype=float)) / fc
    if not np.all(np.isfinite(q)):
        raise InputError('Frequencies must be finite')
    z = np.clip(q, 0., 1.)
    return np.where(q < 1., 2 / np.pi * (np.arccos(z) - z * np.sqrt(1 - z*z)), 0.)


def afocal_fft(boundary_field, spacing_mm, wavelength_mm, *, index=1., padding=4):
    """Scalar far field in exact direction cosines and unit-solid-angle density.

    h_solid_angle is sr^-1, while cell_power includes dl*dm/gamma. Frequencies
    from this DFT are NOT angular cycles/rad. Use ``angular_otf`` for that.
    Only the represented propagating angular window is normalized.
    """
    wavelength = _positive(wavelength_mm, 'vacuum wavelength')
    n = _positive(index, 'index'); lm = wavelength / n
    du, dv = _spacing(spacing_mm)
    f, nu, nv, padded = _fft(boundary_field, (du,dv), padding)
    l = lm * nu; m = lm * nv
    ll, mm = np.meshgrid(l,m)
    prop = ll*ll + mm*mm < 1.
    gamma = np.sqrt(np.maximum(0., 1 - ll*ll - mm*mm))
    dl = float(l[1]-l[0]); dm = float(m[1]-m[0])
    squared = np.abs(f)**2
    mass = np.where(prop, gamma*squared*dl*dm, 0.)
    raw = float(mass.sum())
    if raw <= 0 or not np.isfinite(raw):
        raise NumericalError('No represented propagating wave power')
    cell_power = mass / raw
    h = np.where(prop, gamma**2*squared/raw, 0.)
    alpha = np.arctan2(ll,gamma); beta = np.arctan2(mm,gamma)
    j = angular_chart_jacobian(alpha,beta)
    # All lengths below are mm: Fourier field has amplitude*mm^2 and
    # the final 1e-6 converts area mm^2 to m^2 for physical electric fields.
    flux_scale = n/(2*VACUUM_IMPEDANCE_OHM*lm**2)*1e-6
    return dict(fourier_field=f,l=l,m=m,gamma=gamma,propagating_mask=prop,
                h_solid_angle=h,cell_power=cell_power,alpha_rad=alpha,beta_rad=beta,
                angular_chart_jacobian=j,h_angle_density=h*j,
                dl=dl,dm=dm,raw_angular_integral=raw,
                electric_field_wave_power_w=raw*flux_scale,
                sampled_nonpropagating_spectral_fraction=float(squared[~prop].sum()/squared.sum()),
                edge_energy_fraction=_edge_fraction(cell_power),
                direction_cosine_period=(lm/du,lm/dv),
                padded_shape=padded.shape,normalization_domain='represented_propagating_window',
                frequency_units='cycles/rad (angular_otf only)',
                model='scalar_RS_far_field_obliquity',wavelength_vacuum_mm=wavelength,index=n)


def angular_chart_jacobian(alpha,beta):
    """dOmega/dalpha/dbeta for alpha=atan(dx/dz), beta=atan(dy/dz)."""
    a,b=np.broadcast_arrays(np.asarray(alpha,float),np.asarray(beta,float))
    # sin/cos form avoids overflow of tan near hemisphere edges.
    ca=np.cos(a);cb=np.cos(b);sa=np.sin(a);sb=np.sin(b)
    denominator=(ca*ca*cb*cb+sa*sa*cb*cb+sb*sb*ca*ca)**1.5
    return np.divide(np.abs(ca*cb),denominator,out=np.zeros_like(denominator),where=denominator>0)


def _angular_cell_width(result, axis):
    """Conservative exact-angle span of each direction-cosine cell.

    Include corners and the minimum-absolute cross coordinate; an angle's
    extremum may occur at m=0 (or l=0), not just at a rectangle corner.
    Cells touching the propagating disk edge use gamma=0 as a conservative
    angular bound. This is not a linear direction-cosine frequency cutoff.
    """
    l, m = np.meshgrid(result['l'], result['m'])
    a, b = (l, m) if axis == 'x' else (m, l)
    da, db = (result['dl'], result['dm']) if axis == 'x' else (result['dm'], result['dl'])
    low = np.full(a.shape, np.inf); high = np.full(a.shape, -np.inf)
    cross = (b-db/2, b+db/2, np.maximum(b-db/2, np.minimum(0., b+db/2)))
    for aa in (a-da/2, a+da/2):
        for bb in cross:
            gamma = np.sqrt(np.maximum(0., 1-aa*aa-bb*bb))
            angle = np.arctan2(aa, gamma)
            low = np.minimum(low, angle); high = np.maximum(high, angle)
    return high-low


def angular_otf(result, frequencies_cycles_rad, *, axis='x', settings=None,
                return_diagnostics=False):
    """Power-weighted midpoint quadrature in exact angular coordinates.

    A frequency is unresolved when cells exceeding the configured phase-span
    limit contain more than the configured fractional power tolerance. The
    test is necessary, not a quadrature error/convergence guarantee. Low-level
    calls reject unresolved requests; diagnostic calls return NaN at those
    frequencies so high-level results cannot publish aliased numbers as MTF.
    No cells are removed and no frequency is silently replaced by zero.
    """
    if axis not in ('x','y'):
        raise InputError('axis must be x or y')
    cfg = (settings or EngineSettings()).sampling
    f = np.asarray(frequencies_cycles_rad, float)
    if not np.all(np.isfinite(f)):
        raise InputError('Angular frequencies must be finite')
    angle = np.asarray(result['alpha_rad' if axis=='x' else 'beta_rad'], float)
    weights = np.asarray(result['cell_power'], float)
    if (angle.shape != weights.shape or weights.ndim != 2 or
            np.any(~np.isfinite(weights)) or np.any(weights < 0) or
            np.any(~np.isfinite(angle)) or not np.isclose(weights.sum(), 1., rtol=1e-10, atol=1e-12)):
        raise InputError('Angular quadrature needs finite coordinates and normalized nonnegative cell power')
    width = _angular_cell_width(result, axis)
    use = weights > 0
    angle = angle[use]; weights = weights[use]; width = width[use]
    mass = np.array([weights[2*np.pi*abs(v)*width > cfg.angular_otf_phase_step_limit_rad].sum()
                     for v in f.ravel()]).reshape(f.shape)
    resolved = mass <= cfg.angular_otf_unresolved_power_tolerance
    if np.any(~resolved) and not return_diagnostics:
        raise ModelApplicabilityError('Angular OTF frequency quadrature is underresolved; increase FFT padding and verify convergence')
    values = np.full(f.shape, np.nan+1j*np.nan, complex)
    for i, v in enumerate(f.ravel()):
        if resolved.ravel()[i]:
            values.ravel()[i] = np.sum(weights*np.exp(-2j*np.pi*v*angle))
    if not return_diagnostics:
        return values
    return dict(otf=values, resolved_mask=resolved,
                unresolved_cell_power_fraction=mass,
                max_cell_angle_span_rad=float(width.max()),
                phase_step_limit_rad=cfg.angular_otf_phase_step_limit_rad,
                unresolved_power_tolerance=cfg.angular_otf_unresolved_power_tolerance,
                necessary_check_passed=bool(np.all(resolved)),
                convergence_verified=False,
                criterion='exact angular cell phase span; not an MTF error bound')


def _angular_mtf_axes(result, reference, frequencies, settings):
    """Apply the same output-quadrature gate to actual and reference MTF."""
    diagnostics = {}
    for axis, label in [('x', 's'), ('y', 't')]:
        actual = angular_otf(result, frequencies, axis=axis, settings=settings, return_diagnostics=True)
        ideal = angular_otf(reference, frequencies, axis=axis, settings=settings, return_diagnostics=True)
        result['otf_'+label] = actual['otf']
        result['mtf_'+label] = np.abs(actual['otf'])
        result['diffraction_mtf_'+label] = np.abs(ideal['otf'])
        diagnostics[label] = {k:v for k,v in actual.items() if k != 'otf'}
        diagnostics['diffraction_'+label] = {k:v for k,v in ideal.items() if k != 'otf'}
    diagnostics['necessary_check_passed'] = all(d['necessary_check_passed'] for d in diagnostics.values())
    result['angular_frequency_sampling'] = diagnostics


def direct_rs(boundary_field, u_mm, v_mm, targets_xyz_mm, wavelength_mm, *, index=1., area_weights_mm2=None):
    """Direct Rayleigh-Sommerfeld I scalar half-space quadrature (W-2).

    Targets are (x,y,L) relative to the boundary plane. With 1-D u/v axes,
    regular-grid midpoint cell areas are inferred. With arbitrary coordinate
    arrays, explicit positive ``area_weights_mm2`` are required. No FFT or
    near-focus expansion is used, but boundary-field and quadrature errors remain.
    """
    field=np.asarray(boundary_field,complex)
    u=np.asarray(u_mm,float);v=np.asarray(v_mm,float)
    if u.ndim==v.ndim==1 and field.shape==(len(v),len(u)):
        if len(u)<2 or len(v)<2:
            raise InputError('Need two or more samples per coordinate')
        uu,vv=np.meshgrid(u,v)
        if area_weights_mm2 is None:
            du=np.diff(u);dv=np.diff(v)
            if not np.allclose(du,du[0]) or not np.allclose(dv,dv[0]) or du[0]<=0 or dv[0]<=0:
                raise InputError('Nonuniform coordinates require explicit quadrature weights')
            area_weights_mm2=du[0]*dv[0]
    else:
        uu,vv=np.broadcast_arrays(u,v)
        if uu.shape!=field.shape or area_weights_mm2 is None:
            raise InputError('Coordinate arrays require matching field and explicit cell areas')
    weights=np.broadcast_to(np.asarray(area_weights_mm2,float),field.shape)
    if (np.any(weights<0) or not np.all(np.isfinite(weights)) or not np.all(np.isfinite(field))
            or not np.all(np.isfinite(uu)) or not np.all(np.isfinite(vv))):
        raise InputError('Field and nonnegative cell areas must be finite')
    target=np.asarray(targets_xyz_mm,float)
    if target.ndim==1: target=target[None,:]
    if target.ndim!=2 or target.shape[1]!=3 or np.any(target[:,2]<=0) or not np.all(np.isfinite(target)):
        raise InputError('Targets must be finite downstream (x,y,L) points')
    ki=2*np.pi*_positive(index,'index')/_positive(wavelength_mm,'wavelength')
    nonzero=np.abs(field)>0
    uu=uu[nonzero];vv=vv[nonzero];fw=field[nonzero]*weights[nonzero]
    out=[]
    for x,y,L in target:
        r=np.sqrt((x-uu)**2+(y-vv)**2+L*L)
        kernel=L/(2*np.pi*r**3)*(1-1j*ki*r)*np.exp(1j*ki*r)
        out.append(np.sum(fw*kernel))
    return np.array(out)


def focal_direct_reference(pupil, chi, eta, distance_mm, wavelength_mm, targets_xy_mm, *, index=1.):
    """Compare an existing focal Fourier pupil with exact scalar RS propagation.

    Inverts B=A*L/gamma^2 and the OPD reference phase, then applies the exact
    physical-area Jacobian. The same sampled boundary is shared with the FFT;
    this is independent propagation, not independent ray/OPD validation.
    """
    c=np.asarray(chi,float);e=np.asarray(eta,float);p=np.asarray(pupil,complex)
    if c.ndim!=1 or e.ndim!=1 or p.shape!=(len(e),len(c)) or len(c)<2 or len(e)<2:
        raise InputError('Pupil must match 1-D chi/eta coordinates')
    dc=np.diff(c);de=np.diff(e)
    if not np.allclose(dc,dc[0]) or not np.allclose(de,de[0]) or dc[0]<=0 or de[0]<=0:
        raise InputError('Uniform ascending chi/eta grid required')
    L=_positive(distance_mm,'distance');lam=_positive(wavelength_mm,'wavelength');n=_positive(index,'index')
    cc,ee=np.meshgrid(c,e);g2=1-cc*cc-ee*ee
    valid=g2>0
    if np.any((np.abs(p)>0)&~valid):
        raise ModelApplicabilityError('Nonzero pupil outside propagating reference chart')
    gamma=np.sqrt(np.where(valid,g2,1.))
    u=L*cc/gamma;v=L*ee/gamma;r=L/gamma
    boundary=p*gamma**2/L*np.exp(-2j*np.pi*n/lam*((u*u+v*v)/(r+L)))
    area=L*L/gamma**4*dc[0]*de[0]
    xy=np.asarray(targets_xy_mm,float)
    if xy.ndim==1: xy=xy[None,:]
    if xy.ndim!=2 or xy.shape[1]!=2: raise InputError('Targets must be Nx2')
    return direct_rs(boundary,u,v,np.column_stack([xy,np.full(len(xy),L)]),lam,index=n,area_weights_mm2=area)


@dataclass(frozen=True)
class BoundaryState:
    """Reusable, content-sealed chief-normal ray boundary snapshot.

    Reuse only for moving IMA with unchanged rays, field, wavelength and geometry.
    Base focal lens metric is B=1 at ``reference_distance_mm``; subsequent focus
    evaluations retain the same physical A(u,v)=L0/(L0^2+u^2+v^2).
    """
    origin: np.ndarray
    basis: np.ndarray
    uv_mm: np.ndarray
    opl_mm: np.ndarray
    chief_opl_mm: float
    valid_ray_mask: np.ndarray
    triangles: np.ndarray
    pupil_xy: np.ndarray
    sampled_jacobians_mm2: np.ndarray
    wavelength_nm: float
    output_index: float
    reference_distance_mm: float | None
    chief_position: np.ndarray
    chief_direction: np.ndarray
    fingerprint: str
    plane_bound: float
    clearance_mm: float
    source_amplitudes: np.ndarray | None = None
    trace_provenance_status: str = 'unverified_external_bundle'
    content_fingerprint: str = ''

    def __post_init__(self):
        # Bytes-backed views cannot have WRITEABLE re-enabled by setflags.
        # The producer seals only after construction; dataclasses.replace must
        # preserve the old seal rather than approving its changed values.
        names=('origin','basis','uv_mm','opl_mm','valid_ray_mask','triangles',
               'pupil_xy','sampled_jacobians_mm2','chief_position','chief_direction',
               'source_amplitudes')
        for name in names:
            value=getattr(self,name)
            if value is None: continue
            array=np.ascontiguousarray(value)
            fixed=np.frombuffer(array.tobytes(),dtype=array.dtype).reshape(array.shape)
            object.__setattr__(self,name,fixed)


def _boundary_content_fingerprint(state):
    """Detect accidentally modified cached geometry/amplitude, not authentication."""
    digest=hashlib.sha256()
    for name in ('origin','basis','uv_mm','opl_mm','valid_ray_mask','triangles',
                 'pupil_xy','sampled_jacobians_mm2','chief_position','chief_direction',
                 'source_amplitudes'):
        value=getattr(state,name)
        digest.update(name.encode())
        if value is None: digest.update(b'None');continue
        array=np.asarray(value)
        digest.update(repr((array.dtype.str,array.shape)).encode())
        digest.update(array.tobytes())
    scalars=(state.chief_opl_mm,state.wavelength_nm,state.output_index,
             state.reference_distance_mm,state.fingerprint,state.plane_bound,
             state.clearance_mm,state.trace_provenance_status)
    digest.update(repr(scalars).encode())
    return digest.hexdigest()


def _seal_boundary_state(state):
    return replace(state,content_fingerprint=_boundary_content_fingerprint(state))


def _validate_boundary_snapshot(state):
    if (not isinstance(state,BoundaryState) or not state.content_fingerprint
            or state.content_fingerprint!=_boundary_content_fingerprint(state)):
        raise InputError('Boundary snapshot contents changed or lack an integrity seal; rebuild the boundary')


def _fingerprint(system,bundle,settings=None,source_boundary_amplitudes=None):
    # This is a boundary-cache key, NOT a trace-time provenance seal. Every
    # caller must first verify the existing bundle provenance against current
    # inputs; constructing a new key cannot make stale rays current.
    amp=None if source_boundary_amplitudes is None else np.asarray(source_boundary_amplitudes,dtype='<f8').tobytes()
    payload=(trace_configuration_fingerprint(system,bundle.field,bundle.wavelength_nm,settings),
             bundle_content_fingerprint(bundle),amp)
    return hashlib.sha256(repr(payload).encode()).hexdigest()


def _full_surface_projection_upper(surface,axis):
    """Conservative bound over the entire rotational clear disk, not ray samples."""
    radius=surface.clear_radius;c=surface.curvature;k=surface.conic
    discriminant=1-(1+k)*c*c*radius*radius
    if discriminant<=0:
        raise ModelApplicabilityError('Last clear aperture reaches a singular/undefined sag domain')
    conic_edge=c*radius*radius/(1+np.sqrt(discriminant))
    # Every monomial's positive contribution is bounded independently.
    sag_upper=max(0.,conic_edge)+sum(max(0.,a)*radius**(2*(j+1)) for j,a in enumerate(surface.coefficients))
    if axis[2]<=0:
        raise ModelApplicabilityError('Backward chief unsupported by forward boundary profile')
    return float(axis[2]*(surface.z+sag_upper)+radius*np.linalg.norm(axis[:2]))


def _signed_double_area(points):
    a=points[:,1,:]-points[:,0,:];b=points[:,2,:]-points[:,0,:]
    return a[:,0]*b[:,1]-a[:,1]*b[:,0]


def make_boundary_state(system,bundle,settings=None,*,image_z=None,mode='focal',source_boundary_amplitudes=None,
                        jacobian_tolerance=1e-12,allow_unverified_bundle=False):
    """Build a conservative physical boundary with sampled topology checks.

    Clipping is represented by only triangles whose three STOP vertices trace
    successfully. This is a conservative sampled polygon, not an exact edge
    solve; ray-ring/arm refinement must establish its convergence. Engine rays
    must match their trace-time prescription. Hand-built analytic boundaries
    require ``allow_unverified_bundle=True`` and remain explicitly unverified.
    """
    settings=settings or EngineSettings()
    provenance_status=require_bundle_match(system,bundle,settings,allow_unverified=allow_unverified_bundle)
    if mode not in ('focal','afocal'): raise InputError('mode must be focal or afocal')
    chief=bundle.chief
    if not chief.success or not len(chief.positions):
        raise ModelApplicabilityError('Wave profile requires a transmitted chief ray')
    physical_failures={'aperture_blocked','escaped','tir','blocked','outside_aperture','total_internal_reflection'}
    for trace in bundle.traces:
        if not trace.success and trace.status not in physical_failures:
            raise ModelApplicabilityError(f'Numerical ray failure cannot become zero pupil: {trace.status}')
    valid=bundle.success_mask
    if valid.sum()<4:
        raise ModelApplicabilityError('Too few successful rays to reconstruct a pupil')
    stop=np.asarray(bundle.pupil_xy,float)
    if stop.shape!=(len(bundle.traces),2) or len(np.unique(stop,axis=0))!=len(stop):
        raise InputError('Pupil coordinates must be unique Nx2 STOP samples')
    e3=np.asarray(chief.last_direction,float)
    e1=np.array([1.,0.,0.])-e3[0]*e3
    if np.linalg.norm(e1)<=settings.ray.direction_tol:
        e1=np.array([0.,1.,0.])-e3[1]*e3
    e1=e1/np.linalg.norm(e1);e2=np.cross(e3,e1);basis=np.array([e1,e2,e3])
    bound=_full_surface_projection_upper(system.surfaces[-1],e3)
    clearance=settings.sampling.boundary_clearance_mm
    pchief=np.asarray(chief.last_position,float)
    delta=bound+clearance-np.dot(pchief,e3)
    origin=pchief+delta*e3
    n=system.index(system.surfaces[-1].medium_after,bundle.wavelength_nm)
    chief_s=chief.last_opl+n*delta
    uv=np.full((len(stop),2),np.nan);opls=np.full(len(stop),np.nan)
    for i,trace in enumerate(bundle.traces):
        if not valid[i]: continue
        d=trace.last_direction;p=trace.last_position
        cosine=np.dot(d,e3)
        if cosine<=settings.ray.direction_tol:
            raise ModelApplicabilityError('Grazing or backward boundary crossing')
        t=np.dot(origin-p,e3)/cosine
        if t<-settings.ray.intersection_tol_mm:
            raise ModelApplicabilityError('Canonical boundary is not after every final ray point')
        point=p+t*d
        uv[i]=(basis[:2]@(point-origin))
        opls[i]=trace.last_opl+n*t
    try: topology=Delaunay(stop).simplices
    except Exception as exc: raise ModelApplicabilityError('Degenerate STOP sample topology') from exc
    tri=topology[np.all(valid[topology],axis=1)]
    if not len(tri): raise ModelApplicabilityError('No fully transmitted pupil triangles')
    area0=_signed_double_area(stop[tri]);area1=_signed_double_area(uv[tri])
    jac=area1/area0
    scale=max(1.,float(np.max(np.abs(jac))))
    if np.any(np.abs(jac)<=jacobian_tolerance*scale) or np.any(np.sign(jac)!=np.sign(jac[0])):
        raise ModelApplicabilityError('Sampled ray map has a fold/caustic or singular Jacobian')
    reference=None
    z=system.image_z if image_z is None else image_z
    if mode=='focal':
        if z is None: raise InputError('Focal wave analysis requires image_z')
        if abs(e3[2])<=settings.ray.direction_tol: raise ModelApplicabilityError('Chief parallel to IMA')
        q0=pchief+(z-pchief[2])/e3[2]*e3
        reference=float(np.dot(q0-origin,e3))
        if reference<=0:
            raise ModelApplicabilityError('No regular boundary exists before the requested image plane')
    amp=None
    if source_boundary_amplitudes is not None:
        amp=np.asarray(source_boundary_amplitudes,float).copy()
        if amp.shape!=(len(stop),) or np.any(~np.isfinite(amp[valid])) or np.any(amp[valid]<0):
            raise InputError('Source boundary amplitudes must be nonnegative per traced ray')
        amp.setflags(write=False)
    return _seal_boundary_state(BoundaryState(origin,basis,uv,opls,float(chief_s),valid,tri,stop,jac,
                         float(bundle.wavelength_nm),float(n),reference,pchief,e3,
                         _fingerprint(system,bundle,settings,amp),bound,clearance,amp,provenance_status))


def _mapped_interpolation(vertices,triangles,values,x,y,*,overlap_tolerance=1e-10):
    """Piecewise linear interpolation preserving original STOP triangle topology."""
    xx,yy=np.meshgrid(x,y);points=np.column_stack([xx.ravel(),yy.ravel()])
    out=np.zeros((len(points),values.shape[1]));mask=np.zeros(len(points),bool)
    owner=np.full(len(points),-1,int)
    for ordinal,tri in enumerate(triangles):
        v=vertices[tri];a=v[1]-v[0];b=v[2]-v[0]
        det=a[0]*b[1]-a[1]*b[0]
        if det==0: raise ModelApplicabilityError('Degenerate mapped triangle')
        candidate=np.flatnonzero(np.all((points>=v.min(axis=0)-overlap_tolerance)&(points<=v.max(axis=0)+overlap_tolerance),axis=1))
        d=points[candidate]-v[0]
        w1=(d[:,0]*b[1]-d[:,1]*b[0])/det
        w2=(a[0]*d[:,1]-a[1]*d[:,0])/det
        w0=1-w1-w2
        inside=(w0>=-overlap_tolerance)&(w1>=-overlap_tolerance)&(w2>=-overlap_tolerance)
        inds=candidate[inside]
        strict=(w0[inside]>overlap_tolerance)&(w1[inside]>overlap_tolerance)&(w2[inside]>overlap_tolerance)
        if np.any(mask[inds]&strict):
            raise ModelApplicabilityError('Overlapping pupil map branches detected on wave grid')
        vals=w0[inside,None]*values[tri[0]]+w1[inside,None]*values[tri[1]]+w2[inside,None]*values[tri[2]]
        out[inds]=vals;mask[inds]=True;owner[inds]=ordinal
    return out.reshape(len(y),len(x),-1),mask.reshape(len(y),len(x)),owner.reshape(len(y),len(x))


def analyze_wave(system,bundle,settings=None,*,image_z=None,mode='focal',boundary_state=None,
                 grid_size=None,padding=None,frequencies=None,source_boundary_amplitudes=None,
                 illumination_policy='hold_boundary_field_amplitude',grid_margin=1.05,
                 phase_step_limit_rad=np.pi,direct_check_points=None,allow_unverified_bundle=False):
    """Lens-only real-ray scalar diffraction, with conservative sampled support.

    Return arrays and an explicit ``status``/diagnostic record. ``unresolved``
    means phase sampling or local propagation validity has not been established;
    numerical output is never labelled Zemax-equivalent. Per-ray thermal or RGB
    data are not required for a monochromatic bundle.

    Source amplitudes, if supplied, are physical electric-field amplitudes or
    consistently scaled amplitudes on the physical boundary plane, NOT sqrt of
    ray-packet power. They already include physical transport and transmission.
    Hand-built boundary studies require ``allow_unverified_bundle=True`` on
    every call, including cached-state reuse. This cannot override a mismatched
    engine trace. ``trace_provenance_status`` distinguishes the two origins.
    """
    settings=settings or EngineSettings()
    provenance_status=require_bundle_match(system,bundle,settings,allow_unverified=allow_unverified_bundle)
    N=settings.sampling.wave_grid_size if grid_size is None else grid_size
    pad=settings.sampling.fft_padding if padding is None else padding
    if int(N)!=N or N<8 or N%2: raise InputError('wave grid_size must be even and >=8')
    if grid_margin<=1: raise InputError('grid_margin must be >1 to contain sampled support')
    if illumination_policy not in ('hold_boundary_field_amplitude','uniform_fourier_pupil_each_focus'):
        raise InputError('Unknown focus illumination policy')
    if mode not in ('focal','afocal'): raise InputError('mode must be focal or afocal')
    state=boundary_state or make_boundary_state(system,bundle,settings,image_z=image_z,mode=mode,
                                               source_boundary_amplitudes=source_boundary_amplitudes,
                                               allow_unverified_bundle=allow_unverified_bundle)
    _validate_boundary_snapshot(state)
    if state.fingerprint!=_fingerprint(system,bundle,settings,state.source_amplitudes):
        raise InputError('Boundary snapshot does not match geometry, field, wavelength or traced rays')
    if boundary_state is not None and source_boundary_amplitudes is not None:
        raise InputError('Boundary reuse cannot silently change its illumination')
    if boundary_state is not None and ((mode=='focal')!=(state.reference_distance_mm is not None)):
        raise InputError('Cross-mode boundary reuse would change its illumination convention; rebuild explicitly')
    good=state.valid_ray_mask;uv=state.uv_mm;n=state.output_index;lam=state.wavelength_nm*1e-6
    z=system.image_z if image_z is None else image_z
    if mode=='focal':
        if z is None: raise InputError('Focal calculation requires image_z')
        q0=state.chief_position+(z-state.chief_position[2])/state.chief_direction[2]*state.chief_direction
        L=float(np.dot(q0-state.origin,state.chief_direction))
        if L<=0: raise ModelApplicabilityError('Image plane must remain downstream of the fixed boundary')
        r=np.sqrt(L*L+np.sum(uv*uv,axis=1))
        vertices=uv/r[:,None]
        opd=state.opl_mm-state.chief_opl_mm+n*np.sum(uv*uv,axis=1)/(r+L)
        reference=state.reference_distance_mm
        if reference is None and state.source_amplitudes is None:
            raise InputError('Afocal uniform-A snapshot has no focal B reference; construct a focal snapshot')
    else:
        L=None;q0=None;vertices=uv;opd=state.opl_mm-state.chief_opl_mm
    extent=np.max(np.abs(vertices[good]),axis=0)*grid_margin
    if np.any(extent<=0): raise ModelApplicabilityError('Degenerate pupil extent')
    dx,dy=2*extent/int(N)
    x=(np.arange(int(N))-int(N)//2)*dx;y=(np.arange(int(N))-int(N)//2)*dy
    amp=state.source_amplitudes if state.source_amplitudes is not None else np.ones(len(uv))
    values=np.column_stack([opd,amp])
    interpolated,mask,owner=_mapped_interpolation(vertices,state.triangles,values,x,y)
    w=interpolated[:,:,0];b=interpolated[:,:,1]
    xx,yy=np.meshgrid(x,y)
    if mode=='focal':
        gamma2=1-xx*xx-yy*yy
        if np.any(mask&(gamma2<=0)):
            raise ModelApplicabilityError('Mapped pupil exceeds reference hemisphere')
        g2=np.where(mask,gamma2,1.)
        uu=L*xx/np.sqrt(g2);vv=L*yy/np.sqrt(g2)
        if state.source_amplitudes is not None:
            b=b*L/g2
        elif illumination_policy=='hold_boundary_field_amplitude':
            # Exact analytic A of the original uniform Fourier-pupil convention.
            a=reference/(reference*reference+uu*uu+vv*vv)
            b=a*L/g2
        else:
            b=np.ones_like(w)
    pupil=np.where(mask,b*np.exp(2j*np.pi*w/lam),0j)
    if not np.any(mask): raise ModelApplicabilityError('No transmitted area on reconstructed wave grid')
    if mode=='focal':
        result=fft_psf(pupil,(dx,dy),lam,index=n,padding=pad)
        dl=fft_psf(np.abs(pupil),(dx,dy),lam,index=n,padding=pad)
        f=np.linspace(0,settings.sampling.max_frequency_lp_mm,settings.sampling.frequency_samples) if frequencies is None else np.asarray(frequencies,float)
        support=np.ptp(vertices[good],axis=0)/(lam/n)
        # Grid midpoint boundary rasterization can extend half a cell beyond
        # geometric vertices. Use the discrete pupil's true correlation bound.
        iy,ix=np.where(mask);support=np.array([(ix.max()-ix.min()+1)*dx,(iy.max()-iy.min()+1)*dy])/(lam/n)
        result['frequencies_lp_mm']=f
        result['otf_s']=sample_axis_otf(result,f,'x',support[0]);result['otf_t']=sample_axis_otf(result,f,'y',support[1])
        result['mtf_s']=np.abs(result['otf_s']);result['mtf_t']=np.abs(result['otf_t'])
        result['diffraction_mtf_s']=np.abs(sample_axis_otf(dl,f,'x',support[0]))
        result['diffraction_mtf_t']=np.abs(sample_axis_otf(dl,f,'y',support[1]))
        result['same_pupil_diffraction_psf']=dl['psf']
        # Quantify the actual FFT output window's local-expansion residual,
        # not a universal accuracy verdict or an arbitrary NA threshold.
        corners=np.array([[result['x_mm'][i],result['y_mm'][j]] for i in (0,-1) for j in (0,-1)])
        u=uu[mask];v=vv[mask];R=np.sqrt(L*L+u*u+v*v);chi=u/R;eta=v/R
        phase_residual=0.;amplitude_residual=0.
        for qx,qy in corners:
            rq=np.sqrt((u-qx)**2+(v-qy)**2+L*L)
            phase_residual=max(phase_residual,float(np.max(np.abs(rq-R+chi*qx+eta*qy)))*2*np.pi*n/lam)
            amplitude_residual=max(amplitude_residual,float(np.max(np.abs(R*R/(rq*rq)-1))))
        result['local_expansion_diagnostics']={'full_window_phase_residual_rad':phase_residual,
                                               'full_window_relative_amplitude_residual':amplitude_residual,
                                               'not_an_mtf_error_bound':True}
        if direct_check_points is not None:
            points=np.atleast_2d(np.asarray(direct_check_points,float))
            exact=focal_direct_reference(pupil,x,y,L,lam,points,index=n)
            cc,ee=np.meshgrid(x,y)
            local=np.array([np.sum(pupil*np.exp(-2j*np.pi*n/lam*(qx*cc+qy*ee)))*dx*dy/(1j*(lam/n)) for qx,qy in points])
            result['direct_scalar_check']={'points_xy_mm':points,'exact_intensity':np.abs(exact)**2,
                                          'local_fft_integral_intensity':np.abs(local)**2,
                                          'shared_boundary_input':True}
    else:
        result=afocal_fft(pupil,(dx,dy),lam,index=n,padding=pad)
        dl=afocal_fft(np.abs(pupil),(dx,dy),lam,index=n,padding=pad)
        result['same_pupil_diffraction_h_solid_angle']=dl['h_solid_angle']
        if frequencies is not None:
            f=np.asarray(frequencies,float)
            result['frequencies_cycles_rad']=f
            _angular_mtf_axes(result,dl,f,settings)
    phase_steps=[]
    for axis in (0,1):
        d=np.diff(w,axis=axis)*2*np.pi/lam
        m=np.take(mask,range(mask.shape[axis]-1),axis=axis)&np.take(mask,range(1,mask.shape[axis]),axis=axis)
        if np.any(m): phase_steps.append(float(np.max(np.abs(d[m]))))
    triangle_w=opd[state.triangles]
    ray_phase_span=float(np.max(np.ptp(triangle_w,axis=1))*2*np.pi/lam)
    grid_phase=max(phase_steps,default=0.)
    phase_ok=max(grid_phase,ray_phase_span)<phase_step_limit_rad
    padding_ok=result.get('linear_autocorrelation_supported',True)
    angular_ok=result.get('angular_frequency_sampling',{}).get('necessary_check_passed',True)
    result.update(status='conditional' if phase_ok and padding_ok and angular_ok else 'underresolved',profile_id=PROFILE_ID,
                  model_validity='scalar_profile_requires_convergence_and_local_propagation_check',
                  vendor_equivalence='unverified',mode=mode,boundary_state=state,
                  boundary_origin_mm=state.origin,output_basis=state.basis,
                  axis_convention='S=f_x along projected object-local X; T=f_y along e3 cross e1; not meridional/radial labels',
                  image_reference_point_mm=q0,distance_mm=L,
                  pupil=pupil,pupil_mask=mask,pupil_opd_mm=w,pupil_x=x,pupil_y=y,
                  pupil_spacing=(dx,dy),interpolation_owner=owner,
                  illumination='source_transport_boundary_amplitude' if state.source_amplitudes is not None else 'lens_metric',
                  illumination_policy=illumination_policy,removed_phase_terms=['piston'],
                  sampling_diagnostics={'max_grid_phase_step_rad':grid_phase,
                                        'max_traced_triangle_phase_span_rad':ray_phase_span,
                                        'phase_step_limit_rad':phase_step_limit_rad,
                                        'phase_sampling_necessary_check_passed':phase_ok,
                                        'angular_frequency_sampling_necessary_check_passed':angular_ok,
                                        'linear_autocorrelation_sampling_check_passed':padding_ok,
                                        'actual_traced_ray_count':len(bundle.traces),
                                        'successful_ray_count':int(good.sum()),
                                        'support':'conservative_transmitted_STOP_triangles',
                                        'edge_model':'binary_midpoint_polygon; refine rings/arms and wave grid independently',
                                        'interpolation':'piecewise_linear_unwrapped_OPD',
                                        'mapping_check':'sampled_triangle_orientation_and_grid_overlap; not a global proof',
                                        'sampled_boundary_jacobian_min':float(np.min(state.sampled_jacobians_mm2)),
                                        'sampled_boundary_jacobian_max':float(np.max(state.sampled_jacobians_mm2))})
    result['trace_provenance_status']=provenance_status
    result.update(_mtf_axis_metadata(state.basis))
    result['field']={'x':bundle.field.x,'y':bundle.field.y,'kind':bundle.field.kind}
    return result
