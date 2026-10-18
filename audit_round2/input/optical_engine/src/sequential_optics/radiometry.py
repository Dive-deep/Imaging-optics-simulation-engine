"""Direct-path radiometry with explicit source/measure and loss conventions.

Geometry is mm, vacuum wavelength nm. Absolute radiance is W/(mm^2 sr).
No source radiance is required by geometry or lens-only MTF. Missing absorption
is an explicit optional transparent assumption, never a hidden material value.
Reference: calculation_spec_ko.md R.1--R.5; Fresnel power and L/n^2 invariance.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.stats import qmc
from .errors import InputError, MissingDataError, NumericalError, ModelApplicabilityError
from .models import Ray, OpticalSystem
from .settings import EngineSettings


@dataclass(frozen=True)
class RectangularSource:
    """Continuous emitting panel; pixel metadata never modifies lens MTF.

    radiance=None means a relative active-area radiance of one. If basis is
    'panel_average', fill has already been applied and is not multiplied again.
    Source normal is +z and the emitted domain is the forward hemisphere.
    """
    z_mm: float
    width_mm: float = 12.8
    height_mm: float = 3.2
    fill_factor: float = .64
    radiance: float | None = None
    radiance_basis: str = 'active_area'
    pixels: tuple = (320, 80)

    def __post_init__(self):
        if not all(np.isfinite(v) for v in (self.z_mm,self.width_mm,self.height_mm,self.fill_factor)):
            raise InputError('Source geometry must be finite')
        if self.width_mm<=0 or self.height_mm<=0 or not 0<self.fill_factor<=1:
            raise InputError('Positive source dimensions and 0 < fill <= 1 are required')
        if self.radiance is not None and (not np.isfinite(self.radiance) or self.radiance<0):
            raise InputError('Source radiance must be finite and nonnegative')
        if self.radiance_basis not in ('active_area','panel_average'):
            raise InputError('radiance_basis must be active_area or panel_average')
        if len(self.pixels)!=2 or any(int(v)!=v or v<=0 for v in self.pixels):
            raise InputError('Source pixel metadata must contain two positive integers')

    @property
    def area_mm2(self): return self.width_mm*self.height_mm
    @property
    def panel_radiance(self):
        value=1. if self.radiance is None else self.radiance
        return value*self.fill_factor if self.radiance_basis=='active_area' else value
    @property
    def hemisphere_power(self): return np.pi*self.area_mm2*self.panel_radiance
    @property
    def is_absolute(self): return self.radiance is not None


def fresnel_power(n1,n2,cos_incidence):
    """Real lossless uncoated interface; power T, not field amplitude t.

    Returned gs/gp are positive power-normalized transmitted amplitudes.
    TIR is a physical zero transmitted branch, not numerical failure.
    """
    if not all(np.isfinite(v) for v in (n1,n2,cos_incidence)) or n1<=0 or n2<=0 or not 0<=cos_incidence<=1:
        raise InputError('Fresnel requires positive real indices and cosine in [0,1]')
    ci=float(cos_incidence)
    if n1==n2: return dict(Ts=1.,Tp=1.,Rs=0.,Rp=0.,gs=1.,gp=1.,cos_transmission=ci,tir=False)
    st2=(n1/n2)**2*max(0.,1-ci*ci)
    if st2>=1:
        return dict(Ts=0.,Tp=0.,Rs=1.,Rp=1.,gs=0.,gp=0.,cos_transmission=0.,tir=True)
    ct=math.sqrt(1-st2)
    if ci==0:
        if n1==n2: return dict(Ts=1.,Tp=1.,Rs=0.,Rp=0.,gs=1.,gp=1.,cos_transmission=0.,tir=False)
        return dict(Ts=0.,Tp=0.,Rs=1.,Rp=1.,gs=0.,gp=0.,cos_transmission=ct,tir=False)
    rs=(n1*ci-n2*ct)/(n1*ci+n2*ct)
    rp=(n2*ci-n1*ct)/(n2*ci+n1*ct)
    ts=1-rs*rs;tp=1-rp*rp
    return dict(Ts=ts,Tp=tp,Rs=rs*rs,Rp=rp*rp,gs=math.sqrt(ts),gp=math.sqrt(tp),cos_transmission=ct,tir=False)


def _basis(direction):
    d=np.asarray(direction,dtype=float)
    reference=np.eye(3)[np.argmin(np.abs(d))]
    a=np.cross(d,reference);a/=np.linalg.norm(a)
    return np.column_stack((a,np.cross(d,a)))


def _absorption(system,medium,wavelength,missing_absorption,assumptions):
    try:
        mu=system.materials[medium].absorption(wavelength,system.temperature_c)
    except MissingDataError:
        if missing_absorption!='transparent': raise
        assumptions.add(f'absorption_missing_assumed_zero:{medium}')
        return 0.
    if not np.isfinite(mu) or mu<0: raise InputError('Absorption must be finite nonnegative per mm')
    return float(mu)


def trace_transmission(system,ray,trace,wavelength_nm=None,settings=None,
                       mode='unpolarized_jones_power',missing_absorption='transparent'):
    """Physical attenuation along a *successful full* geometry trace.

    Jones columns represent two incoherent unit-power input states. The
    3-vector representation transports their S/P basis at every interface.
    This computes power only, not vector coherent diffraction or coating phase.
    The returned reflection/absorption ledger is for this successful path.
    """
    from .geometry import gradient
    settings=settings or EngineSettings();wl=settings.wavelength(wavelength_nm)
    modes=('unpolarized_jones_power','scalar_interface_average','lossless_scalar')
    if mode not in modes or missing_absorption not in ('transparent','error'):
        raise InputError('Unknown transmission/absorption policy')
    if not trace.success or len(trace.positions)!=len(system.surfaces):
        raise InputError('trace_transmission requires a successful full system trace')
    direction=ray.direction.copy();point=ray.origin.copy();medium=system.ambient_material
    electric=_basis(direction).astype(complex);tau=1.;absorbed=0.;reflected=0.;assumptions=set()
    if mode=='scalar_interface_average': assumptions.add('interface_average_resets_polarization_each_surface')
    if mode=='lossless_scalar': assumptions.add('surface_and_bulk_loss_disabled')
    for surface,hit,out_direction in zip(system.surfaces,trace.positions,trace.directions):
        distance=float(np.linalg.norm(hit-point))
        if mode!='lossless_scalar':
            mu=_absorption(system,medium,wl,missing_absorption,assumptions)
            survival=math.exp(-mu*distance);absorbed+=tau*(1-survival)
            tau*=survival;electric*=math.sqrt(survival)
            g=np.asarray(gradient(surface,hit[0],hit[1]),dtype=float)
            normal=np.array([-g[0],-g[1],1.]);normal/=np.linalg.norm(normal)
            ci=abs(float(direction@normal));coeff=fresnel_power(system.index(medium,wl),system.index(surface.medium_after,wl),min(1.,ci))
            before=tau
            if mode=='scalar_interface_average':
                tau*=.5*(coeff['Ts']+coeff['Tp'])
            else:
                s=np.cross(direction,normal)
                if np.linalg.norm(s)<=settings.ray.direction_tol:
                    # At normal incidence Ts=Tp: minimal-change transverse basis.
                    s=_basis(direction)[:,0]
                else: s/=np.linalg.norm(s)
                pi=np.cross(direction,s);po=np.cross(out_direction,s)
                es=s@electric;ep=pi@electric
                electric=np.outer(s,coeff['gs']*es)+np.outer(po,coeff['gp']*ep)
                tau=float(.5*np.sum(np.abs(electric)**2))
            reflected+=max(0.,before-tau)
        point=np.asarray(hit);direction=np.asarray(out_direction);medium=surface.medium_after
    return dict(transmission=tau,absorbed_fraction=absorbed,reflected_fraction=reflected,
                ledger_residual=1-tau-absorbed-reflected,assumptions=sorted(assumptions),
                mode=mode,wavelength_nm=wl)


def _trace_outcome(system,ray,wl,settings,mode,missing_absorption):
    """Only explicitly physical geometry statuses may contribute zero power."""
    from .geometry import trace_ray
    trace=trace_ray(system,ray,wl,settings=settings)
    if trace.success:
        transport=trace_transmission(system,ray,trace,wl,settings,mode,missing_absorption)
        return trace,transport['transmission'],transport['assumptions']
    physical=('aperture_blocked','escaped','tir','blocked','vignetted','missed','miss','total_internal_reflection','no_intersection')
    if trace.status not in physical:
        raise NumericalError(f'Radiometry refuses unresolved trace {trace.status}: {trace.message}')
    return trace,0.,[]


def cosine_hemisphere_samples(count,seed=0):
    """Deterministic scrambled Sobol directions, p(omega)=cos(theta)/pi.

    Any count is supported by generating the next power of two then slicing;
    use powers of two for the normal Sobol balance guarantee.
    """
    if int(count)!=count or count<=0: raise InputError('sample_count must be a positive integer')
    uv=qmc.Sobol(2,scramble=True,seed=seed).random_base2(int(math.ceil(math.log2(count))))[:count]
    r=np.sqrt(uv[:,0]);phi=2*np.pi*uv[:,1]
    return np.column_stack((r*np.cos(phi),r*np.sin(phi),np.sqrt(1-r*r)))


def source_direction_samples(source_xy_mm,source_z_mm,count,seed=0,*,system=None,
                             launch_plane_z_mm=None,launch_radius_mm=None):
    """Directions and projected-solid-angle quadrature weights per source point.

    Without a plane: cosine-weighted hemisphere, weights pi/N. With a plane:
    uniform disk before all first-surface material, weights
    A_disk cos_source cos_plane/(N r^2). This is not STOP sampling.
    The disk must cover the convex hull of source point and first clear aperture,
    guaranteeing omitted directions cannot hit the first optical clear surface.
    """
    xy=np.asarray(source_xy_mm,dtype=float)
    if xy.shape!=(2,) or not np.all(np.isfinite(xy)) or not np.isfinite(source_z_mm): raise InputError('Invalid source point')
    if launch_plane_z_mm is None and launch_radius_mm is None:
        return cosine_hemisphere_samples(count,seed),np.full(count,np.pi/count),['cosine_hemisphere']
    if launch_plane_z_mm is None or launch_radius_mm is None or system is None:
        raise InputError('Plane importance requires system, launch_plane_z_mm and launch_radius_mm')
    from .geometry import sag
    first=system.surfaces[0];a=first.clear_radius
    # Triangle bound valid for radial conic sag and finite even polynomial.
    conic_only=replace(first,coefficients=())
    excursion=abs(float(sag(conic_only,a,0.)))+sum(abs(c)*a**(2*j+2) for j,c in enumerate(first.coefficients))
    min_z_bound=first.z-excursion
    if not source_z_mm<launch_plane_z_mm<min_z_bound:
        raise InputError('Launch plane must lie between source and conservative first-surface lower bound')
    if not np.isfinite(launch_radius_mm) or launch_radius_mm<max(a,float(np.linalg.norm(xy))):
        raise InputError('Launch disk must enclose source point and first clear aperture for complete coverage')
    if int(count)!=count or count<=0: raise InputError('sample_count must be positive integer')
    uv=qmc.Sobol(2,scramble=True,seed=seed).random_base2(int(math.ceil(math.log2(count))))[:count]
    radius=launch_radius_mm*np.sqrt(uv[:,0]);phi=2*np.pi*uv[:,1]
    vectors=np.column_stack((radius*np.cos(phi)-xy[0],radius*np.sin(phi)-xy[1],np.full(count,launch_plane_z_mm-source_z_mm)))
    distance=np.linalg.norm(vectors,axis=1);directions=vectors/distance[:,None]
    weights=np.pi*launch_radius_mm**2/count*directions[:,2]**2/distance**2
    return directions,weights,['entrance_plane_importance_complete_clear_aperture_coverage']


def projected_solid_angle_aperture(image_xyz,aperture_z_mm,aperture_radius_mm,radial_order=16,azimuth_count=64):
    """Integrate cos(theta_image)dOmega for a coaxial circular plane aperture.

    A deterministic reference fixture, not a general lens cos^4 model.
    Gauss integration in r^2 and periodic azimuth; image normal parallel z.
    """
    x=np.asarray(image_xyz,dtype=float)
    if x.shape!=(3,) or not np.all(np.isfinite(x)) or aperture_radius_mm<=0 or x[2]==aperture_z_mm:
        raise InputError('Invalid aperture or image point')
    if int(radial_order)!=radial_order or radial_order<=0 or int(azimuth_count)!=azimuth_count or azimuth_count<=0:
        raise InputError('Quadrature orders must be positive integers')
    nodes,weights=np.polynomial.legendre.leggauss(radial_order)
    r=aperture_radius_mm*np.sqrt((nodes+1)/2)
    phi=2*np.pi*(np.arange(azimuth_count)+.5)/azimuth_count
    xx=(r[:,None]*np.cos(phi)[None,:]-x[0]);yy=(r[:,None]*np.sin(phi)[None,:]-x[1])
    dz=x[2]-aperture_z_mm;rr=xx*xx+yy*yy+dz*dz
    weights_area=weights[:,None]*aperture_radius_mm**2/4*(2*np.pi/azimuth_count)
    return float(np.sum(dz*dz/rr**2*weights_area))


def _reverse_system(system):
    """Mirror z after reversing surface order to retain forward tracer convention."""
    media=[system.ambient_material]+[s.medium_after for s in system.surfaces]
    surfaces=tuple(replace(s,z=-s.z,curvature=-s.curvature,
                          coefficients=tuple(-a for a in s.coefficients),medium_after=media[j])
                   for j,s in reversed(list(enumerate(system.surfaces))))
    return replace(system,surfaces=surfaces,object_z=None if system.image_z is None else -system.image_z,
                   image_z=None if system.object_z is None else -system.object_z,
                   stop_index=len(surfaces)-1-system.stop_index,ambient_material=media[-1])


def uniform_scene_ri(system,image_points_xy,wavelength_nm=None,settings=None,*,
                     image_z_mm=None,source_radiance=None,direction_grid_size=96,
                     direction_cosine_bounds=(-1.,1.,-1.,1.),mode='unpolarized_jones_power',
                     missing_absorption='transparent',axis_index=None,source_bounds_mm=None):
    """Actual reverse-ray irradiance integral E=(ni/ns)^2 int tau L dl dm.

    Input evaluation positions are physical image-plane coordinates, usually
    chief-ray intersections. Full direction disk is the safe domain. A smaller
    rectangle is permitted but is explicitly reported as truncated angular
    domain, and its coverage must be established by the caller. A finite source
    rectangle changes the output meaning to finite_source_irradiance.
    No extra image cosine and no cos^4 correction are applied.
    """
    settings=settings or EngineSettings();wl=settings.wavelength(wavelength_nm)
    xy=np.atleast_2d(np.asarray(image_points_xy,dtype=float))
    if xy.shape[1]!=2 or not np.all(np.isfinite(xy)): raise InputError('Image points must be finite N x 2')
    z=system.image_z if image_z_mm is None else image_z_mm
    if z is None: raise MissingDataError('Finite image_z is required for image-plane RI')
    if int(direction_grid_size)!=direction_grid_size or direction_grid_size<2: raise InputError('RI grid requires integer size >= 2')
    l0,l1,m0,m1=map(float,direction_cosine_bounds)
    if not (-1<=l0<l1<=1 and -1<=m0<m1<=1): raise InputError('Direction cosine rectangle must be within [-1,1]')
    L=1. if source_radiance is None else float(source_radiance)
    if not np.isfinite(L) or L<0: raise InputError('Radiance must be finite nonnegative')
    if source_bounds_mm is not None and system.object_z is None: raise MissingDataError('Finite source mask requires object_z')
    reverse=_reverse_system(system)
    ll=l0+(np.arange(direction_grid_size)+.5)*(l1-l0)/direction_grid_size
    mm=m0+(np.arange(direction_grid_size)+.5)*(m1-m0)/direction_grid_size
    l,m=np.meshgrid(ll,mm);inside=l*l+m*m<1
    directions=np.column_stack((l[inside],m[inside],np.sqrt(1-l[inside]**2-m[inside]**2)))
    measure=(l1-l0)*(m1-m0)/direction_grid_size**2
    index_factor=(system.index(system.surfaces[-1].medium_after,wl)/system.index(system.ambient_material,wl))**2
    irradiance=[];accepted=[];assumptions=set();counts={}
    if (l0,l1,m0,m1)!=(-1.,1.,-1.,1.): assumptions.add('truncated_image_direction_domain_requires_coverage_check')
    for x,y in xy:
        total=0.;hits=0
        for d in directions:
            ray=Ray(np.array([x,y,-z]),d)
            trace,tau,notes=_trace_outcome(reverse,ray,wl,settings,mode,missing_absorption)
            counts[trace.status]=counts.get(trace.status,0)+1;assumptions.update(notes)
            if not trace.success: continue
            if system.object_z is not None:
                distance=(-system.object_z-trace.last_position[2])/trace.last_direction[2]
                if distance<0: raise ModelApplicabilityError('Object plane is not upstream of all traced surfaces')
                q=trace.last_position+distance*trace.last_direction
                if source_bounds_mm is not None:
                    xmin,xmax,ymin,ymax=source_bounds_mm
                    if not (xmin<=q[0]<=xmax and ymin<=q[1]<=ymax): continue
                if mode!='lossless_scalar':
                    mu=_absorption(reverse,reverse.surfaces[-1].medium_after,wl,missing_absorption,assumptions)
                    tau*=math.exp(-mu*distance)
            total+=tau;hits+=1
        irradiance.append(L*index_factor*measure*total);accepted.append(hits)
    irradiance=np.asarray(irradiance);peak=float(np.max(irradiance))
    relative=np.full_like(irradiance,np.nan) if peak==0 else irradiance/peak
    result=dict(irradiance=irradiance,relative_illumination=relative,normalization='sampled_peak',
                sampled_peak=peak,image_points_xy=xy,image_z_mm=z,wavelength_nm=wl,
                units='W/mm^2' if source_radiance is not None else 'relative_irradiance',
                status='ok' if peak>0 else 'undefined_zero_illumination',accepted_rays=np.array(accepted),
                direction_samples=len(directions),direction_cell_measure=measure,trace_status_counts=counts,
                assumptions=sorted(assumptions),analysis='uniform_scene_ri' if source_bounds_mm is None else 'finite_source_irradiance')
    if axis_index is not None:
        if not 0<=axis_index<len(xy): raise InputError('axis_index is outside image points')
        denominator=irradiance[axis_index]
        result['relative_to_axis']=irradiance/denominator if denominator>0 else np.full_like(irradiance,np.nan)
    return result


def spherical_bin_solid_angles(theta_edges,phi_edges):
    theta=np.asarray(theta_edges,dtype=float);phi=np.asarray(phi_edges,dtype=float)
    if theta.ndim!=1 or phi.ndim!=1 or len(theta)<2 or len(phi)<2 or not np.all(np.diff(theta)>0) or not np.all(np.diff(phi)>0):
        raise InputError('Angular edges must be increasing 1D arrays')
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi)) or theta[0]<0 or theta[-1]>np.pi or phi[-1]-phi[0]>2*np.pi+1e-12:
        raise InputError('Angular edges exceed a single sphere')
    return (np.cos(theta[:-1])-np.cos(theta[1:]))[:,None]*np.diff(phi)[None,:]


def slope_rectangle_solid_angle(p0,p1,q0,q1):
    """Exact sr of a rectangular forward gnomonic bin, p=dx/dz,q=dy/dz."""
    if not all(np.isfinite(v) for v in (p0,p1,q0,q1)) or p0>=p1 or q0>=q1: raise InputError('Invalid slope rectangle')
    def primitive(p,q): return np.arctan2(p*q,np.sqrt(1+p*p+q*q))
    return float(primitive(p1,q1)-primitive(p0,q1)-primitive(p1,q0)+primitive(p0,q0))


def angle_chart_jacobian(alpha,beta):
    """dOmega/(dalpha dbeta), alpha=atan(dx/dz), beta=atan(dy/dz)."""
    a=np.asarray(alpha);b=np.asarray(beta)
    if np.any(np.abs(a)>=np.pi/2) or np.any(np.abs(b)>=np.pi/2): raise InputError('Forward angle chart requires |angle| < pi/2')
    p=np.tan(a);q=np.tan(b)
    return (1+p*p)*(1+q*q)/(1+p*p+q*q)**1.5


def far_field(system,source,theta_edges,phi_edges,wavelength_nm=None,settings=None,*,
              sample_count=4096,seed=0,mode='unpolarized_jones_power',missing_absorption='transparent',
              launch_plane_z_mm=None,launch_radius_mm=None):
    """Whole continuous finite-source geometric radiant intensity.

    A joint 4D Sobol sample covers source area and cosine-weighted hemisphere;
    every launched packet has pi A L/N power. system=None is a no-optics fixture.
    Optical power is never multiplied by an index-square factor. Spherical bin
    power is divided by exact sr, with no detector cosine. No diffraction blur.
    """
    settings=settings or EngineSettings();wl=settings.wavelength(wavelength_nm)
    solid=spherical_bin_solid_angles(theta_edges,phi_edges)
    ambient_notes=set()
    if system is not None and mode!='lossless_scalar':
        final_mu=_absorption(system,system.surfaces[-1].medium_after,wl,missing_absorption,ambient_notes)
        if final_mu>0:
            raise ModelApplicabilityError('Far-field at infinity requires nonabsorbing output medium; use a finite detector plane')
    if int(sample_count)!=sample_count or sample_count<=0: raise InputError('sample_count must be positive integer')
    uv=qmc.Sobol(4,scramble=True,seed=seed).random_base2(int(math.ceil(math.log2(sample_count))))[:sample_count]
    positions=np.column_stack(((uv[:,0]-.5)*source.width_mm,(uv[:,1]-.5)*source.height_mm,np.full(sample_count,source.z_mm)))
    r=np.sqrt(uv[:,2]);phi=2*np.pi*uv[:,3]
    dirs=np.column_stack((r*np.cos(phi),r*np.sin(phi),np.sqrt(1-r*r)))
    launch_weights=np.full(sample_count,source.hemisphere_power/sample_count)
    assumptions=set(ambient_notes);statuses={}
    if launch_plane_z_mm is not None or launch_radius_mm is not None:
        # Validate free-space and complete coverage at the source corner.
        _,_,notes=source_direction_samples((source.width_mm/2,source.height_mm/2),source.z_mm,1,seed,
                    system=system,launch_plane_z_mm=launch_plane_z_mm,launch_radius_mm=launch_radius_mm)
        assumptions.update(notes)
        target=np.column_stack((launch_radius_mm*r*np.cos(phi),launch_radius_mm*r*np.sin(phi),np.full(sample_count,launch_plane_z_mm)))
        vectors=target-positions;distances=np.linalg.norm(vectors,axis=1);dirs=vectors/distances[:,None]
        launch_weights=source.area_mm2*source.panel_radiance*np.pi*launch_radius_mm**2/sample_count*dirs[:,2]**2/distances**2
    output_directions=[];powers=[]
    for origin,d,packet in zip(positions,dirs,launch_weights):
        if system is None:
            output_directions.append(d);powers.append(packet);continue
        trace,tau,notes=_trace_outcome(system,Ray(origin,d),wl,settings,mode,missing_absorption)
        statuses[trace.status]=statuses.get(trace.status,0)+1;assumptions.update(notes)
        if trace.success:
            output_directions.append(trace.last_direction);powers.append(packet*tau)
    directions=np.asarray(output_directions).reshape(-1,3);powers=np.asarray(powers)
    theta=np.arccos(np.clip(directions[:,2],-1,1))
    # Wrap phi into the edge interval starting at phi_edges[0].
    phis=(np.arctan2(directions[:,1],directions[:,0])-phi_edges[0])%(2*np.pi)+phi_edges[0]
    histogram=np.histogram2d(theta,phis,bins=(theta_edges,phi_edges),weights=powers)[0]
    intensity=histogram/solid;peak=float(np.max(intensity));total_exit=float(np.sum(powers))
    return dict(bin_power=histogram,solid_angle_sr=solid,intensity=intensity,
                relative_far_field_intensity=intensity/peak if peak>0 else np.full_like(intensity,np.nan),
                source_hemisphere_power=source.hemisphere_power,exit_power=total_exit,
                nondirect_power=source.hemisphere_power-total_exit,
                launched_domain_power=float(np.sum(launch_weights)),
                not_launched_hemisphere_power=source.hemisphere_power-float(np.sum(launch_weights)),
                outside_reported_angular_roi_power=total_exit-float(np.sum(histogram)),
                direct_efficiency=total_exit/source.hemisphere_power if source.hemisphere_power>0 else np.nan,
                units='W/sr' if source.is_absolute else 'relative_intensity_per_sr',
                status='ok' if peak>0 else 'undefined_zero_intensity',wavelength_nm=wl,
                assumptions=sorted(assumptions),trace_status_counts=statuses,
                sample_count=sample_count,seed=seed,model='geometric_continuous_lambertian',
                ledger_note='nondirect aggregates reflection, absorption, missed and blocked; not an absorption estimate')
