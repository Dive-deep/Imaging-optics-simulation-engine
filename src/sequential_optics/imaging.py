"""Power-conserving finite-source image synthesis on a common physical grid.

The kernel is a per-cell *probability mass*, not a sampled PSF density. A kernel
whose stored mass sums below one explicitly loses its tail outside the ROI; no
renormalization hides that loss. Kernels are source-conditioned. Uniform lens
metric PSFs are not accepted as physical illumination without separate proof.
"""
from collections.abc import Mapping
from dataclasses import dataclass
import numpy as np
from .errors import InputError, MissingDataError, ModelApplicabilityError
from .models import Ray
from .settings import EngineSettings
from .radiometry import (RectangularSource, source_direction_samples, _trace_outcome,
                         _absorption)


_MASS_SUM_TOLERANCE = 1e-12  # dimensionless floating-point probability ledger tolerance


def _immutable_float_array(value):
    """Own an immutable byte-buffer copy; setflags(write=True) cannot unlock it.

    A read-only owning ndarray is insufficient: its owner can restore the write
    flag. A bytes-backed view also protects all ndarray base views exposed here.
    """
    a = np.asarray(value, dtype=float)
    return np.frombuffer(a.tobytes(order='C'), dtype=a.dtype).reshape(a.shape)


@dataclass(frozen=True)
class PixelGrid:
    x_edges_mm: np.ndarray
    y_edges_mm: np.ndarray

    def __post_init__(self):
        for key in ('x_edges_mm','y_edges_mm'):
            object.__setattr__(self,key,_immutable_float_array(getattr(self,key)))
        self.validate()

    def validate(self):
        """Check physical pixel measure before consuming an external grid."""
        for key in ('x_edges_mm','y_edges_mm'):
            a=np.asarray(getattr(self,key),dtype=float)
            if a.ndim!=1 or len(a)<2 or not np.all(np.isfinite(a)) or not np.all(np.diff(a)>0):
                raise InputError('Pixel edges must be finite strictly increasing 1D arrays')
        with np.errstate(over='ignore',invalid='ignore',under='ignore'):
            areas=self.cell_areas_mm2
        if not np.all(np.isfinite(areas)) or np.any(areas<=0):
            raise InputError('Pixel cell areas must be finite and strictly positive')

    @property
    def cell_areas_mm2(self): return np.diff(self.y_edges_mm)[:,None]*np.diff(self.x_edges_mm)[None,:]
    @property
    def shape(self): return (len(self.y_edges_mm)-1,len(self.x_edges_mm)-1)
    @property
    def x_centers_mm(self): return .5*(self.x_edges_mm[:-1]+self.x_edges_mm[1:])
    @property
    def y_centers_mm(self): return .5*(self.y_edges_mm[:-1]+self.y_edges_mm[1:])


@dataclass(frozen=True)
class SourceKernel:
    grid: PixelGrid
    probability_mass: np.ndarray
    wavelength_nm: float
    provenance: str = 'physical_source'
    center_convention: str = 'global_physical_coordinates'
    model: str = 'provided_source_conditioned_kernel'

    def __post_init__(self):
        object.__setattr__(self,'probability_mass',_immutable_float_array(self.probability_mass))
        self.validate()

    def validate(self):
        """Recheck kernel/grid invariants at the image-synthesis boundary."""
        if not isinstance(self.grid,PixelGrid):
            raise InputError('A source kernel requires a PixelGrid')
        self.grid.validate()
        mass=np.asarray(self.probability_mass,dtype=float)
        if mass.shape!=self.grid.shape or not np.all(np.isfinite(mass)) or np.any(mass<0):
            raise InputError('Kernel probability mass must match grid and be finite nonnegative')
        if np.sum(mass)>1+_MASS_SUM_TOLERANCE: raise InputError('Kernel cell probabilities sum above one; missing normalization/measure')
        if self.provenance!='physical_source':
            raise InputError('Image kernels must be source-conditioned; lens-metric pupil is not interchangeable')
        if self.center_convention!='global_physical_coordinates':
            raise InputError('Register kernels on common physical grid before synthesis')
        if not np.isfinite(self.wavelength_nm) or self.wavelength_nm<=0: raise InputError('Kernel wavelength must be positive nm')

    @property
    def captured_fraction(self): return float(np.sum(self.probability_mass))


def _same_grid(a,b):
    return np.array_equal(a.x_edges_mm,b.x_edges_mm) and np.array_equal(a.y_edges_mm,b.y_edges_mm)


def synthesize_image(kernels,collected_power_density,source_area_weights_mm2,*,bitmap_weights=None):
    """Sum Q(q) h_cell(q) dA_source; no additional warp Jacobian or Fresnel.

    Q already contains angular collection, source radiance, fill and attenuation.
    bitmap_weights is an optional *linear* radiance modulation, applied once.
    Different wavelengths require registered_spectral_sum, not hidden recentering.
    Only roundoff-scale positive power residuals are tolerated, and reported as
    power_balance_residual. Larger violations fail instead of clipping the loss.
    """
    kernels=tuple(kernels)
    if not kernels: raise InputError('At least one source kernel is required')
    for kernel in kernels:
        if not isinstance(kernel,SourceKernel): raise InputError('Image synthesis requires SourceKernel objects')
        kernel.validate()
    q=np.asarray(collected_power_density,dtype=float);areas=np.asarray(source_area_weights_mm2,dtype=float)
    if q.shape!=(len(kernels),) or areas.shape!=q.shape or not np.all(np.isfinite(q)) or not np.all(np.isfinite(areas)) or np.any(q<0) or np.any(areas<=0):
        raise InputError('Q and positive source-area weights must match kernel count')
    weights=np.ones_like(q) if bitmap_weights is None else np.asarray(bitmap_weights,dtype=float)
    if weights.shape!=q.shape or not np.all(np.isfinite(weights)) or np.any(weights<0): raise InputError('Bitmap modulation must be finite nonnegative')
    grid=kernels[0].grid;wl=kernels[0].wavelength_nm
    if any(not _same_grid(grid,k.grid) or k.wavelength_nm!=wl for k in kernels):
        raise InputError('Single-wavelength synthesis requires common physical grid and wavelength')
    with np.errstate(over='ignore',invalid='ignore'):
        input_power=q*areas*weights
        cell_power=sum((p*k.probability_mass for p,k in zip(input_power,kernels)),np.zeros(grid.shape))
        total=float(np.sum(input_power));inside=float(np.sum(cell_power))
    if not np.isfinite(total) or not np.isfinite(inside) or not np.all(np.isfinite(cell_power)):
        raise InputError('Image power arithmetic must remain finite; rescale input powers')
    if inside>total*(1.+_MASS_SUM_TOLERANCE):
        raise InputError('Detected ROI power exceeds collected power beyond probability tolerance')
    outside=max(0.,total-inside)
    with np.errstate(over='ignore',invalid='ignore'):
        irradiance=cell_power/grid.cell_areas_mm2
    if not np.all(np.isfinite(irradiance)):
        raise InputError('Image irradiance arithmetic must remain finite; rescale powers or pixel areas')
    return dict(grid=grid,cell_power=cell_power,irradiance=irradiance,
                collected_power=total,detected_roi_power=inside,outside_roi_power=outside,
                power_balance_residual=inside+outside-total,
                wavelength_nm=wl,model='incoherent_source_conditioned_sum',
                normalization='none; raw collected power and cropped tail retained')


def geometric_source_kernel(system,source_xy_mm,image_grid,wavelength_nm=None,settings=None,*,
                            source_z_mm=None,source_radiance=None,fill_factor=1.,
                            sample_count=1024,seed=0,image_z_mm=None,
                            launch_plane_z_mm=None,launch_radius_mm=None,
                            mode='unpolarized_jones_power',missing_absorption='transparent'):
    """Actual-ray geometric kernel and matched Q for one finite source point.

    Geometry directions retain the actual coma, distortion and lateral color.
    Q is per source area, includes fill exactly once, and includes final-medium
    absorption to the physical image plane. The sample histogram is geometric;
    no diffraction is claimed or silently added to an FFT lens metric PSF.
    """
    settings=settings or EngineSettings();wl=settings.wavelength(wavelength_nm)
    if not isinstance(image_grid,PixelGrid): raise InputError('Image detector requires a PixelGrid')
    image_grid.validate()
    z=system.object_z if source_z_mm is None else source_z_mm
    image_z=system.image_z if image_z_mm is None else image_z_mm
    if z is None or image_z is None: raise MissingDataError('Finite object_z and image_z required for geometric source image')
    L=1. if source_radiance is None else float(source_radiance)
    if not np.isfinite(L) or L<0 or not np.isfinite(fill_factor) or not 0<fill_factor<=1:
        raise InputError('Radiance must be nonnegative, and fill in (0,1]')
    d,measure,notes=source_direction_samples(source_xy_mm,z,sample_count,seed,system=system,
                                launch_plane_z_mm=launch_plane_z_mm,launch_radius_mm=launch_radius_mm)
    launch=L*fill_factor*measure;origin=np.array([*source_xy_mm,z],dtype=float)
    hits=[];powers=[];assumptions=set(notes);statuses={}
    for direction,packet in zip(d,launch):
        trace,tau,notes=_trace_outcome(system,Ray(origin,direction),wl,settings,mode,missing_absorption)
        assumptions.update(notes);statuses[trace.status]=statuses.get(trace.status,0)+1
        if not trace.success: continue
        if trace.last_direction[2]<=0: raise ModelApplicabilityError('Output ray does not propagate toward forward image plane')
        distance=(image_z-trace.last_position[2])/trace.last_direction[2]
        if distance<0: raise ModelApplicabilityError('Geometric detector must be downstream of the whole final traced surface')
        if mode!='lossless_scalar':
            mu=_absorption(system,system.surfaces[-1].medium_after,wl,missing_absorption,assumptions)
            tau*=np.exp(-mu*distance)
        hit=trace.last_position+distance*trace.last_direction
        hits.append(hit[:2]);powers.append(packet*tau)
    powers=np.asarray(powers);hits=np.asarray(hits).reshape(-1,2);Q=float(np.sum(powers))
    histogram=np.histogram2d(hits[:,1],hits[:,0],bins=(image_grid.y_edges_mm,image_grid.x_edges_mm),weights=powers)[0]
    mass=histogram/Q if Q>0 else np.zeros(image_grid.shape)
    kernel=SourceKernel(image_grid,mass,wl,model='geometric_real_ray')
    return dict(kernel=kernel,collected_power_density=Q,source_hemisphere_power_density=np.pi*L*fill_factor,
                launched_domain_power_density=float(np.sum(launch)),
                captured_fraction=kernel.captured_fraction,source_xy_mm=np.asarray(source_xy_mm),
                assumptions=sorted(assumptions),trace_status_counts=statuses,
                status='ok' if Q>0 else 'zero_collected_power',sample_count=sample_count,seed=seed,
                units='W/mm^2_source' if source_radiance is not None else 'relative_power_per_source_mm2')


def geometric_image(system,bitmap,source,image_grid,wavelength_nm=None,settings=None,*,
                    samples_per_source=512,seed=0,image_z_mm=None,
                    launch_plane_z_mm=None,launch_radius_mm=None,
                    mode='unpolarized_jones_power',missing_absorption='transparent',
                    bitmap_encoding='linear'):
    """Finite panel bitmap via actual-ray macro-cell midpoint source quadrature.

    Bitmap shape is a numerical macro-scene grid, independent of LED pixel
    metadata. Increase macro resolution as well as rays/cell for convergence.
    Coordinates: columns x increasing left-to-right, rows y increasing bottom-to-top.
    No unspecified gamma, exposure, or normalization is guessed.
    """
    if bitmap_encoding not in ('linear','srgb'): raise InputError('Explicit bitmap encoding must be linear or srgb')
    values=np.asarray(bitmap,dtype=float)
    if values.ndim!=2 or values.size==0 or not np.all(np.isfinite(values)) or np.any(values<0): raise InputError('Provide finite nonnegative 2D bitmap')
    if bitmap_encoding=='srgb': values=srgb_to_linear(values)
    ny,nx=values.shape
    xs=((np.arange(nx)+.5)/nx-.5)*source.width_mm
    ys=((np.arange(ny)+.5)/ny-.5)*source.height_mm
    kernels=[];Q=[];notes=set();counts={}
    # A single panel-average radiance is supplied, so fill is not repeated.
    for j,y in enumerate(ys):
        for i,x in enumerate(xs):
            item=geometric_source_kernel(system,(x,y),image_grid,wavelength_nm,settings,
                         source_z_mm=source.z_mm,source_radiance=source.panel_radiance,fill_factor=1.,
                         sample_count=samples_per_source,seed=seed+j*nx+i,image_z_mm=image_z_mm,
                         launch_plane_z_mm=launch_plane_z_mm,launch_radius_mm=launch_radius_mm,
                         mode=mode,missing_absorption=missing_absorption)
            kernels.append(item['kernel']);Q.append(item['collected_power_density']);notes.update(item['assumptions'])
            for key,v in item['trace_status_counts'].items(): counts[key]=counts.get(key,0)+v
    result=synthesize_image(kernels,Q,np.full(values.size,source.area_mm2/values.size),bitmap_weights=values.ravel())
    result.update(units='W/mm^2' if source.is_absolute else 'relative_irradiance',
                  assumptions=sorted(notes),trace_status_counts=counts,
                  model='geometric_finite_continuous_lambertian_source',
                  source_sampling_shape=values.shape,rays_per_source=samples_per_source,
                  orientation='x columns increasing, y rows increasing',
                  source_pixel_metadata=source.pixels)
    return result


def registered_spectral_sum(results,response=None):
    """Add intensities on the same *global* grid, preserving lateral color.

    results already contain their source spectral power. Optional response has
    shape (channels,wavelengths), e.g. sensor sensitivity; it is not another copy
    of source spectrum. Raw wavelength results remain separately available.
    Cell power is the authoritative quantity: it must match the physical grid
    exactly, and irradiance is recomputed from it rather than trusting an input
    cache. Malformed data and nonfinite arithmetic fail at this public boundary.
    """
    results=tuple(results)
    if not results: raise InputError('Need spectral image results')
    powers=[];wavelengths=[]
    for r in results:
        if not isinstance(r,Mapping) or not all(k in r for k in ('grid','cell_power','wavelength_nm')):
            raise InputError('Spectral images require grid, cell_power and wavelength_nm')
        if not isinstance(r['grid'],PixelGrid): raise InputError('Spectral images require PixelGrid objects')
        r['grid'].validate()
        try:
            if np.iscomplexobj(r['cell_power']): raise InputError('Spectral cell power must be real')
            if np.iscomplexobj(r['wavelength_nm']): raise InputError('Spectral wavelength must be real')
            power=np.asarray(r['cell_power'],dtype=float)
            wavelength=float(r['wavelength_nm'])
        except (TypeError,ValueError,OverflowError) as exc:
            raise InputError('Spectral cell power and wavelength must be real numeric data') from exc
        if power.shape!=r['grid'].shape or not np.all(np.isfinite(power)) or np.any(power<0):
            raise InputError('Spectral cell power must match grid and be finite nonnegative')
        if not np.isfinite(wavelength) or wavelength<=0:
            raise InputError('Spectral wavelength must be finite positive nm')
        powers.append(power);wavelengths.append(wavelength)
    grid=results[0]['grid']
    if any(not _same_grid(grid,r['grid']) for r in results): raise InputError('Wavelength images must share physical pixel edges')
    stack=np.stack(powers)
    if response is not None:
        try:
            if np.iscomplexobj(response): raise InputError('Spectral response must be real')
            response=np.asarray(response,dtype=float)
        except (TypeError,ValueError,OverflowError) as exc:
            raise InputError('Spectral response must be real numeric data') from exc
        if response.ndim!=2 or response.shape[0]<1 or response.shape[1]!=len(results) or not np.all(np.isfinite(response)) or np.any(response<0):
            raise InputError('Response must be finite nonnegative channels x wavelengths')
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        summed=np.sum(stack,axis=0) if response is None else np.einsum('cw,wyx->cyx',response,stack)
        total=np.sum(summed,axis=(-2,-1))
        irradiance=summed/grid.cell_areas_mm2
    if not np.all(np.isfinite(summed)) or not np.all(np.isfinite(total)):
        raise InputError('Spectral image power arithmetic must remain finite; rescale input powers or response')
    if not np.all(np.isfinite(irradiance)):
        raise InputError('Spectral irradiance arithmetic must remain finite; rescale powers or pixel areas')
    return dict(grid=grid,cell_power=summed,irradiance=irradiance,
                wavelengths_nm=wavelengths,
                registration='global_physical_grid_no_per_wavelength_recentering')


def combine_complex_otfs(otfs,weights):
    """Physical registered complex OTF average; MTF is abs *after* averaging.

    Caller supplies common spatial-frequency axes and included registration
    phases. Each monochromatic OTF must have the same normalization convention.
    """
    values=np.asarray(otfs,dtype=complex);w=np.asarray(weights,dtype=float)
    if values.ndim<2 or w.shape!=(values.shape[0],) or not np.all(np.isfinite(values)) or not np.all(np.isfinite(w)) or np.any(w<0) or np.sum(w)<=0:
        raise InputError('Finite OTFs and nonnegative positive-total spectral weights required')
    result=np.einsum('i,i...->...',w/np.sum(w),values)
    return dict(otf=result,mtf=np.abs(result))


def srgb_to_linear(image):
    a=np.asarray(image,dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a<0) or np.any(a>1): raise InputError('sRGB values must lie in [0,1]')
    return np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)
