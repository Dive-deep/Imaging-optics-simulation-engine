"""Scalar pupil-autocorrelation and optional FFT diagnostics for high-NA seeds.

No core source changes. Main curves evaluate the continuous-frequency
Wiener-Khinchin identity of the original scalar Fourier-pupil profile. This is
NOT a claim that an FFT was executed for those curves. Native traces are fresh
external bundles and are never given a forged original-engine provenance seal.
"""
from pathlib import Path
import importlib.util
import sys
import math
import json
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
for path in (ROOT/'src',HERE):
    if str(path) not in sys.path:sys.path.insert(0,str(path))
from model import native_surfaces,index_normal_line,WAVELENGTHS
from native.raytrace import trace
from sequential_optics.models import Surface,OpticalSystem,Field,TraceResult,RayBundle
from sequential_optics.materials import Material
from sequential_optics.settings import EngineSettings,SamplingSettings
from sequential_optics.wave import make_boundary_state,sample_axis_otf
from sequential_optics.errors import ModelApplicabilityError

FREQUENCIES=np.linspace(0.,20.,41)
FIELDS=(0.,.5,.8,1.)
STATUS={0:'ok',1:'aperture_blocked',2:'tir',3:'numerical_failure',4:'unsupported_backward_ray'}


def build_system(design):
    """Include the separate, unpowered front STOP in the original data model."""
    mats={'air':Material('air',n_constant=1.)}
    for name in design['materials']:
        mats[name]=Material(name,wavelengths_nm=tuple(WAVELENGTHS),
                            indices=tuple(index_normal_line(name,WAVELENGTHS)))
    radius=float(design['stop_radius_mm'])
    surfaces=[Surface(-float(design['stop_distance_mm']),0.,radius,medium_after='air',name='STOP')]
    for j,row in enumerate(native_surfaces(design)):
        name=design['materials'][j//2]
        surfaces.append(Surface(float(row[0]),float(row[1]),float(row[8]),conic=float(row[2]),
            coefficients=(0.,*map(float,row[3:8])),medium_after=name if j%2==0 else 'air',
            surface_type='even_asphere' if name in ('PMMA','PC') else 'standard',name=f'S{j+1}'))
    return OpticalSystem(tuple(surfaces),mats,object_z=None,image_z=float(sum(design['axial_distances_mm'])),
                         stop_index=0,stop_radius=radius,name=str(design.get('candidate_id',design.get('seed_id','candidate'))))


def settings(rings=16,arms=48,grid_size=128,padding=4):
    return EngineSettings(wavelengths_nm={'B':486.,'G':588.,'R':656.},primary_wavelength='G',
        sampling=SamplingSettings(pupil_rings=rings,pupil_arms=arms,wave_grid_size=grid_size,
            fft_padding=padding,max_frequency_lp_mm=20.,frequency_samples=41))


def gauss_pupil(rings=16,arms=48):
    nodes,w=np.polynomial.legendre.leggauss(rings)
    points=[(0.,0.)];weights=[0.]
    for i,(node,weight) in enumerate(zip(nodes,w)):
        r=np.sqrt((node+1)/2);a=2*np.pi*(np.arange(arms)+.5*(i%2))/arms
        points.extend(np.c_[r*np.cos(a),r*np.sin(a)].tolist());weights.extend([weight/(2*arms)]*arms)
    a=2*np.pi*np.arange(arms)/arms
    points.extend(np.c_[np.cos(a),np.sin(a)].tolist());weights.extend([0.]*arms)
    return np.asarray(points),np.asarray(weights)


def native_bundle(design,field_fraction,wavelength_nm,*,rings=16,arms=48):
    xy,weights=gauss_pupil(rings,arms)
    theta=np.deg2rad(12.)*float(field_fraction)
    direction=np.array([0.,np.sin(theta),np.cos(theta)])
    origins=np.c_[xy*float(design['stop_radius_mm']),np.full(len(xy),-float(design['stop_distance_mm']))]
    directions=np.broadcast_to(direction,origins.shape).copy()
    indices=np.ones(7)
    for j,name in enumerate(design['materials']):indices[2*j+1]=index_normal_line(name,wavelength_nm)
    raw=trace(native_surfaces(design),origins,directions,indices)
    initial_opl=origins@direction
    traces=[]
    for i,code in enumerate(raw['status']):
        good=np.all(np.isfinite(raw['history'][i]),axis=1)
        pos=np.vstack([origins[i],raw['history'][i,good]])
        dirs=np.vstack([direction,raw['direction_history'][i,good]])
        opl=np.r_[initial_opl[i],raw['opl_history'][i,good]+initial_opl[i]]
        traces.append(TraceResult(pos,dirs,opl,np.ones(len(pos)),STATUS[int(code)],
                                  message='' if code==0 else f'native failed surface {raw["failed_surface"][i]}'))
    bundle=RayBundle(tuple(traces),xy,weights,Field(0.,theta),float(wavelength_nm),traces[0])
    return bundle,raw


class PiecewisePupil:
    """Piecewise-unwrapped-OPD pupil retaining transmitted STOP triangles."""
    def __init__(self,vertices,triangles,opd,amplitude=None):
        from matplotlib.tri import Triangulation
        vertices=np.asarray(vertices,float);triangles=np.asarray(triangles,int);opd=np.asarray(opd,float)
        used=np.unique(triangles)
        remap=np.full(len(vertices),-1,int);remap[used]=np.arange(len(used))
        self.vertices=vertices[used];self.triangles=remap[triangles];self.opd=opd[used]
        self.amplitude=np.ones(len(used)) if amplitude is None else np.asarray(amplitude,float)[used]
        v=self.vertices[self.triangles];a=v[:,1]-v[:,0];b=v[:,2]-v[:,0]
        det=a[:,0]*b[:,1]-a[:,1]*b[:,0]
        if np.any(~np.isfinite(det)) or np.any(det==0):raise ModelApplicabilityError('Invalid pupil triangle')
        # Matplotlib's trapezoid finder requires CCW connectivity.
        self.triangles[det<0]=self.triangles[det<0][:,[0,2,1]]
        v=self.vertices[self.triangles];a=v[:,1]-v[:,0];b=v[:,2]-v[:,0]
        det=a[:,0]*b[:,1]-a[:,1]*b[:,0]
        self.area=float(det.sum()/2)
        self.v0=v[:,0]
        self.values=np.c_[self.opd,self.amplitude][self.triangles]
        d1=self.values[:,1]-self.values[:,0];d2=self.values[:,2]-self.values[:,0]
        self.gradx=(d1*b[:,1,None]-d2*a[:,1,None])/det[:,None]
        self.grady=(d2*a[:,0,None]-d1*b[:,0,None])/det[:,None]
        self.triangulation=Triangulation(self.vertices[:,0],self.vertices[:,1],self.triangles)
        try:self.finder=self.triangulation.get_trifinder()
        except RuntimeError as exc:raise ModelApplicabilityError('Overlapping/nonregular angular-pupil topology') from exc
        self.minimum=self.vertices.min(axis=0);self.maximum=self.vertices.max(axis=0)

    def evaluate(self,x,y):
        x=np.asarray(x);y=np.asarray(y);owner=self.finder(x,y);inside=owner>=0
        q=np.maximum(owner,0)
        values=self.values[q,0]+self.gradx[q]*(x-self.v0[q,0])[...,None]+self.grady[q]*(y-self.v0[q,1])[...,None]
        return values[...,0],values[...,1],inside

    def otf(self,frequencies,wavelength_mm,*,axis='x',grid_size=128):
        """Unwrapped phase-difference midpoint quadrature, arbitrary frequencies."""
        frequencies=np.asarray(frequencies,float)
        span=self.maximum-self.minimum
        step=span/float(grid_size)
        x=self.minimum[0]+(np.arange(grid_size)+.5)*step[0]
        y=self.minimum[1]+(np.arange(grid_size)+.5)*step[1]
        xx,yy=np.meshgrid(x,y);x=xx.ravel();y=yy.ravel()
        _,amp,mask=self.evaluate(x,y)
        denominator=float(np.sum(np.where(mask,amp*amp,0.)))
        if denominator<=0:raise ModelApplicabilityError('No pupil area in quadrature')
        answer=[]
        for freq in frequencies:
            if freq==0:answer.append(1.+0j);continue
            delta=wavelength_mm*float(freq)/2
            if axis=='x':
                wm,am,mm=self.evaluate(x-delta,y);wp,ap,mp=self.evaluate(x+delta,y)
            elif axis=='y':
                wm,am,mm=self.evaluate(x,y-delta);wp,ap,mp=self.evaluate(x,y+delta)
            else:raise ValueError('axis must be x or y')
            good=mm&mp
            # Negative OTF phase: a pupil tilt +a*u yields image displacement +a.
            answer.append(np.sum(am[good]*ap[good]*np.exp(2j*np.pi*(wm[good]-wp[good])/wavelength_mm))/denominator)
        return np.asarray(answer),{'quadrature_grid':grid_size,'sampled_area':float(mask.sum()*np.prod(step)),
            'triangle_area':self.area,'relative_area_error':float(abs(mask.sum()*np.prod(step)-self.area)/self.area),
            'normalization':'same-grid integral |P|²; DC exactly one','phase_interpolation':'unwrapped OPD before phase difference',
            'interpolation_chart':('physical_UV_with_exact_inverse_direction_cosine_query' if isinstance(self,InverseMappedBoundaryPupil) else 'direction_cosine_linear'),
            'area_reference':('7-point-per-UV-triangle projected-area quadrature estimate' if isinstance(self,InverseMappedBoundaryPupil) else 'exact straight-triangle area'),
            'max_mtf_above_one':float(max(0.,np.max(np.abs(answer))-1))}



class InverseMappedBoundaryPupil(PiecewisePupil):
    """Exact angular-to-boundary support query; UV-linear unwrapped OPD.

    Integration remains dchi*deta with uniform Fourier-pupil B=1. No UV
    Jacobian is multiplied into amplitude. The Jacobian is needed only for the
    independent projected-area estimate below, whose integral is in dU*dV.
    """
    def __init__(self,uv,triangles,opd,distance):
        self.boundary=PiecewisePupil(uv,triangles,opd)
        self.distance=float(distance)
        if self.distance<=0:raise ModelApplicabilityError('Nonpositive reference distance')
        self.opd=self.boundary.opd;self.triangles=self.boundary.triangles
        # A conservative bounding rectangle, including extrema along curved
        # mapped edges. Merely transforming vertices can omit such extrema.
        lo=self.boundary.minimum;hi=self.boundary.maximum
        minabs=np.where((lo<=0)&(hi>=0),0.,np.minimum(abs(lo),abs(hi)))
        maxabs=np.maximum(abs(lo),abs(hi))
        lower_other=np.where(lo<=0,minabs[::-1],maxabs[::-1])
        upper_other=np.where(hi>=0,minabs[::-1],maxabs[::-1])
        self.minimum=lo/np.sqrt(self.distance**2+lo**2+lower_other**2)
        self.maximum=hi/np.sqrt(self.distance**2+hi**2+upper_other**2)
        # Degree-five triangle quadrature for projected angular area. This is
        # an independent diagnostic estimate, not the OTF normalization.
        bary=np.array([[1/3,1/3,1/3],
            [.059715871789770,.470142064105115,.470142064105115],
            [.470142064105115,.059715871789770,.470142064105115],
            [.470142064105115,.470142064105115,.059715871789770],
            [.797426985353087,.101286507323456,.101286507323456],
            [.101286507323456,.797426985353087,.101286507323456],
            [.101286507323456,.101286507323456,.797426985353087]])
        weights=np.array([.225]+[.132394152788506]*3+[.125939180544827]*3)
        vertices=self.boundary.vertices[self.triangles]
        a=vertices[:,1]-vertices[:,0];b=vertices[:,2]-vertices[:,0]
        area=.5*(a[:,0]*b[:,1]-a[:,1]*b[:,0])
        nodes=np.einsum('qj,tjk->tqk',bary,vertices)
        jac=self.distance**2/(self.distance**2+np.sum(nodes*nodes,axis=-1))**2
        self.area=float(np.sum(area*np.sum(jac*weights,axis=-1)))

    def evaluate(self,x,y):
        x=np.asarray(x);y=np.asarray(y);g2=1-x*x-y*y
        scale=self.distance/np.sqrt(np.maximum(g2,1e-30))
        w,a,m=self.boundary.evaluate(x*scale,y*scale)
        return w,a,m&(g2>0)



def _down(x):return np.nextafter(x,-np.inf)
def _up(x):return np.nextafter(x,np.inf)
def _mul_interval(al,ah,bl,bh):
    values=(al*bl,al*bh,ah*bl,ah*bh)
    low=np.minimum(np.minimum(values[0],values[1]),np.minimum(values[2],values[3]))
    high=np.maximum(np.maximum(values[0],values[1]),np.maximum(values[2],values[3]))
    return _down(low),_up(high)


def conservative_surface_projection_upper(surface,axis,intervals=32768):
    """Outward-rounded radial interval envelope over the entire clear surface.

    f(r)=axis_z*(z+sag(r))+hypot(axis_x,axis_y)*r is the exact azimuthal
    maximum at radius r. A union of radial intervals bounds its global max.
    Signed aspheric terms use their endpoint power ranges, retaining local
    cancellations rather than summing each term's whole-disk maximum.
    """
    axis=np.asarray(axis,float)
    if axis.shape!=(3,) or not np.all(np.isfinite(axis)) or axis[2]<=0:
        raise ModelApplicabilityError('Forward finite chief axis required')
    radius=float(surface.clear_radius);c=float(surface.curvature);k=float(surface.conic)
    edges=np.linspace(0.,radius,int(intervals)+1);lo=edges[:-1];hi=edges[1:]
    qlo=np.maximum(0.,_down(lo*lo));qhi=_up(hi*hi)
    # Conic sag is monotone with sign(c). Evaluate its maximizing endpoint
    # with outward basic arithmetic, not a sampled maximization of f itself.
    rr=hi if c>=0 else lo
    rl=np.maximum(0.,_down(rr*rr));rh=_up(rr*rr)
    kl,kh=_down(k+1.),_up(k+1.)
    cl,ch=max(0.,float(_down(c*c))),float(_up(c*c))
    tl,th=_mul_interval(kl,kh,cl,ch)
    whole_domain_lower=float(_down(1-_up(max(0.,float(th))*float(_up(radius*radius)))))
    if whole_domain_lower<=0:raise ModelApplicabilityError('Whole clear aperture conic domain cannot be bounded positive')
    tl,th=_mul_interval(tl,th,rl,rh)
    dl,dh=_down(1-th),_up(1-tl)
    # An uncertain hemisphere boundary is rejected rather than clipping it.
    if np.any(dl<=0):raise ModelApplicabilityError('Conic domain cannot be bounded strictly positive')
    sl,sh=_down(np.sqrt(dl)),_up(np.sqrt(dh))
    denl,denh=_down(1+sl),_up(1+sh)
    nl,nh=_mul_interval(c,c,rl,rh)
    invl,invh=_down(1/denh),_up(1/denl)
    saglo,saghi=_mul_interval(nl,nh,invl,invh)
    _,upper=_mul_interval(axis[2],axis[2],saglo,saghi)
    constant=float(_up(axis[2]*surface.z));upper=_up(upper+constant)
    magnitude=np.abs(upper)+abs(constant)
    powerlo=np.ones_like(lo);powerhi=np.ones_like(hi)
    for coefficient in surface.coefficients:
        powerlo=np.maximum(0.,_down(powerlo*qlo));powerhi=_up(powerhi*qhi)
        al,ah=_mul_interval(axis[2],axis[2],float(coefficient),float(coefficient))
        termlo,termhi=_mul_interval(al,ah,powerlo,powerhi)
        upper=_up(upper+termhi);magnitude+=np.maximum(abs(termlo),abs(termhi))
    transverse=float(_up(np.sqrt(_up(_up(axis[0]*axis[0])+_up(axis[1]*axis[1])))))
    radial=_up(transverse*hi);upper=_up(upper+radial);magnitude+=abs(radial)
    if not np.all(np.isfinite(upper)):raise ModelApplicabilityError('Overflow in conservative surface projection bound')
    # Explicit engineering slack beyond outward elementary-operation bounds.
    margin=1e-12+64*np.finfo(float).eps*float(np.max(1+magnitude))
    return float(_up(float(np.max(upper))+margin))


def _private_profile():
    global _PRIVATE_WAVE
    if _PRIVATE_WAVE is None:
        spec=importlib.util.spec_from_file_location('sequential_optics._candidate_wave_profile',ROOT/'src/sequential_optics/wave.py')
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        module._mapped_interpolation=fast_mapped_interpolation
        module._full_surface_projection_upper=conservative_surface_projection_upper
        _PRIVATE_WAVE=module
    return _PRIVATE_WAVE


def pupil_from_bundle(system,bundle,engine_settings):
    state=_private_profile().make_boundary_state(system,bundle,engine_settings,allow_unverified_bundle=True)
    q0=state.chief_position+(system.image_z-state.chief_position[2])/state.chief_direction[2]*state.chief_direction
    L=float(np.dot(q0-state.origin,state.chief_direction))
    uv=state.uv_mm;r=np.sqrt(L*L+np.sum(uv*uv,axis=1))
    vertices=uv/r[:,None]
    opd=state.opl_mm-state.chief_opl_mm+state.output_index*np.sum(uv*uv,axis=1)/(r+L)
    pupil=InverseMappedBoundaryPupil(uv,state.triangles,opd,L)
    return pupil,state,q0


def geometric_otf(bundle,image_z,frequencies=FREQUENCIES):
    weights=bundle.weights.copy();valid=bundle.success_mask
    weights*=valid
    if weights.sum()<=0:return None
    points=np.zeros((len(weights),2))
    for i,tr in enumerate(bundle.traces):
        if valid[i]:
            t=(image_z-tr.last_position[2])/tr.last_direction[2]
            points[i]=tr.last_position[:2]+t*tr.last_direction[:2]
    weights/=weights.sum()
    x=np.sum(weights[:,None]*np.exp(-2j*np.pi*points[:,0,None]*frequencies),axis=0)
    y=np.sum(weights[:,None]*np.exp(-2j*np.pi*points[:,1,None]*frequencies),axis=0)
    center=np.sum(weights[:,None]*points,axis=0)
    rms=np.sqrt(np.sum(weights*np.sum((points-center)**2,axis=1)))
    return {'otf_s':x,'otf_t':y,'rms_radius_mm':float(rms),'centroid_mm':center,'method':'geometric ray-intercept characteristic function; no diffraction'}


def _complex_json(values):return {'real':np.asarray(values).real.tolist(),'imag':np.asarray(values).imag.tolist()}
def _complex_read(value):return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def evaluate_one(design,field_fraction,wavelength_nm,*,rings=16,arms=48,quadrature=(128,256),frequencies=FREQUENCIES,fft_diagnostic=False,fft_grid=1024):
    start=time.perf_counter();system=build_system(design);s=settings(rings,arms)
    bundle,raw=native_bundle(design,field_fraction,wavelength_nm,rings=rings,arms=arms)
    good=bundle.success_mask
    result={'field_fraction':float(field_fraction),'field_angle_deg':12*float(field_fraction),'wavelength_nm':float(wavelength_nm),
      'frequencies_lp_mm':np.asarray(frequencies).tolist(),'traced_ray_count':len(good),'successful_ray_count':int(good.sum()),
      'status_counts':{STATUS[int(k)]:int(v) for k,v in zip(*np.unique(raw['status'],return_counts=True))},
      'geometric_transmission':float(np.sum(bundle.weights*good)),
      'method':'scalar_pupil_autocorrelation_MTF','fft_executed_for_primary_curves':False,
      'scope':'Same scalar Fourier-pupil model as v0.2.1; continuous requested-frequency autocorrelation. High-NA vector/polarization equivalence unverified.',
      'trace_provenance_status':'unverified_external_bundle_fresh_native_trace',
      'illumination':'uniform Fourier-pupil amplitude; no Fresnel/absorption in lens-only MTF',
      'axis_convention':'Y-field: T frequency along global sensor Y; S along global sensor X',
      'sensor_plane_conversion':'local chief-normal projection; T frequency / chief_z_cosine, not exact tilted-plane vector diffraction',
      'interpolation_chart':'unwrapped reference-sphere OPD linear on physical UV triangles; exact angular-to-UV inverse query',
      'interpolation_change':'Original v0.2.1 interpolates linearly on straight direction-cosine triangles; this conditional variant preserves curved support through inverse UV lookup.',
      'integration_measure':'uniform dchi*deta; original B=1 Fourier-pupil amplitude, no duplicate UV Jacobian',
      'boundary_projection_bound':'32768 radial interval upper envelopes, outward nextafter elementary arithmetic plus1e-12mm+64eps magnitude margin'}
    geom=geometric_otf(bundle,system.image_z,np.asarray(frequencies))
    if geom is not None:
        result['geometric_diagnostic']={**{k:v for k,v in geom.items() if k not in ('otf_s','otf_t','centroid_mm')},
            'centroid_mm':geom['centroid_mm'].tolist(),'mtf_s':abs(geom['otf_s']).tolist(),'mtf_t':abs(geom['otf_t']).tolist(),
            'otf_s':_complex_json(geom['otf_s']),'otf_t':_complex_json(geom['otf_t'])}
    try:
        pupil,state,q0=pupil_from_bundle(system,bundle,s)
        cosine=float(state.chief_direction[2]);frequencies=np.asarray(frequencies,float)
        if cosine<=0:raise ModelApplicabilityError('Backward chief')
        if abs(state.chief_direction[0])>1e-10:raise ModelApplicabilityError('Only Y-field sensor projection implemented')
        histories=[];prev=None
        for grid in quadrature:
            os,ds=pupil.otf(frequencies,wavelength_nm*1e-6,axis='x',grid_size=grid)
            ot,dt=pupil.otf(frequencies/cosine,wavelength_nm*1e-6,axis='y',grid_size=grid)
            os*=np.exp(-2j*np.pi*frequencies*q0[0]);ot*=np.exp(-2j*np.pi*frequencies*q0[1])
            entry={'grid':grid,'sagittal':ds,'tangential':dt}
            if prev is not None:
                delta=np.r_[abs(abs(os)-abs(prev[0])),abs(abs(ot)-abs(prev[1]))]
                denom=np.r_[abs(os),abs(ot)]
                entry['max_mtf_absolute_change']=float(delta.max())
                nz=denom!=0
                entry['max_mtf_relative_change_percent_nonzero']=float(100*np.max(delta[nz]/denom[nz])) if nz.any() else 0.
                entry['zero_reference_changes']=int(np.count_nonzero((~nz)&(delta!=0)))
            histories.append(entry);prev=(os,ot)
        last=histories[-1]
        # This is a declared engineering quadrature diagnostic, not a rigorous error bound.
        convergence=(len(histories)>1 and last['max_mtf_relative_change_percent_nonzero']<=.5 and last['zero_reference_changes']==0)
        absolute_stability=(len(histories)>1 and last['max_mtf_absolute_change']<=.005)
        ray_span=float(np.max(np.ptp(pupil.opd[pupil.triangles],axis=1))*2*np.pi/(wavelength_nm*1e-6))
        result.update(status='conditional' if convergence else 'underresolved_quadrature',
          usable_for_final_accuracy_claim=False,quadrature_convergence_strict_relative_005=convergence,
          quadrature_absolute_stability_005=absolute_stability,
          ray_mesh_convergence_checked=False,quadrature_history=histories,
          original_fft_phase_triangle_span_rad=ray_span,
          original_fft_phase_span_below_pi=ray_span<np.pi,
          image_reference_point_mm=q0.tolist(),chief_direction=state.chief_direction.tolist(),
          mtf_s=abs(os).tolist(),mtf_t=abs(ot).tolist(),otf_s=_complex_json(os),otf_t=_complex_json(ot),
          boundary_distance_mm=float(np.dot(q0-state.origin,state.chief_direction)),
          boundary_plane_bound_mm=float(state.plane_bound),boundary_origin_mm=state.origin.tolist(),
          pupil_triangle_count=len(state.triangles),
          necessary_note='Autocorrelation avoids wrapping/undersampled phase-field FFT, but unwrapped OPD interpolation still requires independent traced-mesh refinement.')
        if fft_diagnostic:
            try:
                result['actual_fft_diagnostic']=fft_curves(system,bundle,s,frequencies,grid_size=fft_grid)
            except (ValueError,RuntimeError,ArithmeticError,ModelApplicabilityError) as exc:
                result['actual_fft_diagnostic']={'method':'actual_scalar_FFT','status':'not_applicable',
                    'grid_size':fft_grid,'reason':f'{type(exc).__name__}: {exc}',
                    'note':'Original straight angular-triangle profile failure does not erase separate inverse-UV primary results.'}
    except (ValueError,RuntimeError,ArithmeticError,ModelApplicabilityError) as exc:
        result.update(status='not_applicable',reason=f'{type(exc).__name__}: {exc}',mtf_s=None,mtf_t=None)
    result['seconds']=time.perf_counter()-start
    return result


def fast_mapped_interpolation(vertices,triangles,values,x,y,*,overlap_tolerance=1e-10):
    """Original triangle formula/ownership, local integer grid bounding boxes."""
    x=np.asarray(x);y=np.asarray(y);dx=x[1]-x[0];dy=y[1]-y[0]
    out=np.zeros((len(y),len(x),values.shape[1]));mask=np.zeros((len(y),len(x)),bool);owner=np.full(mask.shape,-1,int)
    for ordinal,tri in enumerate(triangles):
        v=vertices[tri];a=v[1]-v[0];b=v[2]-v[0];det=a[0]*b[1]-a[1]*b[0]
        if det==0:raise ModelApplicabilityError('Degenerate mapped triangle')
        lo=np.floor((v.min(axis=0)-overlap_tolerance-np.array([x[0],y[0]]))/[dx,dy]).astype(int)
        hi=np.ceil((v.max(axis=0)+overlap_tolerance-np.array([x[0],y[0]]))/[dx,dy]).astype(int)
        ix=np.arange(max(0,lo[0]),min(len(x)-1,hi[0])+1);iy=np.arange(max(0,lo[1]),min(len(y)-1,hi[1])+1)
        if not len(ix) or not len(iy):continue
        xx,yy=np.meshgrid(x[ix],y[iy]);rx=xx-v[0,0];ry=yy-v[0,1]
        w1=(rx*b[1]-ry*b[0])/det;w2=(a[0]*ry-a[1]*rx)/det;w0=1-w1-w2
        inside=(w0>=-overlap_tolerance)&(w1>=-overlap_tolerance)&(w2>=-overlap_tolerance)
        yyids,xxids=np.where(inside);gy=iy[yyids];gx=ix[xxids]
        strict=(w0[inside]>overlap_tolerance)&(w1[inside]>overlap_tolerance)&(w2[inside]>overlap_tolerance)
        if np.any(mask[gy,gx]&strict):raise ModelApplicabilityError('Overlapping pupil branches')
        out[gy,gx]=w0[inside,None]*values[tri[0]]+w1[inside,None]*values[tri[1]]+w2[inside,None]*values[tri[2]]
        mask[gy,gx]=True;owner[gy,gx]=ordinal
    return out,mask,owner


_PRIVATE_WAVE=None

def fft_curves(system,bundle,engine_settings,frequencies,*,grid_size=1024):
    profile=_private_profile()
    raw=profile.analyze_wave(system,bundle,engine_settings,grid_size=grid_size,padding=2,
                                  frequencies=frequencies,allow_unverified_bundle=True)
    cosine=float(bundle.chief.last_direction[2])
    os=sample_axis_otf(raw,frequencies,'x');ot=sample_axis_otf(raw,np.asarray(frequencies)/cosine,'y')
    return {'method':'actual_scalar_FFT','grid_size':grid_size,'padding':2,'status':raw['status'],
        'mtf_s':abs(os).tolist(),'mtf_t':abs(ot).tolist(),'sampling_diagnostics':raw['sampling_diagnostics'],
        'native_frequency_step_lp_mm':{'x':float(raw['frequency_x'][1]-raw['frequency_x'][0]),'y':float(raw['frequency_y'][1]-raw['frequency_y'][0])},
        'sensor_plane_conversion':'same local chief-normal projection',
        'local_expansion_diagnostics':raw.get('local_expansion_diagnostics'),
        'warning':'If underresolved, numerical FFT curves are diagnostics only and must not rank candidates.'}


def evaluate_design(design,*,fields=FIELDS,wavelengths=WAVELENGTHS,rings=16,arms=48,quadrature=(128,256),fft_diagnostics=False):
    records=[]
    for field in fields:
        for wavelength in wavelengths:
            records.append(evaluate_one(design,field,wavelength,rings=rings,arms=arms,quadrature=quadrature,
                fft_diagnostic=fft_diagnostics and wavelength==588 and field in (0.,.8)))
    polychromatic=[]
    for field in fields:
        group=[r for r in records if r['field_fraction']==field]
        if not all(r.get('otf_t') for r in group):continue
        os=sum(_complex_read(r['otf_s']) for r in group)/len(group)
        ot=sum(_complex_read(r['otf_t']) for r in group)/len(group)
        polychromatic.append({'field_fraction':float(field),'frequencies_lp_mm':FREQUENCIES.tolist(),
            'mtf_s':abs(os).tolist(),'mtf_t':abs(ot).tolist(),'otf_s':_complex_json(os),'otf_t':_complex_json(ot),
            'spectral_weighting':'equal normalized image energy in each supplied wavelength; detector response unspecified',
            'registration':'common physical sensor origin; wavelength chief displacements retained in complex OTF',
            'status':'conditional' if all(r['status']=='conditional' for r in group) else 'underresolved_quadrature'})
    return {'candidate_id':design.get('candidate_id',design.get('seed_id')),'monochromatic':records,'polychromatic':polychromatic,
        'wavelengths_nm':list(map(float,wavelengths)),'fields':list(map(float,fields)),
        'method_note':'Primary curves are continuous-frequency scalar pupil autocorrelation, not actual FFT outputs. Optional FFT diagnostics are separate.',
        'na_model_limit':'Scalar and local sensor projection only; vector/polarization high-NA corrections not certified.',
        'seconds':sum(x['seconds'] for x in records)}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('design');p.add_argument('--output',required=True)
    p.add_argument('--fft',action='store_true');p.add_argument('--rings',type=int,default=16);p.add_argument('--arms',type=int,default=48)
    args=p.parse_args();d=json.loads(Path(args.design).read_text())
    out=evaluate_design(d,rings=args.rings,arms=args.arms,fft_diagnostics=args.fft)
    Path(args.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
