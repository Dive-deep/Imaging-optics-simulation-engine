"""Independent analytic/scalar controls; never changes the original engine."""
from pathlib import Path
import json,time
import numpy as np
from scipy.spatial import Delaunay
from wave_evaluation import PiecewisePupil,InverseMappedBoundaryPupil,fast_mapped_interpolation
from sequential_optics.wave import circular_mtf,fft_psf,sample_axis_otf,_mapped_interpolation

start=time.perf_counter();checks=[]
def check(name,ok,**details):
 checks.append({'name':name,'passed':bool(ok),**details})
 if not ok:print('FAILED',name,details,flush=True)

# Dense polygon triangulates the exact declared aperture; use many rings so
# piecewise quadratic defocus is independently close to the analytic surface.
a=.7;lam=.000588
points=[(0.,0.)]
for r in np.linspace(a/20,a,20):
 theta=2*np.pi*np.arange(256)/256
 points.extend(np.c_[r*np.cos(theta),r*np.sin(theta)].tolist())
v=np.asarray(points);tri=Delaunay(v).simplices
p=PiecewisePupil(v,tri,np.zeros(len(v)))
f=np.array([0.,.1,.25,.5,.75,.95,1.,1.1])*2*a/lam
z,diag=p.otf(f,lam,grid_size=384)
expected=circular_mtf(f,2*a/lam)
check('clear circular pupil vs analytic diffraction MTF',max(abs(abs(z)-expected))<.0015,
 max_absolute_error=float(max(abs(abs(z)-expected))),frequencies_lp_mm=f.tolist(),actual=abs(z).tolist(),expected=expected.tolist())
check('clear pupil imaginary OTF zero',max(abs(z.imag))<1e-15)
check('DC normalized and beyond cutoff zero',z[0]==1 and abs(z[-1])==0)

# Pure pupil tilt: PSF shift must change OTF phase, not its magnitude.
shift=.0217
pt=PiecewisePupil(v,tri,shift*v[:,0])
zt,_=pt.otf(f,lam,grid_size=192)
z0,_=p.otf(f,lam,grid_size=192)
zt_exact=z0*np.exp(-2j*np.pi*f*shift)
check('pupil tilt gives correct complex OTF translation sign',max(abs(zt-zt_exact))<2e-11,max_complex_error=float(max(abs(zt-zt_exact))))

# Two shifted monochromatic images combine as complex OTF, not |OTF|.
fsmall=np.linspace(0,20,41);d=.025
poly=(np.ones_like(fsmall)+np.exp(-2j*np.pi*fsmall*d))/2
check('polychromatic displacement preserved',max(abs(abs(poly)-abs(np.cos(np.pi*fsmall*d))))<1e-14 and abs(poly[-1])<1e-14,
 mtf_at_20=float(abs(poly[-1])),incorrect_mean_magnitude=1.)

# Compare smooth phase with actual FFT at its own native frequency bins. No
# phase alias is present in this control; both implement the scalar identity.
N=256;extent=a*1.05;du=2*extent/N;u=(np.arange(N)-N//2)*du;xx,yy=np.meshgrid(u,u)
b=.0006
pd=PiecewisePupil(v,tri,b*np.sum(v*v,axis=1))
w,amp,mask=pd.evaluate(xx,yy)
field=np.where(mask,np.exp(2j*np.pi*w/lam),0j)
ff=fft_psf(field,(du,du),lam,padding=2)
step=ff['frequency_x'][1]-ff['frequency_x'][0]
fd=np.arange(0,81,8)*step
zd,_=pd.otf(fd,lam,grid_size=384)
zf=sample_axis_otf(ff,fd,'x')
check('defocus continuous autocorrelation vs actual FFT',max(abs(zd-zf))<.005,
 max_complex_error=float(max(abs(zd-zf))),fft_frequency_step_lp_mm=float(step))

# Preserve original affine triangle interpolation, support and ownership.
rng=np.random.default_rng(82);v2=np.r_[np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]]),rng.uniform(-.9,.9,(25,2))]
t2=Delaunay(v2).simplices;values=rng.normal(size=(len(v2),2));x=np.linspace(-1.2,1.2,64);y=np.linspace(-1.1,1.1,64)
a0,m0,o0=_mapped_interpolation(v2,t2,values,x,y);a1,m1,o1=fast_mapped_interpolation(v2,t2,values,x,y)
check('fast FFT interpolation support and triangle owner unchanged',np.array_equal(m0,m1) and np.array_equal(o0,o1))
check('fast FFT interpolation values unchanged',max(abs(a0.ravel()-a1.ravel()))<1e-12,max_error=float(max(abs(a0.ravel()-a1.ravel()))))
# Inverse chart: angular measure remains unchanged; compare projected UV disk.
R=2.;L=3.;angular_radius=R/np.sqrt(L*L+R*R)
vu=v/a*R;tu=tri
pi=InverseMappedBoundaryPupil(vu,tu,np.full(len(vu),.123),L)
fi=np.array([0.,.1,.25,.5,.75,1.])*2*angular_radius/lam
zi,_=pi.otf(fi,lam,grid_size=384)
ei=circular_mtf(fi,2*angular_radius/lam)
check('inverse UV constant-OPD circular pupil analytic MTF',max(abs(abs(zi)-ei))<.0015,
 max_error=float(max(abs(abs(zi)-ei))))
check('inverse UV projected angular area and Jacobian',abs(pi.area-np.pi*angular_radius**2)/(np.pi*angular_radius**2)<.001,
 relative_error=float(abs(pi.area-np.pi*angular_radius**2)/(np.pi*angular_radius**2)))

# A chi-linear pupil tilt is no longer exactly linear in UV. Verify convergence
# under independent traced-mesh-like UV refinement, rather than claiming exact
# equivalence of the two interpolation charts.
errors=[]
for nr,na in ((6,48),(12,96),(24,192)):
 uvp=[(0.,0.)]
 for rr in np.linspace(R/nr,R,nr):
  aa=2*np.pi*np.arange(na)/na
  uvp.extend(np.c_[rr*np.cos(aa),rr*np.sin(aa)].tolist())
 uvp=np.asarray(uvp);tt=Delaunay(uvp).simplices
 chi=uvp[:,0]/np.sqrt(L*L+np.sum(uvp*uvp,axis=1));shift=.0217
 base=InverseMappedBoundaryPupil(uvp,tt,np.zeros(len(uvp)),L)
 tilted=InverseMappedBoundaryPupil(uvp,tt,shift*chi,L)
 ft=np.array([0.,5.,10.,15.,20.])
 zb,_=base.otf(ft,lam,grid_size=192);zz,_=tilted.otf(ft,lam,grid_size=192)
 errors.append(float(max(abs(zz-zb*np.exp(-2j*np.pi*ft*shift)))))
check('inverse UV tilt phase converges with mesh refinement',errors[2]<errors[1]<errors[0] and errors[2]<.001,
 max_complex_errors_by_mesh=errors,meshes=[[6,48],[12,96],[24,192]])

report={'passed':all(c['passed'] for c in checks),'checks':checks,'seconds':time.perf_counter()-start,
 'scope':'Analytic controls of scalar Fourier identity, translation phase, RGB complex sum, and unchanged interpolation. Not high-NA vector validation or candidate convergence.'}
out=Path(__file__).parent/'results/wave_controls_v2.json';out.parent.mkdir(exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
if not report['passed']:raise SystemExit(1)
