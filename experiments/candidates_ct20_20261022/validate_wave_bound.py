from pathlib import Path
import json,time
import numpy as np
from wave_evaluation import conservative_surface_projection_upper,Surface
from sequential_optics.wave import _full_surface_projection_upper
from sequential_optics.geometry import sag
start=time.perf_counter();rng=np.random.default_rng(371);checks=[]
for j in range(80):
 radius=rng.uniform(2,30);k=rng.uniform(-6,0);c=rng.uniform(-.8,.8)/radius
 coeff=np.r_[0.,rng.uniform(-15,15,5)/radius**np.array([4,6,8,10,12])]
 surface=Surface(rng.uniform(-5,60),c,radius,conic=k,coefficients=tuple(coeff),surface_type='even_asphere')
 angle=rng.uniform(-.7,.7);axis=np.array([0.,np.sin(angle),np.cos(angle)])
 r=np.linspace(0,radius,32769)
 actual=axis[2]*(surface.z+sag(surface,r,np.zeros_like(r)))+abs(axis[1])*r
 bound=conservative_surface_projection_upper(surface,axis)
 old=_full_surface_projection_upper(surface,axis)
 checks.append({'case':j,'passed':bool(bound>=max(actual)),'upper_bound_mm':bound,
    'dense_maximum_mm':float(max(actual)),'bound_excess_mm':float(bound-max(actual)),
    'previous_coarse_bound_mm':float(old),'relative_old_bound_gap':float(old-max(actual))})
# A cancellation pair, whose positive-term whole-disk bound is very loose.
s=Surface(30.,.01,25.,conic=-1.,coefficients=(0.,20/25**4,-20/25**6),surface_type='even_asphere')
b=conservative_surface_projection_upper(s,np.array([0.,0.,1.]));old=_full_surface_projection_upper(s,np.array([0.,0.,1.]))
r=np.linspace(0,25,100001);actual=30+sag(s,r,np.zeros_like(r))
checks.append({'case':'cancelling_high_order_terms','passed':bool(b>=max(actual) and b<old-10),
 'upper_bound_mm':b,'dense_maximum_mm':float(max(actual)),'previous_coarse_bound_mm':float(old)})
# Undefined rim must not be hidden by using the r_lo endpoint of a negative c.
try:
 conservative_surface_projection_upper(Surface(0.,-.1,11.),np.array([0.,0.,1.]))
 rejected=False
except Exception:rejected=True
checks.append({'case':'negative_curvature_undefined_rim_rejected','passed':rejected})
report={'passed':all(c['passed'] for c in checks),'checks':checks,'check_count':len(checks),
 'seconds':time.perf_counter()-start,'scope':'Analytic signed-term interval envelope is conservative in real arithmetic; outward nextafter elementary arithmetic and explicit numerical margin are used. Dense comparisons validate implementation, not a formal proof of all transcendental library rounding.'}
p=Path(__file__).parent/'results/wave_projection_bound_controls.json';p.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
print('failures',[x for x in checks if not x['passed']])
if not report['passed']:raise SystemExit(1)
