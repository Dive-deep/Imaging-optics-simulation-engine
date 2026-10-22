"""Reproducible aperture-continuation least-squares seed refinement.
Every launch ray is retained in residuals. No claim of a global optimum.
"""
from pathlib import Path
import argparse, copy, json, time, sys
import numpy as np
from scipy.optimize import least_squares
from model import native_surfaces,sag,trace_design,first_order,geometry_checks
E=Path(__file__).resolve().parent
class Parameterization:
 def __init__(self,seed,orders=2,free_conic=True):
  self.seed=copy.deepcopy(seed);self.plastic=[j for j in range(6) if seed['materials'][j//2] in ('PMMA','PC')]
  self.orders=orders;self.free_conic=free_conic
  self.names=[];self.lower=[];self.upper=[];self.initial=[]
  def add(name,val,lo,hi):self.names.append(name);self.initial.append(val);self.lower.append(lo);self.upper.append(hi)
  for j in range(6):
   bound=2.5 if j in self.plastic else .98*25/seed['clear_radii_mm'][j]
   add(('c',j),seed['curvatures'][j]*25,-bound,bound)
  if free_conic:
   for j in self.plastic:add(('k',j),seed['conics'][j],-40.,2.)
  for j in self.plastic:
   for order in range(orders):add(('b',j,order),seed['asphere_sag_coefficients_mm'][j][order],-30.,30.)
  for j in range(6):
   lo,hi=(2.02,35.) if j in (0,2,4) else ((1.,20.) if j in (1,3) else (3.,50.))
   add(('d',j),seed['axial_distances_mm'][j],lo,hi)
  add(('stop',),seed['stop_distance_mm'],.1,28.)
  if seed['aperture_mode']=='front_diameter50':add(('a',),seed['stop_radius_mm'],21.6,24.8)
  self.lower=np.array(self.lower);self.upper=np.array(self.upper);self.initial=np.clip(self.initial,self.lower+1e-7,self.upper-1e-7)
 def unpack(self,x):
  d=copy.deepcopy(self.seed)
  for name,val in zip(self.names,x):
   if name[0]=='c':d['curvatures'][name[1]]=float(val/25)
   elif name[0]=='k':d['conics'][name[1]]=float(val)
   elif name[0]=='b':d['asphere_sag_coefficients_mm'][name[1]][name[2]]=float(val)
   elif name[0]=='d':d['axial_distances_mm'][name[1]]=float(val)
   elif name[0]=='stop':d['stop_distance_mm']=float(val)
   else:d['stop_radius_mm']=float(val)
  return d

def residual(design,pupil_fraction=1.,fields=(0.,.5,1.),wavelengths=(486.,588.,656.),rings=3,arms=8,na_target=.71):
 t=trace_design(design,fields,wavelengths,rings=rings,arms=arms,pupil_fraction=pupil_fraction)
 img=np.clip(t['image_xy'],-100.,100.);good=(t['status']==0)&t['image_forward']
 # Absolute desired mapping shared by wavelengths: near-sine mapping limits RI loss.
 target=6.4*np.sin(np.radians(12)*np.array(fields))/np.sin(np.radians(12))
 desired=np.stack([np.zeros_like(target),target],axis=-1)
 rayres=(img-desired[None,:,None,:])
 # Uniform weights across all sampled rays including pupil edge (not performance quadrature).
 rayres=rayres/np.sqrt(img.shape[-2]);rayres=np.nan_to_num(rayres,nan=30.,posinf=30.,neginf=-30.)
 bad=(~good).astype(float)*2.
 chief=(img[:,:,0]-desired[None,:,:])*4.
 # NA constraint via positive/negative y rim, evaluated on axis from actual image directions.
 dirs=t['dir'][:,0];na=np.sqrt(np.sum(dirs[...,:2]**2,axis=-1));actual_na=np.max(np.where(good[:,0],na,0.),axis=1)
 nares=np.maximum(0,na_target*pupil_fraction-actual_na)*8.
 s=native_surfaces(design);geom=[]
 for j in range(6):
  rr=np.linspace(0,s[j,8],17);ss,disc=sag(s[[j]],rr)
  geom.extend(np.maximum(0,.02-disc[0])*2.)
  # Avoid high-order shape excursions that bound numerical but invalid lens geometry.
  geom.extend(np.maximum(0,np.abs(ss[0])-55.)*.5)
 for j in range(5):
  rr=np.linspace(0,min(s[j,8],s[j+1,8]),25);ss,dd=sag(s[[j,j+1]],rr)
  sep=s[j+1,0]+ss[1]-s[j,0]-ss[0]
  minimum=.4 if j%2==0 else 1.
  geom.extend(np.maximum(0,minimum-sep)*2.)
 ss,_=sag(s[[0,5]],np.linspace(0,min(s[0,8],s[5,8]),17))
 geom.append(max(0,.1-design['stop_distance_mm']-np.min(ss[0]))*2.)
 # Full last clear aperture matters for physical sensor clearance.
 back,_=sag(s[[5]],np.linspace(0,s[5,8],17));geom.append(max(0,.1+back.max()-design['axial_distances_mm'][5])*2.)
 # Stop should remain less than paraxial EFL; smooth weak guard where positive.
 fo=first_order(design);geom.append(max(0,design['stop_distance_mm']-max(1.,fo['efl_mm']))*.1)
 # Small high-order sag norm prevents runaway cancellations, without forcing their omission.
 b=np.array(design['asphere_sag_coefficients_mm']);regularizer=b[:,2:].ravel()*.0005
 return np.r_[rayres.ravel(),bad.ravel(),chief.ravel(),nares,geom,regularizer]

def quick_metrics(design,rings=6,arms=16):
 t=trace_design(design,rings=rings,arms=arms);valid=(t['status']==0)&t['image_forward'];w=t['weights'];ww=valid*w
 cent=np.sum(t['image_xy']*ww[...,None],axis=2)/np.maximum(ww.sum(axis=2)[...,None],1e-30)
 rms=np.sqrt(np.sum(np.sum((t['image_xy']-cent[:,:,None,:])**2,axis=-1)*ww,axis=2)/np.maximum(ww.sum(axis=2),1e-30))
 na=np.max(np.where(valid[:,0],np.sqrt(np.sum(t['dir'][:,0,:,:2]**2,axis=-1)),0.),axis=1)
 return {'rms_mm':rms.tolist(),'valid_fraction':(np.sum(ww,axis=2)).tolist(),'NA_rgb_axis':na.tolist(),
         'chief_full_height_mm':t['image_xy'][:,-1,0,1].tolist(),'image_z_mm':float(sum(design['axial_distances_mm'])),
         'geometry':geometry_checks(design),'first_order':first_order(design),'failed_ray_count':int(np.sum(~valid))}

def optimize_seed(seed,*,max_nfev=100,stages=(.45,.65,.85,1.),full_order=False):
 d=copy.deepcopy(seed);log=[];start=time.perf_counter()
 for stage_index,frac in enumerate(stages):
  order=5 if full_order and frac>=.99 else 2
  p=Parameterization(d,orders=order)
  before=time.perf_counter();calls=[0]
  def fun(x):calls[0]+=1;return residual(p.unpack(x),pupil_fraction=frac)
  fit=least_squares(fun,p.initial,bounds=(p.lower,p.upper),max_nfev=max_nfev,
          x_scale='jac',ftol=2e-5,xtol=2e-5,gtol=2e-5,diff_step=2e-5)
  d=p.unpack(fit.x)
  log.append({'pupil_fraction':frac,'cost':float(fit.cost),'nfev':fit.nfev,'actual_residual_calls':calls[0],
              'seconds':time.perf_counter()-before,'status':int(fit.status),'message':fit.message})
  print(json.dumps({'seed':seed['seed_id'],'stage':log[-1]}),flush=True)
 return {'design':d,'optimization_log':log,'total_seconds':time.perf_counter()-start,'quick_metrics':quick_metrics(d)}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seeds',type=Path,default=E/'seeds.json');ap.add_argument('--ids',nargs='*');ap.add_argument('--limit',type=int);ap.add_argument('--max-nfev',type=int,default=100);ap.add_argument('--stages',type=float,nargs='+',default=[.45,.65,.85,1.]);ap.add_argument('--full-order',action='store_true');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 seeds=json.loads(args.seeds.read_text());seeds=seeds.get('seeds',seeds) if isinstance(seeds,dict) else seeds
 if args.ids:seeds=[s for s in seeds if s['seed_id'] in args.ids]
 if args.limit:seeds=seeds[:args.limit]
 args.output.mkdir(parents=True,exist_ok=True)
 for seed in seeds:
  out=args.output/(seed['seed_id']+'.json')
  if out.exists():continue
  try:r=optimize_seed(seed,max_nfev=args.max_nfev,stages=args.stages,full_order=args.full_order);out.write_text(json.dumps(r,indent=2,allow_nan=False))
  except Exception as exc:
   import traceback;traceback.print_exc();(args.output/(seed['seed_id']+'.error.txt')).write_text(str(exc))
if __name__=='__main__':main()
