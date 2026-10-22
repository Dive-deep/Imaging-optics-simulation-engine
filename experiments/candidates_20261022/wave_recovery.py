"""Read-only diagnosis of physical boundary-map folds; no guard bypass."""
from pathlib import Path
import json,hashlib,time
import numpy as np
from scipy.spatial import Delaunay
from wave_evaluation import build_system,native_bundle,conservative_surface_projection_upper

HERE=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cross_area(p):
 a=p[:,1]-p[:,0];b=p[:,2]-p[:,0]
 return a[:,0]*b[:,1]-a[:,1]*b[:,0]

def diagnose(design,field,wavelength=588.,rings=16,arms=48):
 b,raw=native_bundle(design,field,wavelength,rings=rings,arms=arms);system=build_system(design)
 good=b.success_mask;p=raw['pos'];d=raw['dir'];e3=b.chief.last_direction
 e1=np.array([1.,0.,0.])-e3[0]*e3;e1/=np.linalg.norm(e1);e2=np.cross(e3,e1);basis=np.array([e1,e2])
 pc=b.chief.last_position;q0=pc+(system.image_z-pc[2])/e3[2]*e3
 full=conservative_surface_projection_upper(system.surfaces[-1],e3)+.001
 actual=float(np.max(p[good]@e3))+.001
 image=float(q0@e3)
 topology=Delaunay(b.pupil_xy).simplices;tri=topology[np.all(good[topology],axis=1)]
 area0=cross_area(b.pupil_xy[tri]);cosines=d@e3
 def map_check(points):
  jac=cross_area(points[tri])/area0
  return {'positive_triangles':int(np.count_nonzero(jac>0)),
   'negative_triangles':int(np.count_nonzero(jac<0)),
   'jacobian_min':float(jac.min()),'jacobian_max':float(jac.max()),
   'fold_detected':bool(np.any(np.sign(jac)!=np.sign(jac[0])))}
 surface_checks=[{'surface_index_1_based':i+1,**map_check(raw['history'][:,i,:2])}
  for i in range(raw['history'].shape[1])]
 direction_check=map_check(d[:,:2])
 projections=list(dict.fromkeys([actual,full,*np.linspace(actual,full,5),*np.linspace(actual,image,7)]))
 checks=[]
 for proj in projections:
  origin=pc+(proj-pc@e3)*e3
  t=(proj-p@e3)/cosines
  uv=(p+t[:,None]*d-origin)@basis.T
  areas=cross_area(uv[tri]);jac=areas/area0
  checks.append({'projection_mm':float(proj),'reference_distance_mm':float(image-proj),
    'after_all_sampled_exit_points':bool(np.all(t[good]>=-1e-10)),
    'full_surface_conservative_plane':bool(proj>=full-1e-10),
    'positive_triangles':int(np.count_nonzero(jac>0)),
    'negative_triangles':int(np.count_nonzero(jac<0)),
    'jacobian_min':float(jac.min()),'jacobian_max':float(jac.max()),
    'fold_detected':bool(np.any(np.sign(jac)!=np.sign(jac[0]))),
    'minimum_ray_boundary_cosine':float(cosines[good].min())})
 return {'field':field,'wavelength_nm':wavelength,'rings':rings,'arms':arms,'rays':len(good),
  'valid_rays':int(good.sum()),'triangle_count':len(tri),'full_surface_projection_mm':full,
  'sampled_exit_projection_mm':actual,'image_projection_mm':image,'plane_checks':checks,
  'physical_surface_xy_map_checks':surface_checks,
  'output_direction_xy_map_check':direction_check}

if __name__=='__main__':
 start=time.perf_counter();selection=HERE/'results/frozen_selection.json';data=json.loads(selection.read_text())
 case=next(x for x in data['candidates'] if x['candidate_id']=='C01')
 before=digest(HERE/'wave_evaluation.py');out=[]
 for field in (.8,1.):
  for rings,arms in ((16,48),(32,96),(64,192)):
   z=diagnose(case['design'],field,rings=rings,arms=arms);out.append(z)
   print(field,rings,arms,'full',[p for p in z['plane_checks'] if abs(p['projection_mm']-z['full_surface_projection_mm'])<1e-10],flush=True)
 report={'candidate_id':'C01','source_selection_sha256':digest(selection),'wave_source_sha256_before':before,
  'wave_source_sha256_after':digest(HERE/'wave_evaluation.py'),'diagnostics':out,'seconds':time.perf_counter()-start,
  'scope':'No optical prescriptions, frozen wave engine, or final wave results changed. Sampled-footprint planes are diagnostic only, not accepted physical recovery.',
  'conclusion':'C01 0.8F/1F at 588 nm has persistent signed-Jacobian folds in the actual last-surface XY ray map and all tested posterior physical planes. Refining the STOP mesh or choosing a plane after only the sampled exit footprint does not restore a regular spatial map. The existing single-valued boundary profile correctly rejects these bundles. Output-direction-map regularity suggests a possible future momentum-chart approach but does not validate its phase/amplitude/Maslov physics. Preserve current guard failures and show separately labeled geometric MTF diagnostics.'}
 (HERE/'results/wave_recovery_C01_diagnosis.json').write_text(json.dumps(report,indent=2)+'\n')
 print('seconds',report['seconds'],flush=True)
