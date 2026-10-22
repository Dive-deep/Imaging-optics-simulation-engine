"""Small common-plane scan; retain only Pareto-improved qualified prescriptions."""
from pathlib import Path
import copy,json,time
import numpy as np
from scale_focus import normalize_height
from model import trace_design, native_surfaces, ORDERS
from radiometry import evaluate_design
from refine import summary_proxy
from aperture_fit import fit_apertures

E=Path(__file__).resolve().parent
OUT=E/'focus_similarity_best'

def safe_metrics(d):
 m=evaluate_design(d,image_height_tolerance_mm=.001)
 good=m['all_declared_numeric_constraints_and_sampled_rays_passed'] and min(m['NA_rgb_onaxis_conservative_complete_rim'])>=.7
 return m,good

def proxy_from_xy(xy,w):
 # Same equal-incident-RGB common-centroid proxy as refine.summary_proxy.
 total=3*np.sum(w);center=np.sum(xy*w[None,None,:,None],axis=(0,2))/total
 freq=np.linspace(0,20,81);phase=2*np.pi*(xy-center[None,:,None,:])[...,None]*freq
 otf=np.sum(np.exp(-1j*phase)*w[None,None,:,None,None],axis=(0,2))/total
 auc=np.trapezoid(abs(otf),freq,axis=-1)/20
 return {'mean':float(auc.mean()),'min':float(auc.min()),'negative':int(np.sum(otf.real<0))}

def main():
 start=time.perf_counter();OUT.mkdir(exist_ok=True)
 byseed={};review=[]
 for folder in ('qualified_round1','refined'):
  for path in sorted((E/'results'/folder).glob('seed_*.json')):
   raw=json.loads(path.read_text());d=raw['design']
   if d['aperture_mode']!='front_diameter50':continue
   m=raw.get('metrics')
   if m is None:m,good=safe_metrics(d)
   else:good=m['all_declared_numeric_constraints_and_sampled_rays_passed'] and min(m['NA_rgb_onaxis_conservative_complete_rim'])>=.7
   if not good:continue
   if 'aperture_fit'in raw:
    a=raw['aperture_fit']
    if a['original_transmitted_ray_count']!=a['sampled_ray_count']:continue
   else:
    dense=trace_design(d,fields=np.linspace(0,1,13),rings=16,arms=48)
    if not np.all((dense['status']==0)&dense['image_forward']):continue
   q=raw['geometric_mtf_proxy'];score=.5*(q['mean_AUC_normalized']+q['minimum_AUC_normalized'])
   row={'path':str(path),'raw':raw,'metrics':m,'score':score}
   sid=d['seed_id']
   review.append({'source':str(path),'seed_id':sid,'score':score,'mean_AUC':q['mean_AUC_normalized'],'minimum_AUC':q['minimum_AUC_normalized'],'negative_OTF_points':q['negative_OTF_points'],'RI1_bare':m['polychromatic_fields']['1.0']['RI_uncoated_fresnel']})
   if sid not in byseed or score>byseed[sid]['score']:byseed[sid]=row
 chosen=sorted(byseed.values(),key=lambda r:r['score'],reverse=True)[:8]
 results=[]
 for selected in chosen:
  original=selected['raw']['design'];sid=original['seed_id'];base=selected['raw']['geometric_mtf_proxy'];base_m=selected['metrics']
  t=trace_design(original,rings=16,arms=48)
  assert np.all((t['status']==0)&t['image_forward'])
  original_z=float(sum(original['axial_distances_mm']));h=t['image_xy'][1,-1,0,1]
  ratio=t['dir'][...,:2]/t['dir'][...,2,None]
  trial_rows=[]
  for delta in np.linspace(-.5,.7,49):
   shifted=t['image_xy']+delta*ratio
   scale=6.4/shifted[1,-1,0,1]
   distances=np.asarray(original['axial_distances_mm']).copy();distances[5]+=delta;distances*=scale
   if scale<=0 or not .8<scale<1.2 or min(distances[[0,2,4]])<=2 or min(distances[[1,3]])<1 or distances[5]<3:continue
   q=proxy_from_xy(shifted*scale,t['weights'])
   if q['mean']>base['mean_AUC_normalized']+1e-4 and q['min']>=base['minimum_AUC_normalized']-1e-5:
    trial_rows.append({'delta_mm':float(delta),'scale':float(scale),**q,'score':.5*(q['mean']+q['min'])})
  accepted=None
  for row in sorted(trial_rows,key=lambda r:r['score'],reverse=True)[:8]:
   shifted=copy.deepcopy(original);shifted['axial_distances_mm'][5]+=row['delta_mm']
   d,scale=normalize_height(shifted)
   m,good=safe_metrics(d)
   if not good:continue
   ri=m['polychromatic_fields']['1.0']['RI_uncoated_fresnel'];old_ri=base_m['polychromatic_fields']['1.0']['RI_uncoated_fresnel']
   if abs(ri-1)>abs(old_ri-1)+.002:continue
   d,audit=fit_apertures(d,coupled=False)
   if not audit['accepted'] or audit['original_transmitted_ray_count']!=audit['sampled_ray_count']:continue
   m,good=safe_metrics(d)
   q=summary_proxy(d,rings=16,arms=48)
   if not good or not q['all_rays_valid'] or q['mean_AUC_normalized']<=base['mean_AUC_normalized']+1e-4 or q['minimum_AUC_normalized']<base['minimum_AUC_normalized']-1e-5:continue
   out={'source_result':selected['path'],'design':d,'initial_design':original,'aperture_fit':audit,'metrics':m,'geometric_mtf_proxy':q,'focus_similarity_scan':{'method':'One common sensor shift for allfields/RGB, then global similarity toG12degheight6.4; frontclear reset25, allrays retraced','scan_offsets_mm':[-.5,.7],'scan_points':49,'selected_offset_mm':row['delta_mm'],'similarity_scale':scale,'baseline_mean_AUC':base['mean_AUC_normalized'],'baseline_minimum_AUC':base['minimum_AUC_normalized'],'baseline_RI1_bare':old_ri,'RI_deviation_allowed_extra':.002,'full_dense_summary_recomputed':True,'not_FFT_or_global_optimum':True}}
   path=OUT/(sid+'.json');path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
   accepted={'seed_id':sid,'path':str(path),'delta_mm':row['delta_mm'],'scale':scale,'old_mean':base['mean_AUC_normalized'],'new_mean':q['mean_AUC_normalized'],'old_min':base['minimum_AUC_normalized'],'new_min':q['minimum_AUC_normalized'],'RI1_bare':ri}
   break
  results.append(accepted or {'seed_id':sid,'improved_qualified_result':False,'eligible_proposals':len(trial_rows)})
  print(json.dumps(results[-1]),flush=True)
 summary={'seconds':time.perf_counter()-start,'source_review':review,'tested_distinct_seed_count':len(chosen),'results':results,'retention_rule':'Mean geometricMTFAUC improves>1e-4, minimumfield/axisAUC notlower>1e-5, RI1 deviation increases<=.002; allnumericconstraints, allRGBNA>=.7,13fieldallraytransmission pass. No per-color/fieldrefocus.'}
 (OUT/'scan_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 print('complete_seconds',summary['seconds'],flush=True)

if __name__=='__main__':main()
