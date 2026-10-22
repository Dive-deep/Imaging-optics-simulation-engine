"""Representative-only final-stage ray-mesh refinement of scalar MTF.

Uses already generated candidate snapshots without changing prescriptions.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
from wave_evaluation import evaluate_one

HERE=Path(__file__).resolve().parent
INPUTS=['results/refine_pilot/seed_0012.json','results/qualified_round1/seed_0108.json']
SOURCES=['wave_evaluation.py','model.py','native/raytrace.py','native/raytrace.c','native/_raytrace_kernel.so','probe_wave_mesh.py']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def frozen():return {p:sha(HERE/p) for p in SOURCES}
def curve_summary(record):
    if record.get('mtf_t') is None:return {'status':record['status'],'reason':record.get('reason')}
    f=np.asarray(record['frequencies_lp_mm'])
    out={'status':record['status'],'quadrature_strict_relative_005':record['quadrature_convergence_strict_relative_005'],
         'quadrature_max_absolute_change':record['quadrature_history'][-1]['max_mtf_absolute_change'],
         'quadrature_max_relative_change_percent':record['quadrature_history'][-1]['max_mtf_relative_change_percent_nonzero'],
         'original_fft_status':record.get('actual_fft_diagnostic',{}).get('status')}
    for axis in ('t','s'):
        v=np.asarray(record['mtf_'+axis]);out[axis]={'mtf_20':float(v[-1]),'normalized_auc_0_20':float(np.trapezoid(v,f)/20),
                                                  'minimum_0_20':float(v.min())}
    return out

def compare(coarse,fine):
    if coarse.get('mtf_t') is None or fine.get('mtf_t') is None:return {'comparable':False,'coarse_status':coarse['status'],'fine_status':fine['status']}
    result={'comparable':True,'relative_reference':'finer32x96 mesh; every nonzero point, no absolute tolerance floor','axes':{}}
    passed=True
    for axis in ('t','s'):
        a=np.asarray(coarse['mtf_'+axis]);b=np.asarray(fine['mtf_'+axis]);error=abs(b-a);nz=b!=0
        relative=error[nz]/abs(b[nz]);fail=(relative>.005).sum();zero=((~nz)&(a!=0)).sum()
        passed=passed and fail==0 and zero==0
        sf=curve_summary(fine)[axis];sc=curve_summary(coarse)[axis]
        result['axes'][axis]={'max_absolute_change':float(error.max()),
          'max_relative_change_percent':float(100*relative.max()) if len(relative) else 0.,
          'strict_relative_failure_points':int(fail),'changed_zero_points':int(zero),
          'coarse_summary':sc,'fine_summary':sf,
          'auc_absolute_change':abs(sf['normalized_auc_0_20']-sc['normalized_auc_0_20']),
          'mtf20_absolute_change':abs(sf['mtf_20']-sc['mtf_20'])}
    result['strict_relative_005_all_points_passed']=bool(passed)
    return result

started=time.perf_counter();before=frozen();records=[]
for item in INPUTS:
    path=HERE/item;data=json.loads(path.read_text());design=data['design']
    case={'input':item,'input_sha256':sha(path),'seed_id':design['seed_id'],'fields':[]}
    for field in (0.,.8):
        pair=[]
        for rings,arms in ((16,48),(32,96)):
            wave=evaluate_one(design,field,588.,rings=rings,arms=arms,quadrature=(128,256),
                              fft_diagnostic=(rings==16),fft_grid=1024)
            pair.append(wave)
            print(design['seed_id'],field,rings,arms,wave['status'],'time',wave['seconds'],flush=True)
        case['fields'].append({'field_fraction':field,'wavelength_nm':588.,
            'coarse_summary':curve_summary(pair[0]),'fine_summary':curve_summary(pair[1]),
            'comparison':compare(*pair),'coarse_wave':pair[0],'fine_wave':pair[1]})
    records.append(case)
after=frozen()
report={'schema':'representative-wave-ray-mesh-refinement/1','representative_only':True,
 'scope':'Two saved design snapshots,588nm,0/.8F. Final selection may differ. This is not convergence validation of every delivered candidate or RGB.',
 'settings':{'coarse_rings_arms':[16,48],'fine_rings_arms':[32,96],'angular_quadrature_grids':[128,256],
             'actual_fft_coarse_only':{'grid':1024,'padding':2},'frequency_lp_mm':[0,20,.5]},
 'sources_before':before,'sources_after':after,'source_unchanged':before==after,
 'cases':records,'seconds':time.perf_counter()-started,
 'interpretation':'Scalar and local-sensor-projection model only. Underresolved curves can be displayed as numerical screening diagnostics, not trusted training labels. A stable AUC does not waive pointwise failures.'}
(HERE/'results/wave_mesh_probe.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'seconds':report['seconds'],'source_unchanged':report['source_unchanged'],
 'comparisons':[{'seed':c['seed_id'],'field':f['field_fraction'],**f['comparison']} for c in records for f in c['fields']]},indent=2))
assert before==after
