"""Select distinct qualified starts and export explicit physical prescriptions.

Ranking is an engineering screen using geometric OTF and approximate RI, not a
claim that underresolved scalar wave curves establish the globally best lens.
"""
from pathlib import Path
import argparse, copy, hashlib, json
import numpy as np
from model import MATERIALS, WAVELENGTHS, COEFFS, index_normal_line, export_prescription

HERE = Path(__file__).resolve().parent


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def eligible(record):
    m, a = record['metrics'], record['aperture_fit']
    return (m['all_declared_numeric_constraints_and_sampled_rays_passed']
            and min(m['NA_rgb_onaxis_conservative_complete_rim']) >= .7
            and a['original_transmitted_ray_count'] == a['sampled_ray_count'])


def score(record):
    g = record['geometric_mtf_proxy']
    ri = record['metrics']['polychromatic_fields']['1.0']['RI_uncoated_fresnel']
    return .6*g['mean_AUC_normalized'] + .4*g['minimum_AUC_normalized'] - .25*abs(ri-1)


def history(path):
    steps=[];visited=set()
    for _ in range(12):
        p=Path(path)
        if not p.is_absolute():p=HERE.parents[1]/p
        if not p.exists() or str(p) in visited:break
        visited.add(str(p));raw=json.loads(p.read_text())
        step={'file':str(p.relative_to(HERE.parents[1])),
              'file_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
        for key in ('optimization_log','total_seconds','refinement','focus_similarity_scan'):
            if key in raw:step[key]=raw[key]
        if 'initial_seed' in raw:step['initial_seed']=raw['initial_seed']
        steps.append(step)
        if not raw.get('source_result'):break
        path=raw['source_result']
    return list(reversed(steps))


def select(folders, count=24):
    pool=[]
    for folder in folders:
        for path in sorted(Path(folder).glob('*.json')):
            r=json.loads(path.read_text())
            if not all(k in r for k in ('metrics','aperture_fit','geometric_mtf_proxy','design')):
                continue
            if eligible(r):
                pool.append((score(r), str(path), r))
    # A different optimization checkpoint is not counted as a new design family.
    best={}
    for value, path, r in sorted(pool, key=lambda x:(-x[0],x[1])):
        best.setdefault(r['design']['seed_id'], (value,path,r))
    chosen=list(best.values())[:count]
    if len(chosen)<count:
        raise ValueError(f'Only {len(chosen)} distinct eligible parents, requested {count}')
    output=[]
    for rank,(value,path,r) in enumerate(chosen,1):
        d=copy.deepcopy(r['design']); cid=f'C{rank:02d}'
        out={'candidate_id':cid,'source_seed_id':d['seed_id'],'design':d,
            'design_sha256':canonical_hash(d),
            'input':{'object':'infinity','field_angles_deg':[0.,6.,9.6,12.],
                     'full_symmetric_FOV_deg':24.,'image_diameter_mm':12.8,
                     'wavelengths_nm':WAVELENGTHS.tolist(),
                     'stop':{'z_mm':-d['stop_distance_mm'],'radius_mm':d['stop_radius_mm'],
                             'aperture_mode':d['aperture_mode']},
                     'image_plane_z_mm':float(sum(d['axial_distances_mm'])),
                     'surface_sag_convention':'z=c*r^2/(1+sqrt(1-(1+k)*c^2*r^2))+sum(A_m*r^m), m=4,6,8,10,12',
                     'radius_null_means':'plane (curvature=0)',
                     'surfaces':export_prescription(d)},
            'metrics':r['metrics'],'aperture_fit':r['aperture_fit'],
            'geometric_mtf_proxy':r['geometric_mtf_proxy'],
            'selection':{'rank':rank,'score':value,'source_result':path,
                         'score_formula':'0.6*mean_geometric_MTF_AUC/20 + 0.4*minimum_geometric_MTF_AUC/20 - 0.25*abs(RI_edge_uncoated-1)',
                         'interpretation':'preliminary engineering screen; not globally optimal or a wave-MTF accuracy ranking'},
            'provenance':{'seed_generation':d.get('provenance',{}),'checkpoint':r.get('source_result',path),
                          'generation_history':history(path)},
            'qualification':{'geometric_sampled_constraints':'passed',
                             'full_13_field_pupil_samples_all_transmitted':True,
                             'all_three_wavelength_onaxis_complete_rim_NA_at_least_0_7':True,
                             'height_reference_wavelength_nm':588.,'height_tolerance_mm':.001,
                             'manufacturing_qualified':False,'Zemax_validated':False,
                             'surrogate_training_label_ready':False}}
        output.append(out)
    return {'schema':'three_element_initial_candidates_v1','created_date':'2026-10-22',
            'units':{'length':'mm','wavelength':'nm','angle':'degree','spatial_frequency':'line pairs/mm'},
            'purpose':'24 initial prescriptions for subsequent local/global design-space exploration',
            'scope':'Geometric constraints are sampled and qualified. Wave output has separate conditional/underresolved status. No global-optimum, manufacturing, or Zemax equivalence claim.',
            'materials':{name:{'nd':nd,'Vd':vd,'nd_reference_nm':587.5618,
                               'dispersion':'normal_line_three_term_Cauchy',
                               'cauchy_coefficients_lambda_um':COEFFS[name].tolist(),
                               'n_at_requested_wavelengths':{str(int(w)):float(index_normal_line(name,w)) for w in WAVELENGTHS}}
                         for name,(nd,vd) in MATERIALS.items()},
            'selection_pool':{'eligible_checkpoints':len(pool),'distinct_eligible_parents':len(best),
                              'selected_distinct_parents':len(output),'folders':[str(p) for p in folders]},
            'candidates':output}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folders',nargs='+',type=Path)
    p.add_argument('--count',type=int,default=24);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();result=select(a.folders,a.count)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False))
    print(json.dumps(result['selection_pool']))
