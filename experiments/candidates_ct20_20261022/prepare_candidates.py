"""Preserve all distinct CT20-qualified parents as an unranked candidate set.

The engineering screen only chooses among checkpoints of the same parent.
Presentation order and persistent IDs do not express an optical quality rank.
"""
from pathlib import Path
import argparse, copy, hashlib, json, re
import numpy as np
from model import MATERIALS, WAVELENGTHS, COEFFS, index_normal_line, export_prescription

HERE = Path(__file__).resolve().parent


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def eligible(record):
    m, a = record['metrics'], record['aperture_fit']
    ct=np.asarray(record['design']['axial_distances_mm'],float)[[0,2,4]]
    constraints=m.get('constraint_results',{})
    na=np.asarray(m['NA_rgb_onaxis_conservative_complete_rim'],float)
    return bool(m['all_declared_numeric_constraints_and_sampled_rays_passed']
            and np.all(np.isfinite(ct)&(ct>2)&(ct<=20))
            and constraints.get('center_thickness_strictly_over2mm',{}).get('passed') is True
            and constraints.get('center_thickness_at_most20mm',{}).get('passed') is True
            and all(item.get('passed') is True for item in constraints.values())
            and na.shape==(3,) and np.all(np.isfinite(na)&(na>=.7))
            and a['original_transmitted_ray_count'] == a['sampled_ray_count']
            and a['sampled_ray_count']>0)


def score(record):
    g = record['geometric_mtf_proxy']
    ri = record['metrics']['polychromatic_fields']['1.0']['RI_uncoated_fresnel']
    return .6*g['mean_AUC_normalized'] + .4*g['minimum_AUC_normalized'] - .25*abs(ri-1)


def history(path):
    """Walk source and embedded parent links, oldest first, without losing CT repair."""
    steps=[];visited=set();repo=HERE.parents[1]
    def walk(value,relative_to=repo):
        p=Path(value)
        if not p.is_absolute():
            root_path=repo/p
            p=root_path if root_path.exists() else relative_to/p
        p=p.resolve()
        if p in visited:return
        visited.add(p)
        label=str(p.relative_to(repo)) if p.is_relative_to(repo) else str(p)
        if not p.exists():
            steps.append({'file':label,'status':'referenced_source_missing'});return
        raw=json.loads(p.read_text());links=[]
        if raw.get('source_result'):links.append(raw['source_result'])
        for key in ('parent_input','old_baseline','baseline'):
            nested=raw.get(key)
            if isinstance(nested,dict) and nested.get('source_result'):links.append(nested['source_result'])
            elif isinstance(nested,str):links.append(nested)
        if raw.get('old_baseline_source_result'):links.append(raw['old_baseline_source_result'])
        for source in links:walk(source,p.parent)
        step={'file':label,'file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
              'source_results':list(dict.fromkeys(links)),
              'phase':'ct20_revision' if p.is_relative_to(HERE) else 'preserved_prior_history'}
        for key in ('optimization_log','total_seconds','refinement','focus_similarity_scan',
                    'ct20_continuation','ct20_polish','ct20_rescue','ct20_stop_scan',
                    'common_image_plane_delta_mm','initial_seed'):
            if key in raw:step[key]=raw[key]
        if isinstance(raw.get('parent_input'),dict):
            parent=raw['parent_input'];step['parent_input']={k:parent[k] for k in ('parent_seed_id','source_result','source_sha256') if k in parent}
        steps.append(step)
    walk(path)
    return steps


def persistent_id(seed_id):
    match=re.fullmatch(r'seed_(\d+)',str(seed_id))
    suffix=f'{int(match.group(1)):04d}' if match else re.sub(r'[^a-zA-Z0-9_.-]','_',str(seed_id))
    return 'CT20_'+suffix


def presentation_key(item):
    record=item[2];materials=record['design']['materials']
    glass=next(i for i,name in enumerate(materials) if name not in ('PMMA','PC'))
    return glass,tuple(materials),str(record['design']['seed_id'])


def select(folders, count=20):
    # Legacy --count is retained as a minimum, never as a truncation limit.
    if count<20:raise ValueError('At least 20 distinct qualified parents are required')
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
    chosen=sorted(best.values(),key=presentation_key)
    if len(chosen)<count:
        raise ValueError(f'Only {len(chosen)} distinct eligible parents, requested {count}')
    output=[]
    for value,path,r in chosen:
        d=copy.deepcopy(r['design']); cid=persistent_id(d['seed_id'])
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
            'selection':{'score':value,'score_role':'same-parent checkpoint screening diagnostic only','source_result':path,
                         'score_formula':'0.6*mean_geometric_MTF_AUC/20 + 0.4*minimum_geometric_MTF_AUC/20 - 0.25*abs(RI_edge_uncoated-1)',
                         'interpretation':'all eligible distinct parents retained; score chooses only a same-parent checkpoint, not an overall ranking',
                         'presentation_order':'glass position, material order, source seed ID'},
            'provenance':{'seed_generation':d.get('provenance',{}),'checkpoint':r.get('source_result',path),
                          'generation_history':history(path)},
            'qualification':{'geometric_sampled_constraints':'passed',
                             'each_lens_center_thickness_strictly_over2_and_at_most20_mm':True,
                             'full_13_field_pupil_samples_all_transmitted':True,
                             'all_three_wavelength_onaxis_complete_rim_NA_at_least_0_7':True,
                             'height_reference_wavelength_nm':588.,'height_tolerance_mm':.001,
                             'manufacturing_qualified':False,'Zemax_validated':False,
                             'surrogate_training_label_ready':False}}
        output.append(out)
    if len({r['candidate_id'] for r in output})!=len(output):raise ValueError('Persistent candidate ID collision')
    return {'schema':'three_element_initial_candidates_ct20_v2','created_date':'2026-10-22',
            'units':{'length':'mm','wavelength':'nm','angle':'degree','spatial_frequency':'line pairs/mm'},
            'purpose':f'All {len(output)} distinct eligible CT20 initial prescriptions for subsequent design-space exploration; no overall winner',
            'center_thickness_requirement_mm':{'exclusive_lower':2.,'inclusive_upper':20.},
            'scope':'Geometric constraints are sampled and qualified. Wave output has separate conditional/underresolved status. No global-optimum, manufacturing, or Zemax equivalence claim.',
            'materials':{name:{'nd':nd,'Vd':vd,'nd_reference_nm':587.5618,
                               'dispersion':'normal_line_three_term_Cauchy',
                               'cauchy_coefficients_lambda_um':COEFFS[name].tolist(),
                               'n_at_requested_wavelengths':{str(int(w)):float(index_normal_line(name,w)) for w in WAVELENGTHS}}
                         for name,(nd,vd) in MATERIALS.items()},
            'selection_pool':{'eligible_checkpoints':len(pool),'distinct_eligible_parents':len(best),
                              'selected_distinct_parents':len(output),'minimum_required_parents':count,
                              'all_distinct_eligible_parents_retained':len(output)==len(best),
                              'presentation_order':'glass position, material order, source seed ID',
                              'folders':[str(p) for p in folders]},
            'candidates':output}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folders',nargs='+',type=Path)
    p.add_argument('--count',type=int,default=20,help='Minimum required distinct parents; all eligible parents are retained')
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();result=select(a.folders,a.count)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False))
    print(json.dumps(result['selection_pool']))
