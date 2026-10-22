"""Independent original-engine audit of the frozen CT20 candidate prescriptions (all distinct eligible parents).

This script changes no engine/model/native source. It samples chief, eight rim
rays, and sixteen interior rays at every RGB wavelength and four requested
fields, and applies the original engine's continuous radial interval geometry
check. The report is a sampled ray equivalence check, not wave accuracy proof.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time,traceback
import numpy as np

HERE=Path(__file__).resolve().parent
EXPERIMENT=HERE.parent
ROOT=EXPERIMENT.parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(EXPERIMENT))
from model import trace_design,index_normal_line,polar_samples
from sequential_optics.models import Surface,OpticalSystem,Ray
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_ray,validate_geometry

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def validate_selection(selection):
    """Reject stale labels or a truncated pool; never filter failed candidates."""
    records=selection['candidates']
    assert len(records)>=20, 'At least 20 frozen CT20 candidates are required'
    assert len({r['candidate_id'] for r in records})==len(records), 'Duplicate candidate ID'
    assert len({r['source_seed_id'] for r in records})==len(records), 'Duplicate parent'
    pool=selection['selection_pool']
    assert pool['all_distinct_eligible_parents_retained'] is True, 'Truncated parent pool'
    assert pool['distinct_eligible_parents']==pool['selected_distinct_parents']==len(records), 'Parent pool count mismatch'
    for record in records:
        design=record['design'];cid=record['candidate_id']
        assert record['source_seed_id']==design['seed_id'], cid+' parent identity mismatch'
        ct=np.asarray(design['axial_distances_mm'],float)[[0,2,4]]
        assert np.all(np.isfinite(ct)&(ct>2.)&(ct<=20.)), cid+' actual CT outside (2,20] mm'
        metric=record['metrics'];constraints=metric['constraint_results'];fit=record['aperture_fit']
        assert constraints['center_thickness_strictly_over2mm']['passed'] is True, cid+' missing fresh lower CT gate'
        assert constraints['center_thickness_at_most20mm']['passed'] is True, cid+' missing fresh upper CT gate'
        assert all(item.get('passed') is True for item in constraints.values()), cid+' failed constraint'
        assert metric['all_declared_numeric_constraints_and_sampled_rays_passed'] is True, cid+' failed sampled metric gate'
        na=np.asarray(metric['NA_rgb_onaxis_conservative_complete_rim'],float)
        assert na.shape==(3,) and np.all(np.isfinite(na)&(na>=.7)), cid+' invalid conservative RGB NA'
        assert fit['sampled_ray_count']>0 and fit['original_transmitted_ray_count']==fit['sampled_ray_count'], cid+' original clipping'
        digest=hashlib.sha256(json.dumps(design,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        assert digest==record['design_sha256'], cid+' design hash mismatch'
    return records

def original_system(design):
    """Construct original types independently from published design values."""
    media={'air':Material('air',n_constant=1.)}
    wavelengths=(486.,588.,656.)
    for name in design['materials']:
        media[name]=Material(name,wavelengths_nm=wavelengths,
                            indices=tuple(float(index_normal_line(name,w)) for w in wavelengths))
    surfaces=[];z=0.
    for j in range(6):
        material=design['materials'][j//2]
        coeff=tuple(float(design['asphere_sag_coefficients_mm'][j][k])/25.**order
                    for k,order in enumerate((4,6,8,10,12)))
        plastic=material in ('PMMA','PC')
        if not plastic and (design['conics'][j]!=0 or any(coeff)):
            raise AssertionError('Non-spherical glass prescription')
        surfaces.append(Surface(z=z,curvature=design['curvatures'][j],
                                conic=design['conics'][j],coefficients=(0.,*coeff),
                                clear_radius=design['clear_radii_mm'][j],
                                medium_after=material if j%2==0 else 'air',
                                surface_type='even_asphere' if plastic else 'standard'))
        z+=design['axial_distances_mm'][j]
    # External plane STOP has already selected launch locations. The original
    # first surface clear aperture is still physically checked by trace_ray.
    return OpticalSystem(tuple(surfaces),media,image_z=z,stop_index=0,
                         stop_radius=surfaces[0].clear_radius)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection',type=Path,default=EXPERIMENT/'results/frozen_selection.json')
    parser.add_argument('--output',type=Path,default=EXPERIMENT/'results/native_final_candidates')
    args=parser.parse_args();selection_file=args.selection.resolve();output=args.output.resolve()
    selection=json.loads(selection_file.read_text())
    files=validate_selection(selection)
    output.mkdir(parents=True,exist_ok=False)
    snapshots={str(selection_file.relative_to(ROOT)) if selection_file.is_relative_to(ROOT) else str(selection_file):sha(selection_file)}
    tracked=[ROOT/'src/sequential_optics/geometry.py',ROOT/'src/sequential_optics/models.py',
             ROOT/'src/sequential_optics/materials.py',EXPERIMENT/'model.py',
             HERE/'raytrace.py',HERE/'raytrace.c',HERE/'_raytrace_kernel.so',Path(__file__).resolve()]
    frozen={str(p.relative_to(ROOT)):sha(p) for p in tracked}
    xy,_=polar_samples(16,48)
    azimuths=(0,6,12,18,24,30,36,42)
    picked=np.array([0]+[1+ring*48+a for ring in (3,11) for a in azimuths]+
                    [1+16*48+a for a in azimuths])
    sample=xy[picked];fields=(0.,.5,.8,1.);waves=(486.,588.,656.)
    plan={'schema':'independent-original-engine-ct20-final-candidates-ray-audit/1','input_files_sha256':snapshots,
          'frozen_sources_sha256':frozen,'script_sha256':sha(Path(__file__)),
          'selection': f'All {len(files)} frozen distinct CT20 candidates, no refiltering or failed-case removal',
          'candidate_design_sha256':{r['candidate_id']:r['design_sha256'] for r in files},
          'candidate_count':len(files),'minimum_candidate_count':20,'fields':fields,'wavelengths_nm':waves,
          'base_sampling':'16rings48arms+chief+rim','selected_original_ray_indices':picked.tolist(),
          'normalized_pupil_points':sample.tolist(),'rays_per_candidate':len(sample)*len(fields)*len(waves),
          'position_and_opl_abs_limit_mm':1e-8,'direction_abs_limit':2e-10,
          'scope':'Original v0.2.1 trace_ray versus native model trace. Every selected ray retained. Original validate_geometry radialinterval across commonclearapertures. No proof of wave/MTF convergence or untested-ray equivalence.'}
    (output/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    case_reports=[];started=time.perf_counter()
    for case_index,record in enumerate(files):
        design=record['design'];sid=record['candidate_id']
        result={'candidate_id':sid,'source_seed_id':record['source_seed_id'],
                'source_file':(str(selection_file.relative_to(ROOT)) if selection_file.is_relative_to(ROOT) else str(selection_file))+'#'+sid,'source_sha256':sha(selection_file),
                'design_sha256':record['design_sha256'],
                'passed':False,'geometry_passed':False,'original_valid_rays':0,'native_valid_rays':0,
                'ray_count':len(sample)*len(fields)*len(waves),'failures':[],'max_position_difference_mm':0.,
                'max_direction_difference':0.,'max_opl_difference_mm':0.,
                'max_image_xy_difference_mm':0.,'worst_points':{},'rim_NA_comparison':{}}
        begin=time.perf_counter()
        try:
            system=original_system(design)
            result['geometry_original_interval']=validate_geometry(system)
            result['geometry_passed']=result['geometry_original_interval']['status']=='valid'
            native=trace_design(design,fields=fields,wavelengths=waves,xy=sample)
            result['native_valid_rays']=int(np.sum(native['status']==0))
            for wi,wave in enumerate(waves):
                original_rim_na=[];native_rim_na=[]
                for fi,field in enumerate(fields):
                    theta=np.deg2rad(12)*field;direction=np.array([0.,np.sin(theta),np.cos(theta)])
                    for pi,pupil in enumerate(sample):
                        origin=np.array([pupil[0]*design['stop_radius_mm'],pupil[1]*design['stop_radius_mm'],-design['stop_distance_mm']])
                        reference=trace_ray(system,Ray(origin,direction),wavelength_nm=wave)
                        location={'wave_nm':wave,'field_fraction':field,'selected_ray':pi,
                                  'original_ray_index':int(picked[pi]),'normalized_pupil_xy':pupil.tolist()}
                        ns=int(native['status'][wi,fi,pi])
                        if reference.success:result['original_valid_rays']+=1
                        if not reference.success or ns!=0:
                            result['failures'].append({**location,'original_status':reference.status,
                                                       'original_message':reference.message,'native_status':ns})
                            continue
                        vals={
                            'max_position_difference_mm':(native['history'][wi,fi,pi]-reference.positions,1e-8),
                            'max_direction_difference':(native['direction_history'][wi,fi,pi]-reference.directions,2e-10),
                            'max_opl_difference_mm':(native['opl_history'][wi,fi,pi]-reference.opl,1e-8),
                        }
                        dist=(system.image_z-reference.last_position[2])/reference.last_direction[2]
                        image=reference.last_position[:2]+dist*reference.last_direction[:2]
                        if not np.isfinite(dist) or dist < -1e-9 or not native['image_forward'][wi,fi,pi]:
                            result['failures'].append({**location,'reason':'non-forward or invalid image projection','original_projection_distance_mm':float(dist),'native_image_forward':bool(native['image_forward'][wi,fi,pi])})
                        vals['max_image_xy_difference_mm']=(native['image_xy'][wi,fi,pi]-image,1e-8)
                        for key,(difference,limit) in vals.items():
                            value=float(np.max(np.abs(difference)))
                            if value>result[key]:
                                result[key]=value;result['worst_points'][key]={**location,'maximum_absolute_difference':value}
                            if not np.isfinite(value) or value>limit:
                                result['failures'].append({**location,'quantity':key,'difference':value,'limit':limit})
                        if fi==0 and pi>=1+2*len(azimuths):
                            original_rim_na.append(float(np.linalg.norm(reference.last_direction[:2])))
                            native_rim_na.append(float(np.linalg.norm(native['dir'][wi,fi,pi,:2])))
                result['rim_NA_comparison'][str(wave)]={'original_min':min(original_rim_na,default=None),
                                                       'native_min':min(native_rim_na,default=None),
                                                       'sampled_rim_count':len(original_rim_na)}
            result['passed']=result['geometry_passed'] and not result['failures'] and result['original_valid_rays']==result['ray_count'] and result['native_valid_rays']==result['ray_count']
        except Exception as ex:
            result['failures'].append({'exception':str(ex),'traceback':traceback.format_exc()})
        result['seconds']=time.perf_counter()-begin
        (output/f'{sid}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        case_reports.append(result)
        print(json.dumps({'completed':case_index+1,'total':len(files),'candidate_id':sid,
                          'passed':result['passed'],'seconds':result['seconds']}),flush=True)
    after={str(p.relative_to(ROOT)):sha(p) for p in tracked}
    changed=[path for path,h in frozen.items() if after[path]!=h]
    changed_inputs=[path for path,h in snapshots.items() if sha(ROOT/path)!=h]
    summary={'schema':'independent-original-engine-ct20-final-candidates-ray-audit-summary/1','candidate_count':len(files),
             'total_rays':sum(r['ray_count'] for r in case_reports),
             'passed':all(r['passed'] for r in case_reports) and not changed and not changed_inputs,
             'failed_candidates':[r['candidate_id'] for r in case_reports if not r['passed']],
             'candidate_ids':[r['candidate_id'] for r in case_reports],
             'original_interval_geometry_passed_candidates':sum(r['geometry_passed'] for r in case_reports),
             'status_or_numeric_failure_count':sum(len(r['failures']) for r in case_reports),
             'changed_sources':changed,'changed_inputs':changed_inputs,'seconds':time.perf_counter()-started,
             'scope':plan['scope']}
    for key in ('max_position_difference_mm','max_direction_difference','max_opl_difference_mm','max_image_xy_difference_mm'):
        worst=max(case_reports,key=lambda r:r[key]);summary[key]=worst[key]
        summary[key+'_source']={'candidate_id':worst['candidate_id'],**worst['worst_points'].get(key,{})}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    if not summary['passed']:raise SystemExit(1)

if __name__=='__main__':main()
