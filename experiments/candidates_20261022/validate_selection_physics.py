"""Fresh independent physics/serialization audit of a frozen candidate selection.

This does not trust stored qualification booleans. Re-traces all13 field points
and RGB, examines every marginal ray, independently reconstructs the published
radius/coefficient surface table, and samples geometry4097 times per interval.
It preserves color-dependent image heights and approximate RI limitations.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

from model import trace_design, native_surfaces, ORDERS
from radiometry import evaluate_design


def digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def direct_sag(row, radii):
    q=radii*radii;c=row[1];k=row[2]
    with np.errstate(invalid='ignore'):
        result=c*q/(1+np.sqrt(1-(1+k)*c*c*q))
    for coefficient,order in zip(row[3:8],ORDERS):
        result=result+coefficient*radii**order
    return result


def reconstruct_export(candidate):
    published=candidate['input'];surfaces=published['surfaces'];d=candidate['design']
    rows=[];fail=[];z=0.;distances=[]
    for j,s in enumerate(surfaces):
        radius=s['radius_mm'];c=0. if radius is None else 1/float(radius)
        plastic=d['materials'][j//2] in ('PMMA','PC')
        co=s.get('aspheric_coefficients',{})
        coefficients=[float(co.get(f'A{order}_mm_pow_{1-order}',0.)) for order in ORDERS]
        expected_keys={f'A{order}_mm_pow_{1-order}' for order in ORDERS} if plastic else set()
        expected_after=d['materials'][j//2] if j%2==0 else 'air'
        if set(co)!=expected_keys:fail.append(f'surface{j+1}:unsupported_or_missing_coefficient_keys')
        if s['material_after']!=expected_after:fail.append(f'surface{j+1}:material_after')
        if s['type']!=('even_asphere' if plastic else 'standard'):fail.append(f'surface{j+1}:surface_type')
        if not plastic and (s.get('conic',0)!=0 or any(coefficients)):fail.append(f'surface{j+1}:glass_not_standard')
        if abs(s['vertex_z_mm']-z)>1e-11:fail.append(f'surface{j+1}:vertex_thickness_chain')
        if abs(s['curvature_per_mm']-c)>1e-14:fail.append(f'surface{j+1}:radius_curvature_disagreement')
        rows.append([s['vertex_z_mm'],c,s.get('conic',0),*coefficients,s['clear_semi_diameter_mm']])
        distances.append(s['thickness_to_next_surface_mm']);z+=s['thickness_to_next_surface_mm']
    rows=np.asarray(rows,float);expected=native_surfaces(d)
    error=float(np.max(np.abs(rows-expected)))
    if error>1e-13:fail.append('exported_radius_coefficients_do_not_roundtrip')
    if np.max(np.abs(np.asarray(distances)-d['axial_distances_mm']))>1e-13:fail.append('distance_array_roundtrip')
    if abs(z-published['image_plane_z_mm'])>1e-11:fail.append('image_plane_thickness_chain')
    if abs(published['stop']['z_mm']+d['stop_distance_mm'])>1e-13:fail.append('STOP_z_roundtrip')
    if abs(published['stop']['radius_mm']-d['stop_radius_mm'])>1e-13:fail.append('STOP_radius_roundtrip')
    return {'passed':not fail,'failures':fail,'maximum_native_row_difference':error,
            'surface_count':len(surfaces),'reconstructed_image_plane_z_mm':z},rows


def audit(path):
    raw=json.loads(path.read_text());candidates=raw['candidates'];records=[]
    global_fail=[]
    if len(candidates)<20:global_fail.append('fewer_than20_candidates')
    parents=[c['source_seed_id'] for c in candidates]
    hashes=[digest(c['design']) for c in candidates]
    if len(set(parents))!=len(candidates):global_fail.append('duplicate_parent_seed')
    if len(set(hashes))!=len(candidates):global_fail.append('duplicate_numeric_design')
    if len({c['candidate_id'] for c in candidates})!=len(candidates):global_fail.append('duplicate_candidate_id')
    scores=[]
    for c,sha in zip(candidates,hashes):
        d=c['design'];issues=[]
        if sha!=c['design_sha256']:issues.append('design_hash_mismatch')
        roundtrip,surfaces=reconstruct_export(c)
        issues.extend(roundtrip['failures'])
        if len(d['materials'])!=3 or sum(m in ('PMMA','PC') for m in d['materials'])!=2 or sum(m in ('FDS18W','TAFD45') for m in d['materials'])!=1:
            issues.append('three_element_material_contract')
        distances=np.asarray(d['axial_distances_mm'])
        if not np.all(distances[[0,2,4]]>2):issues.append('center_thickness_not_strict_over2')
        if not np.all(distances[[1,3]]>=1):issues.append('axis_gap_below1')
        if not distances[5]>=3:issues.append('vertex_BFL_below3')
        gaps=[];domains=[]
        for j,row in enumerate(surfaces):
            radii=np.linspace(0,row[8],4097)
            domain=1-(1+row[2])*row[1]**2*radii**2
            domains.append(float(domain.min()))
            if not np.all(domain>0):issues.append(f'surface{j+1}_conic_domain')
        for j in range(5):
            r=np.linspace(0,min(surfaces[j,8],surfaces[j+1,8]),4097)
            gap=surfaces[j+1,0]+direct_sag(surfaces[j+1],r)-surfaces[j,0]-direct_sag(surfaces[j],r)
            gaps.append(float(np.min(gap)))
            if not np.all(np.isfinite(gap)) or not np.all(gap>0):issues.append(f'physical_adjacent_surface_crossing{j+1}')
        # Fresh full13field/RGB traces; failed rays are retained and counted.
        t=trace_design(d,fields=np.linspace(0,1,13),rings=16,arms=48)
        valid=(t['status']==0)&t['image_forward']&np.isfinite(t['image_xy']).all(axis=-1)
        if not np.all(valid):issues.append('failed_ray_in_fresh13field_RGB_trace')
        if valid.size!=31863:issues.append('wrong_fresh_dense_ray_count')
        rim=np.isclose(np.linalg.norm(t['pupil_xy'],axis=1),d['stop_radius_mm'],rtol=1e-11,atol=1e-9)
        axis_rim=np.sqrt(np.sum(t['dir'][:,0,rim,:2]**2,axis=-1))
        nas=axis_rim.min(axis=1)
        if not np.all(valid[:,0,rim]) or not np.all(nas>=.7):issues.append('RGB_complete_rim_NA_below0p7_or_failed')
        heights=t['image_xy'][:,-1,0,1]
        if abs(heights[1]-6.4)>.001:issues.append('588nm_fullfield_height_error_over0p001mm')
        epd=2*d['stop_radius_mm'];front=2*d['clear_radii_mm'][0]
        if not (abs(epd-50)<1e-8 or abs(front-50)<1e-8):issues.append('neither_aperture50_alternative')
        # Fresh metric calculation checks stale/exported metrics and independent
        # chief derivative/Fresnel paths, not stored booleans alone.
        fresh=evaluate_design(d,image_height_tolerance_mm=.001)
        if not fresh['all_declared_numeric_constraints_and_sampled_rays_passed']:issues.append('fresh_radiometry_constraints_fail')
        stored=c['metrics'];metric_error=0.
        for wl in ('486.0','588.0','656.0'):
            for field in ('0.0','0.5','0.8','1.0'):
                a=stored['by_wavelength_nm'][wl]['fields'][field];b=fresh['by_wavelength_nm'][wl]['fields'][field]
                for key in ('chief_height_mm','RI_geometric','RI_uncoated_fresnel','distortion_percent_paraxial_EFL','geometric_transmitted_pupil_fraction'):
                    if a[key] is None or b[key] is None:issues.append('undefined_metric_'+wl+'_'+field+'_'+key)
                    else:metric_error=max(metric_error,abs(a[key]-b[key]))
                metric_error=max(metric_error,abs(a['spot_geometric']['rms_radius_mm']-b['spot_geometric']['rms_radius_mm']))
        if metric_error>1e-9:issues.append('stored_metrics_do_not_match_fresh_design')
        proxy=c['geometric_mtf_proxy']
        if not proxy['all_rays_valid']:issues.append('geometric_score_contains_failed_rays')
        f=np.asarray(proxy['frequencies_lp_mm']);curves=np.asarray([*proxy['MTF_t'],*proxy['MTF_s']])
        auc=np.trapezoid(curves,f,axis=-1)/20
        if np.max(np.abs([auc.mean()-proxy['mean_AUC_normalized'],auc.min()-proxy['minimum_AUC_normalized']]))>1e-12:
            issues.append('stored_AUC_not_consistent_with_curves')
        ri=fresh['polychromatic_fields']['1.0']['RI_uncoated_fresnel']
        score=.6*auc.mean()+.4*auc.min()-.25*abs(ri-1)
        scores.append(score)
        if abs(score-c['selection']['score'])>1e-12:issues.append('selection_score_formula_mismatch')
        ap=c['aperture_fit']
        if ap['sampled_ray_count']!=31863 or ap['original_transmitted_ray_count']!=31863:
            issues.append('incomplete_original_aperture_fit_trace')
        if ap.get('result_design_sha256')!=sha:
            issues.append('aperture_fit_result_hash_does_not_match_design')
        if ap.get('accepted'):
            if ap.get('lost_transmitted_ray_count') or ap.get('changed_ray_status_count'):
                issues.append('accepted_aperture_fit_did_not_preserve_transmitted_rays')
        elif sha!=ap.get('original_design_sha256'):
            issues.append('rejected_aperture_fit_did_not_preserve_original_design')
        records.append({'candidate_id':c['candidate_id'],'source_seed_id':c['source_seed_id'],
                        'passed':not issues,'failures':issues,'design_sha256':sha,
                        'export_roundtrip':roundtrip,'fresh_dense_ray_count':int(valid.size),
                        'fresh_failed_ray_count':int(np.sum(~valid)),
                        'fresh_native_status_counts':{str(int(k)):int(v) for k,v in zip(*np.unique(t['status'],return_counts=True))},
                        'NA_rgb_axis_complete_rim_minimum':nas.tolist(),
                        'fullfield_chief_height_RGB_mm':heights.tolist(),
                        'fullfield_chief_height_RGB_minus6p4_mm':(heights-6.4).tolist(),
                        'RGB_chief_beyond6p4_on_positive_field':(heights>6.4+1e-10).tolist(),
                        'EPD_mm':epd,'first_front_clear_diameter_mm':front,
                        'vertex_BFL_mm':float(distances[5]),
                        'paraxial_BFL588_mm':fresh['first_order_primary588nm']['paraxial_bfl_mm'],
                        'minimum_adjacent_surface_separations4097_mm':gaps,
                        'minimum_conic_domains':domains,
                        'maximum_stored_vs_fresh_metric_difference':metric_error,
                        'aperture_fit_accepted':bool(ap.get('accepted')),
                        'aperture_fit_applied':bool(ap.get('applied')),
                        'RI_RGB_weighted_edge_geometric':fresh['polychromatic_fields']['1.0']['RI_geometric'],
                        'RI_RGB_weighted_edge_uncoated':ri,
                        'independent_recomputed_selection_score':float(score)})
        print(c['candidate_id'],'PASS' if not issues else issues,flush=True)
    if any(scores[i]<scores[i+1]-1e-12 for i in range(len(scores)-1)):global_fail.append('rank_not_descending_score')
    return {'schema':'frozen-selection-independent-physics-audit/1','passed':not global_fail and all(r['passed'] for r in records),
            'candidate_count':len(candidates),'distinct_parent_count':len(set(parents)),
            'fresh_dense_ray_count':sum(r['fresh_dense_ray_count'] for r in records),
            'fresh_failed_ray_count':sum(r['fresh_failed_ray_count'] for r in records),
            'global_failures':global_fail,'candidates':records,
            'input_path':str(path),'input_file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'definitions':{'NA_acceptance':'onaxis_complete_rim_minimum_allRGB;notfullfieldNAguarantee',
                           'image_height':'588nm6.4mmwith0.001mmnumericproximity;RGBheightsreportednotforcedequal;nodetectorcropping',
                           'RI':'conditionalchiefmappinguniformradiance;uncoatedscalarFresneloptional;notexactfullPSFdetectorirradiance',
                           'geometry':'4097radialsamplesnotcontinuousmathematicalextremumproof',
                           'ranking':'declaredgeometricMTFAUC+RIheuristic;notglobaloptimumorconvergedwaveMTFrank'}}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=audit(args.input)
    with args.output.open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('candidates','definitions')},indent=2))
    if not report['passed']:raise SystemExit(1)
