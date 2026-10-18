"""Read-only independent optical engineering audit. Outputs findings.json.

Use ENGINE_ROOT=/path/to/optical_engine if the reviewed release is elsewhere.
These are characterization assertions: confirming an issue is NOT physics validation.
"""
from pathlib import Path
import os, sys, json, math
from collections import Counter
from dataclasses import replace
import numpy as np
ROOT=Path(os.environ.get('ENGINE_ROOT',str(Path(__file__).resolve().parents[1]/'input'/'optical_engine')))
sys.path[:0]=[str(ROOT/'src'),str(ROOT)]
from sequential_optics.models import Surface,OpticalSystem,Field,Ray
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_bundle,trace_ray,launch_plane_z,refract,sag
from sequential_optics.metrics import spot_diagram,longitudinal_spherical_aberration,ray_plane_intercept
from sequential_optics.paraxial import paraxial_properties
from sequential_optics.thermal import thermal_system
from examples.sample_systems import synthetic_triplet

def spot_outcome(s,b):
 try:
  r=spot_diagram(s,b)
  return {'status':r['status'],'rms_radius_mm':r['rms_radius'],'transmitted_quadrature_fraction':r['transmitted_quadrature_fraction']}
 except Exception as e:return {'status':'exception','type':type(e).__name__,'message':str(e)}

def counts(b):return dict(Counter(t.status for t in b.traces))

air=Material('air',n_constant=1,index_temperature_invariant=True,cte_per_k=0)
s=OpticalSystem((Surface(0,0,2),Surface(10,0,2)),{'air':air},stop_index=1,image_z=20)
wide=replace(s,surfaces=(replace(s.surfaces[0],mechanical_radius=4),s.surfaces[1]))
f=Field(.1,0)
b=trace_bundle(s,f,rings=8,arms=24,include_boundary=False)
bw=trace_bundle(wide,f,rings=8,arms=24,include_boundary=False)
d=10*math.tan(.1);R=2
expected_mask=np.linalg.norm(2*b.pupil_xy-np.array([d,0]),axis=1)<=2
expected_fraction=float(b.weights[expected_mask].sum())
expected_xy=2*b.pupil_xy[expected_mask]+np.array([d,0])
w=b.weights[expected_mask];cen=(w/w.sum())@expected_xy
expected_rms=float(np.sqrt((w/w.sum())@np.sum((expected_xy-cen)**2,axis=1)))
exact_fraction=(2*R*R*math.acos(d/(2*R))-.5*d*math.sqrt(4*R*R-d*d))/(math.pi*R*R)
result={'engine_root':str(ROOT),'LD1_pre_stop_vignetting':{
 'prescription':'air planes (z=0,r=2),(z=10,r=2 STOP), IMA20; infinite angle(.1,0)rad',
 'rim_samples_disabled_to_isolate_cause':True,'actual_counts':counts(b),'actual_spot':spot_outcome(s,b),
 'expected_sample_count':int(expected_mask.sum()),'expected_sample_fraction':expected_fraction,'expected_sample_rms_mm':expected_rms,
 'continuous_disk_overlap_fraction':exact_fraction,
 'control_only_front_mechanical_radius_to_4':{'counts':counts(bw),'spot':spot_outcome(wide,bw)}}}
assert result['LD1_pre_stop_vignetting']['actual_spot']['status']=='exception'
assert np.array_equal(bw.success_mask,expected_mask)
assert abs(spot_diagram(wide,bw)['rms_radius']-expected_rms)<1e-13

bb=trace_bundle(wide,f,rings=8,arms=24)
failed=[{'pupil':p.tolist(),'weight':float(w),'status':t.status,'message':t.message} for p,w,t in zip(bb.pupil_xy,bb.weights,bb.traces) if t.status=='aim_unresolved']
pupil=np.array([math.cos(3*math.pi/2),math.sin(3*math.pi/2)])
target=2*pupil;z=launch_plane_z(wide);direction=np.array([math.tan(.1),0,1]);direction/=np.linalg.norm(direction)
origin=np.array([*(target-(10-z)*direction[:2]/direction[2]),z])
tr=trace_ray(wide,Ray(origin,direction),stop_at=0,ignore_apertures=True)
o,direction=tr.last_position,tr.last_direction
hit=o+(10-o[2])*direction/direction[2]
a=float(direction[:2]@direction[:2]);bcoef=2*float(o[:2]@direction[:2]);c=float(o[:2]@o[:2])-4
disc=bcoef*bcoef-4*a*c
result['LD2_false_boundary_miss']={'counts':counts(bb),'spot':spot_outcome(wide,bb),'unresolved_samples':failed,
 'analytic_plane_hit_mm':hit.tolist(),'analytic_hit_radius_mm':float(np.linalg.norm(hit[:2])),
 'plane_normal_dot_direction':float(direction[2]),'cylinder_discriminant':disc,'b_squared':bcoef*bcoef,'four_ac':4*a*c}
assert disc<0 and abs(np.linalg.norm(hit[:2])-2)<1e-13 and len(failed)==1 and failed[0]['weight']==0

p=synthetic_triplet();zf=paraxial_properties(p)['paraxial_image_z'];ref=zf+10
lsa=longitudinal_spherical_aberration(p,pupil_coordinates=[0,1e-3,.1],reference_focus_z=ref)
result['LD3_lsa_custom_reference']={'actual_paraxial_focus_mm':zf,'requested_reference_mm':ref,
 'reported_axis_crossings_mm':lsa['axis_crossing_z_mm'].tolist(),'reported_lsa_mm':lsa['lsa_mm'].tolist(),
 'expected_q_zero_axis_crossing_mm':zf,'expected_q_zero_lsa_mm':-10.}
assert lsa['lsa_mm'][0]==0 and abs(lsa['lsa_mm'][1]+10)<1e-6

rows=[]
for rings in (1,8):
 bb=trace_bundle(p,Field(),rings=rings,arms=24);out=spot_diagram(p,bb)
 xs=np.array([ray_plane_intercept(t,p.image_z)[0] for t in bb.traces if t.success])
 mx=float(np.linalg.norm(xs-out['centroid'],axis=1).max())
 rows.append({'rings':rings,'arms':24,'reported_geo_radius_sampled_mm':out['geo_radius_sampled'],
 'max_all_transmitted_sample_rays_mm':mx,'reported_coordinate_count':len(out['coordinates_mm']),
 'transmitted_trace_count':len(xs),'understatement_percent':100*(1-out['geo_radius_sampled']/mx)})
result['LD4_geo_radius_omits_traced_rim']=rows
assert all(r['max_all_transmitted_sample_rays_mm']>r['reported_geo_radius_sampled_mm'] for r in rows)

# Positive controls: an independent closed-form thick biconvex lens and uniform scaling.
glass=Material('glass',n_constant=1.5,index_temperature_invariant=True,cte_per_k=1e-4)
lens=OpticalSystem((Surface(0,1/50,5,medium_after='glass'),Surface(5,-1/50,5)),{'air':air,'glass':glass},image_z=60)
phi1=(1.5-1)/50;phi2=(1-1.5)/(-50);t_over_n=5/1.5
power=phi1+phi2-t_over_n*phi1*phi2
expected_efl=1/power;expected_bfl=(1-t_over_n*phi1)/power
actual=paraxial_properties(lens)
hot=thermal_system(lens,120);hotp=paraxial_properties(hot)
angle=.3;d0=np.array([math.sin(angle),0,math.cos(angle)]);r=refract(d0,[0,0,1],1,1.5)
scale=1.01
surface=Surface(0,.02,3,conic=-.5,coefficients=(1e-3,-2e-5,3e-7,1e-10,-1e-13,2e-16),surface_type='even_asphere')
from sequential_optics.thermal import scale_surface
radii=np.linspace(0,3,17)
shape_error=float(np.max(np.abs(sag(scale_surface(surface,scale),scale*radii,0)-scale*sag(surface,radii,0))))
result['positive_analytic_controls']={
 'thick_lens_efl_expected_mm':expected_efl,'thick_lens_efl_actual_mm':actual['efl_image_mm'],
 'thick_lens_bfl_expected_mm':expected_bfl,'thick_lens_bfl_actual_mm':actual['bfl_mm'],
 'uniform_thermal_efl_ratio':hotp['efl_image_mm']/actual['efl_image_mm'],'expected_thermal_ratio':scale,
 'snell_transmitted_sin_actual':float(r[0]),'snell_transmitted_sin_expected':math.sin(angle)/1.5,
 'all_A2_A12_thermal_sag_identity_error_mm':shape_error}
assert abs(actual['efl_image_mm']-expected_efl)<1e-12
assert abs(actual['bfl_mm']-expected_bfl)<1e-12
assert abs(hotp['efl_image_mm']/actual['efl_image_mm']-scale)<1e-13
assert abs(r[0]-math.sin(angle)/1.5)<1e-15
assert shape_error<1e-14
out=Path(__file__).with_name('findings.json');out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
