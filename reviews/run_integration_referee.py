"""Run from project root: PYTHONPATH=src python reviews/run_integration_referee.py.
Independent integration diagnostics; synthetic lens, no vendor oracle.
"""
import sys,json,time
from dataclasses import replace
from pathlib import Path
import numpy as np
sys.path.insert(0,'src');sys.path.insert(0,'.')
from tests.test_adversarial import synthetic_triplet
from sequential_optics.models import Field
from sequential_optics.settings import EngineSettings,SamplingSettings
from sequential_optics.geometry import trace_bundle,validate_geometry
from sequential_optics.metrics import spot_diagram,common_geometric_refocus,angular_spot,chief_mapping
from sequential_optics.wave import analyze_wave
from sequential_optics.thermal import thermal_system
from sequential_optics.radiometry import RectangularSource,far_field

out={'kind':'independent numerical integration, no vendor reference','cases':[]};t0=time.perf_counter()
import atexit
atexit.register(lambda: Path('reviews/integration_referee_results.json').write_text(json.dumps(out,indent=2)))
def settings(r=8,a=32,n=64):return EngineSettings(sampling=SamplingSettings(pupil_rings=r,pupil_arms=a,wave_grid_size=n))
def wavesum(w):
 d=w['sampling_diagnostics']
 return {'status':w['status'],'mtf_s':np.asarray(w.get('mtf_s',[])).tolist(),'mtf_t':np.asarray(w.get('mtf_t',[])).tolist(),'grid_phase_step':d['max_grid_phase_step_rad'],'triangle_phase_span':d['max_traced_triangle_phase_span_rad'],'mask_count':int(w['pupil_mask'].sum()),'edge_fraction':w.get('edge_energy_fraction'),'local':w.get('local_expansion_diagnostics')}
base=synthetic_triplet();validate_geometry(base)
for r,a,n in [(6,24,64),(10,36,64),(14,44,64),(14,44,80)]:
 cfg=settings(r,a,n);bs=[trace_bundle(base,Field(f,0),settings=cfg) for f in [0,.03]]
 common=common_geometric_refocus(base,bs)
 row={'type':'focal_sampling','rings':r,'arms':a,'grid':n,'fixed_rms_mm':[spot_diagram(base,b)['rms_radius'] for b in bs],'common_z':common['common_image_z_mm'],'common_rms_mm':[s['rms_radius'] for s in common['common_focus_spots']]}
 w=analyze_wave(base,bs[1],cfg,frequencies=[0,10,30,60,120])
 row['fixed_wave']=wavesum(w);row['common_wave']=wavesum(analyze_wave(base,bs[1],cfg,image_z=common['common_image_z_mm'],boundary_state=w['boundary_state'],frequencies=[0,10,30,60,120]))
 out['cases'].append(row)
 print('sample',r,a,n,'done',flush=True)
# Larger field intentionally clips; this is not a design acceptance target.
for f in [.15,.45,.6]:
 cfg=settings(8,32,64);b=trace_bundle(base,Field(f,0),settings=cfg)
 row={'type':'high_field','field_rad':f,'status_counts':{s:sum(t.status==s for t in b.traces) for s in {t.status for t in b.traces}},'chief_status':b.chief.status}
 try:row['spot_rms_mm']=spot_diagram(base,b)['rms_radius'];row['wave']=wavesum(analyze_wave(base,b,cfg,frequencies=[0,10,30,60,120]))
 except Exception as e:row['analysis_exception']={'type':type(e).__name__,'message':str(e)}
 out['cases'].append(row);print('field',f,'done',flush=True)
# All coefficients contribute a controlled three micron at rim each.
a=6.;coeff=tuple((-1)**(m+1)*.003/a**(2*m) for m in range(1,7))
surf=list(base.surfaces);surf[4]=replace(surf[4],surface_type='even_asphere',coefficients=coeff)
asph=replace(base,surfaces=tuple(surf));validate_geometry(asph)
b=trace_bundle(asph,Field(.03,0),settings=settings())
out['cases'].append({'type':'all_asphere_terms','coefficients':coeff,'trace_success':int(b.success_mask.sum()),'trace_total':len(b.traces),'wave':wavesum(analyze_wave(asph,b,settings(),frequencies=[0,10,30,60,120]))})
# Finite source camera complement, local angular lens metric.
coll=replace(base,object_z=-78.52124279501913,image_z=None)
for h in [0.,3.2,6.4]:
 b=trace_bundle(coll,Field(h,0,'height'),settings=settings(10,36,80))
 row={'type':'finite_collimator','source_x_mm':h,'trace_success':int(b.success_mask.sum()),'angular_rms_rad':angular_spot(b)['rms_radius'],'chief_direction':b.chief.last_direction.tolist()}
 row['wave']=wavesum(analyze_wave(coll,b,settings(10,36,80),mode='afocal',frequencies=[0,50,100,200]))
 out['cases'].append(row);print('collimator',h,'done',flush=True)
# Explicit synthetic thermal laws; nominal materials had none.
materials={key:replace(m,index_temperature_invariant=True) if key=='air' else replace(m,cte_per_k=70e-6 if key=='positive' else 8e-6,dn_dT_per_k=-1e-4 if key=='positive' else 3e-6) for key,m in base.materials.items()}
thermal=replace(base,materials=materials)
for temp in [20.,80.,0.,20.]:
 s=thermal_system(thermal,temp);cfg=settings(10,36,80);bs=[trace_bundle(s,Field(f,0),settings=cfg) for f in [0,.03]]
 common=common_geometric_refocus(s,bs);w=analyze_wave(s,bs[1],cfg,frequencies=[0,10,30,60,120]);wr=analyze_wave(s,bs[1],cfg,image_z=common['common_image_z_mm'],boundary_state=w['boundary_state'],frequencies=[0,10,30,60,120])
 out['cases'].append({'type':'thermal','temperature_c':temp,'fixed_image_z':s.image_z,'common_z':common['common_image_z_mm'],'fixed_rms_mm':[spot_diagram(s,b)['rms_radius'] for b in bs],'common_rms_mm':[t['rms_radius'] for t in common['common_focus_spots']],'chief_map':None,'fixed_wave':wavesum(w),'common_wave':wavesum(wr)})
 print('temp',temp,'done',flush=True)
# Focused source importance integrates finite source transport, not MTF pixel fill.
ff=far_field(coll,RectangularSource(coll.object_z),np.linspace(0,.2,9),np.linspace(-np.pi,np.pi,17),sample_count=1024,launch_plane_z_mm=-.5,launch_radius_mm=6.6)
out['cases'].append({'type':'finite_source_far_field','source_power':ff['source_hemisphere_power'],'exit_power':ff['exit_power'],'efficiency':ff['direct_efficiency'],'outside_roi':ff['outside_reported_angular_roi_power'],'status_counts':ff['trace_status_counts'],'assumptions':ff['assumptions']})
# Independent propagation checks share the traced boundary (not a new oracle).
cfg=settings(14,44,80);b=trace_bundle(base,Field(.03,0),settings=cfg)
w=analyze_wave(base,b,cfg,frequencies=[0,10,30,60,120],direct_check_points=[[0,0],[.005,0],[.02,0],[.04,0]])
d=w['direct_scalar_check'];exact=d['exact_intensity'];approx=d['local_fft_integral_intensity']
out['cases'].append({'type':'near_focus_direct_RS','points_xy_mm':d['points_xy_mm'].tolist(),'exact_intensity':exact.tolist(),'local_intensity':approx.tolist(),'relative_intensity_difference':((approx-exact)/exact).tolist(),'shared_boundary':True})
surfaces=list(base.surfaces);surfaces[0]=replace(surfaces[0],curvature=0,surface_type='even_asphere',coefficients=(.01,1e-6,-2e-9,3e-12,-5e-15,7e-18))
pure=replace(base,surfaces=tuple(surfaces));validate_geometry(pure)
b=trace_bundle(pure,Field(.03,0),settings=cfg);w=analyze_wave(pure,b,cfg,frequencies=[0,10,30,60,120])
out['cases'].append({'type':'pure_asphere_c0','first_vertex_curvature':.02,'success':int(b.success_mask.sum()),'total':len(b.traces),'spot_rms_mm':spot_diagram(pure,b)['rms_radius'],'wave_status':w['status'],'mtf_s':w['mtf_s'].tolist(),'mtf_t':w['mtf_t'].tolist()})
for c in out['cases']:
 if c['type']=='thermal':
  v=chief_mapping(thermal_system(thermal,c['temperature_c']),[Field(0,0),Field(.03,0)])
  c['chief_map']=v['coordinates'].tolist()
out['seconds']=time.perf_counter()-t0
Path('reviews/integration_referee_results.json').write_text(json.dumps(out,indent=2))
print('FINISHED',out['seconds'],flush=True)
