"""No mutations to release: stale real-ray bundle silently accepted by wave APIs."""
import sys,json
from pathlib import Path
from dataclasses import replace
ROOT=Path(__file__).resolve().parents[1]/'input'/'optical_engine'
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'examples')]
import numpy as np
from sequential_optics.geometry import trace_bundle
from sequential_optics.models import Field
from sequential_optics.wave import analyze_wave
from sequential_optics.source_wave import source_diffraction_kernel
from sequential_optics.settings import EngineSettings,SamplingSettings
from sample_systems import synthetic_triplet
cfg=EngineSettings(sampling=SamplingSettings(pupil_rings=6,pupil_arms=24,wave_grid_size=64,fft_padding=4))
s=synthetic_triplet();f=Field();bundle=trace_bundle(s,f,settings=cfg)
surfaces=list(s.surfaces);surfaces[-1]=replace(surfaces[-1],curvature=surfaces[-1].curvature*.995)
changed=replace(s,surfaces=tuple(surfaces))
frequencies=[0.,10.,30.,60.,100.]
old=analyze_wave(s,bundle,cfg,frequencies=frequencies)
stale=analyze_wave(changed,bundle,cfg,frequencies=frequencies)
fresh=analyze_wave(changed,trace_bundle(changed,f,settings=cfg),cfg,frequencies=frequencies)
result={'lens_wave':{'frequencies_lp_mm':frequencies,'nominal_mtf_s':old['mtf_s'].tolist(),'stale_changed_system_mtf_s':stale['mtf_s'].tolist(),'retrace_changed_system_mtf_s':fresh['mtf_s'].tolist(),'stale_status':stale['status'],'fresh_status':fresh['status'],'stale_exactly_matches_nominal':bool(np.array_equal(stale['mtf_s'],old['mtf_s']))}}
# Same issue on physical-source bridge; first surface is untouched.
finite=replace(s,object_z=-200.,image_z=160.)
changed=replace(changed,object_z=-200.,image_z=160.)
b=trace_bundle(finite,Field(kind='height'),settings=cfg)
r=source_diffraction_kernel(changed,[0.,0.],.001,1.,bundle=b,settings=cfg,frequencies=frequencies)
result['source_wave']={'stale_bundle_accepted':True,'status':r['status'],'source_cell_collected_power':r['source_cell_collected_power']}
# A cached state also fails to detect changing STOP radius itself.
small_stop=replace(s,stop_radius=2.)
cached=analyze_wave(small_stop,bundle,cfg,boundary_state=old['boundary_state'],frequencies=frequencies)
fresh_small=analyze_wave(small_stop,trace_bundle(small_stop,f,settings=cfg),cfg,frequencies=frequencies)
result['cached_stop_radius']={'old_radius':4.,'new_radius':2.,'cached_state_accepted':True,'cached_mtf_s':cached['mtf_s'].tolist(),'fresh_mtf_s':fresh_small['mtf_s'].tolist(),'cached_status':cached['status'],'fresh_status':fresh_small['status'],'fingerprints_identical':cached['boundary_state'].fingerprint==old['boundary_state'].fingerprint}
print(json.dumps(result,indent=2))
