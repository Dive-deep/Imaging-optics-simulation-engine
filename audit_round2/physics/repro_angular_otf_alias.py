"""Requested angular frequency can be far above resolved angular quadrature."""
import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'input'/'optical_engine'
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from sequential_optics.wave import afocal_fft,angular_otf
N=64;du=.01;lam=.00055
u=(np.arange(N)-N//2)*du;x,y=np.meshgrid(u,u)
# Gaussian with insignificant edge field; same input samples in every run.
pupil=np.exp(-(x*x+y*y)/.06**2)
f_alias=(N*4*du)/lam
rows=[]
for padding in (4,8,16):
 r=afocal_fft(pupil,du,lam,padding=padding)
 freq=np.array([0.,100.,1000.,f_alias,10000.])
 mtf=abs(angular_otf(r,freq))
 rows.append(dict(padding=padding,direction_step=float(r['dl']),frequencies_cycles_rad=freq.tolist(),mtf=mtf.tolist(),edge_energy_fraction=r['edge_energy_fraction']))
# Small-angle analytic Gaussian reference: PSF angular variance lambda²/(4*pi²*w²),
# therefore angular MTF = exp[-lambda² f²/(2*w²)]. At alias f, vanishingly small.
result=dict(boundary_grid=N,spacing_mm=du,wavelength_mm=lam,gaussian_waist_mm=.06,alias_frequency=f_alias,analytic_paraxial_mtf_at_alias=float(np.exp(-lam**2*f_alias**2/(2*.06**2))),rows=rows)
from sequential_optics.models import OpticalSystem,Surface,Field
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_bundle
from sequential_optics.wave import analyze_wave
from sequential_optics.settings import EngineSettings,SamplingSettings
system=OpticalSystem((Surface(0.,0.,.3),),{'air':Material('air',n_constant=1.)},image_z=None)
cfg=EngineSettings(sampling=SamplingSettings(pupil_rings=8,pupil_arms=32,wave_grid_size=64,fft_padding=4))
bundle=trace_bundle(system,Field(),settings=cfg)
# Explicit physical boundary E amplitude independent of geometrical ray weights.
amp=np.exp(-np.sum((.3*bundle.pupil_xy)**2,axis=1)/.06**2)
f_alias=4*(2*.3*1.05)/lam
default_axis=np.linspace(0,cfg.sampling.max_angular_frequency_cycles_rad,cfg.sampling.angular_frequency_samples)
nearest_default=float(default_axis[np.argmin(abs(default_axis-f_alias))])
rows=[]
for padding in (4,8,16):
 r=analyze_wave(system,bundle,cfg,mode='afocal',padding=padding,source_boundary_amplitudes=amp,frequencies=[0.,100.,f_alias,nearest_default])
 rows.append({'padding':padding,'frequencies':[0.,100.,f_alias,nearest_default],'mtf_s':r['mtf_s'].tolist(),'status':r['status'],'phase_max':r['sampling_diagnostics']['max_grid_phase_step_rad'],'edge_energy_fraction':r['edge_energy_fraction']})
result['full_analyze_wave_rows']=rows
print(json.dumps(result,indent=2))
