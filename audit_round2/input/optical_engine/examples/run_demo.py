"""Reproducible synthetic optical scenarios; writes results, never asserts Zemax agreement."""
from pathlib import Path
import sys
import argparse
import time
from dataclasses import replace, asdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sample_systems import synthetic_triplet,synthetic_collimator
from sequential_optics import EngineSettings,SamplingSettings,Field
from sequential_optics.geometry import trace_bundle,sag_bounds
from sequential_optics.paraxial import paraxial_properties
from sequential_optics.metrics import spot_diagram,common_geometric_refocus,distortion,grid_mapping,chief_mapping
from sequential_optics.wave import analyze_wave
from sequential_optics.thermal import thermal_system
from sequential_optics.radiometry import RectangularSource,far_field,uniform_scene_ri
from sequential_optics.imaging import PixelGrid,geometric_image
from sequential_optics.io import save_system,save_results

def wave_summary(w):
    f=np.array([0,10,30,60,100,120.])
    return {'status':w['status'],'frequency_lp_mm':f,'mtf_x':np.interp(f,w['frequencies_lp_mm'],w['mtf_s']),
      'mtf_y':np.interp(f,w['frequencies_lp_mm'],w['mtf_t']),
      'diffraction_mtf_x':np.interp(f,w['frequencies_lp_mm'],w['diffraction_mtf_s']),
      'sampling_diagnostics':w['sampling_diagnostics'],'local_expansion_diagnostics':w['local_expansion_diagnostics']}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--settings',default=str(ROOT/'examples/settings.json'))
    p.add_argument('--output',default=str(ROOT/'results/demo'))
    p.add_argument('--source-samples',type=int,default=2048)
    p.add_argument('--image-rays-per-cell',type=int,default=64)
    p.add_argument('--ri-grid',type=int,default=32)
    a=p.parse_args();cfg=EngineSettings.from_json(a.settings)
    destination=Path(a.output);destination.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();system=synthetic_triplet(thermal=True,spectral=True)
    save_system(system,ROOT/'examples/triplet.json');save_system(synthetic_triplet(),ROOT/'examples/triplet_G_only.json')
    fields=[Field(0,0),Field(0,.015),Field(0,.03)]
    report={'scope':'Synthetic unoptimized lenses; actual engine execution, not vendor/measurement correlation',
            'settings':asdict(cfg),'source_discretization':vars(a),'thermal':[]}
    baseline_map=None
    for temperature in (20.,60.,85.):
        current=thermal_system(system,temperature,layout='fixed_anchors',anchor_fraction=.5,
                               stop_policy='fixed_radius',image_policy='fixed')
        bundles=[trace_bundle(current,f,settings=cfg) for f in fields]
        focus=common_geometric_refocus(current,bundles)
        waves=[analyze_wave(current,b,settings=cfg) for b in bundles]
        mapping=chief_mapping(current,fields,settings=cfg)['coordinates']
        if baseline_map is None:baseline_map=mapping.copy()
        report['thermal'].append({'temperature_c':temperature,'first_order':paraxial_properties(current,settings=cfg),
                 'spot_rms_mm':[spot_diagram(current,b)['rms_radius'] for b in bundles],
                 'wave':[wave_summary(w) for w in waves],'common_focus':focus,
                 'mapping_mm':mapping,'nominal_mapping_change_mm':mapping-baseline_map})
        for i,w in enumerate(waves):
            np.savez_compressed(destination/f'wave_T{temperature:g}_field{i}.npz',
                 psf=w['psf'],x_mm=w['x_mm'],y_mm=w['y_mm'],otf=w['otf'],
                 frequency_lp_mm=w['frequencies_lp_mm'],mtf_x=w['mtf_s'],mtf_y=w['mtf_t'],
                 diffraction_mtf_x=w['diffraction_mtf_s'],pupil=w['pupil'],opd_mm=w['pupil_opd_mm'])
    report['distortion']=distortion(system,fields,settings=cfg)
    # RI truncation bound: any admitted ray must connect the detector to the
    # final clear disk, therefore |direction cosine| <= |transverse slope|.
    last=system.surfaces[-1];zmax=last.z+sag_bounds(last)[1]
    xy=baseline_map;bound=float(np.max(np.abs(xy))+last.clear_radius)/(system.image_z-zmax)
    cosine_bound=min(1.,bound*1.05)
    report['ri']=uniform_scene_ri(system,xy,settings=cfg,direction_grid_size=a.ri_grid,
                direction_cosine_bounds=(-cosine_bound,cosine_bound,-cosine_bound,cosine_bound),axis_index=0)
    report['ri']['domain_coverage_bound']={'direction_cosine_abs_bound':bound,
           'method':'entire last-surface disk/sag bound at each finite image point'}
    # Separate projector branch: finite emitter, angular output, finite EFL.
    projector=synthetic_collimator()
    source=RectangularSource(projector.object_z)
    report['far_field']=far_field(projector,source,np.linspace(0,.14,15),np.linspace(-np.pi,np.pi,33),
                  settings=cfg,sample_count=a.source_samples,seed=2718,
                  launch_plane_z_mm=-5.,launch_radius_mm=7.)
    afocal=analyze_wave(projector,trace_bundle(projector,Field(kind='height'),settings=cfg),settings=cfg,mode='afocal',
                  frequencies=np.linspace(0,cfg.sampling.max_angular_frequency_cycles_rad,cfg.sampling.angular_frequency_samples))
    report['afocal']={k:afocal[k] for k in ('status','frequencies_cycles_rad','mtf_s','mtf_t','sampling_diagnostics')}
    report['angular_grid']=grid_mapping(projector,[-6.4,0,6.4],[-1.6,0,1.6],angular=True,settings=cfg)
    # Finite-conjugate camera image. A coarse bitmap is quadrature of a
    # continuous source scene, not explicit tracing of all 25,600 LED pixels.
    camera=replace(synthetic_triplet(),object_z=-200.)
    camera=replace(camera,image_z=float(paraxial_properties(camera,settings=cfg)['paraxial_image_z']))
    bitmap=np.zeros((4,8));bitmap[:,1:3]=1.;bitmap[1:3,5:7]=.5
    detector=PixelGrid(np.linspace(-8,8,65),np.linspace(-3,3,25))
    rendered=geometric_image(camera,bitmap,RectangularSource(camera.object_z),detector,settings=cfg,
                samples_per_source=a.image_rays_per_cell,seed=31415,launch_plane_z_mm=-5.,launch_radius_mm=7.)
    report['image']=rendered
    np.save(destination/'source_bitmap.npy',bitmap)
    report['seconds']=time.perf_counter()-start
    save_results(report,destination/'summary.json')
    print(f'Synthetic demo written to {destination}; {report["seconds"]:.2f}s')

if __name__=='__main__':main()
