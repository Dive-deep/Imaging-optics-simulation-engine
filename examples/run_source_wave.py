"""Matched physical-source diffraction example; separate from lens-only MTF."""
from pathlib import Path
from dataclasses import replace
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sample_systems import synthetic_triplet,synthetic_collimator
from sequential_optics import EngineSettings,SamplingSettings
from sequential_optics.paraxial import paraxial_properties
from sequential_optics.source_wave import source_diffraction_kernel,bin_source_angular_response
from sequential_optics.imaging import PixelGrid,synthesize_image
from sequential_optics.io import save_results

def main():
    cfg=EngineSettings(sampling=SamplingSettings(pupil_rings=6,pupil_arms=32,wave_grid_size=48,fft_padding=2))
    camera=replace(synthetic_triplet(),object_z=-200.)
    camera=replace(camera,image_z=float(paraxial_properties(camera,settings=cfg)['paraxial_image_z']))
    detector=PixelGrid(np.linspace(-8.,8.,65),np.linspace(-2.,2.,25))
    points=[(-3.2,0.),(0.,0.),(3.2,0.)]
    # Three quadrature cells are a sample scene, not the full continuous panel.
    area=.04
    kernels=[];powers=[];records=[]
    for point in points:
        result=source_diffraction_kernel(camera,point,area,settings=cfg,
                fill_factor=.64,image_grid=detector,allow_local_detector_projection=True)
        kernels.append(result['kernel']);powers.append(result['collected_power_density'])
        records.append({k:result[k] for k in ('source_xy_mm','status','collected_power_density',
                       'raw_power_diagnostics','detector_projection','assumptions')})
    image=synthesize_image(kernels,powers,[area]*len(points),bitmap_weights=[1.,.5,1.])
    angular=source_diffraction_kernel(synthetic_collimator(),(0.,0.),.001,settings=cfg,fill_factor=.64,
             mode='afocal',frequencies=np.linspace(0,cfg.sampling.max_angular_frequency_cycles_rad,
                                                  cfg.sampling.angular_frequency_samples))
    binned=bin_source_angular_response(angular,np.linspace(0,.01,21),np.linspace(-np.pi,np.pi,33))
    out={'scope':'Three finite source quadrature cells and one collimator source cell; synthetic prescription, scalar approximation',
         'focal_cells':records,'focal_image':image,'afocal':binned,
         'afocal_raw_power':angular['raw_power_diagnostics'],
         'warning':'Off-axis focal detector uses explicitly selected local affine projection, not exact tilted-detector diffraction'}
    target=ROOT/'results/source_wave_demo.json';save_results(out,target);print(target)

if __name__=='__main__':main()
