"""Run from project root to reproduce the physical-source triplet diagnostic."""
import json
import sys
from pathlib import Path
from dataclasses import replace
import numpy as np

sys.path.insert(0, 'src')
sys.path.insert(0, '.')
from tests.test_adversarial import synthetic_triplet
from sequential_optics.source_wave import source_diffraction_kernel, bin_source_angular_response
from sequential_optics.settings import EngineSettings, SamplingSettings
from sequential_optics.models import Field
from sequential_optics.geometry import trace_bundle

system = replace(synthetic_triplet(), object_z=-78.52124279501913, image_z=None)
settings = EngineSettings(sampling=SamplingSettings(pupil_rings=8, pupil_arms=32, wave_grid_size=48))
bundle = trace_bundle(system, Field(0, 0, 'height'), settings=settings)
result = source_diffraction_kernel(system, [0, 0], .001, 1., fill_factor=.64,
                                  settings=settings, bundle=bundle, mode='afocal')
bins = bin_source_angular_response(result, np.linspace(0, np.pi, 25), np.linspace(-np.pi, np.pi, 33))
record = {
    'case': 'actual synthetic refracting triplet finite-source afocal scalar diffraction',
    'status': result['status'],
    'source_cell_collected_power': result['source_cell_collected_power'],
    'collected_power_density': result['collected_power_density'],
    'bin_power_sum': float(bins['bin_power'].sum()),
    'outside_roi_power': bins['outside_roi_power'],
    'registration': bins['registration'],
    'raw_power_diagnostics': result['raw_power_diagnostics'],
    'sampling_diagnostics': result['wave']['sampling_diagnostics'],
}
Path('reviews/source_wave_camera_check.json').write_text(json.dumps(record, indent=2))
print(json.dumps(record, indent=2))
