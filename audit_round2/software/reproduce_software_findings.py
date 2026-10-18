"""Independent review reproductions; never writes or imports the editable release.

Run: python reproduce_software_findings.py [--output report.json]
No patches are applied. Exit code 0 means the known release behavior reproduced.
"""
from pathlib import Path
import argparse
import json
import platform
import sys

RELEASE = Path(__file__).resolve().parents[1] / 'input' / 'optical_engine'
sys.path.insert(0, str(RELEASE / 'src'))
sys.dont_write_bytecode = True
import numpy as np
import scipy
import sequential_optics
from sequential_optics import EngineSettings, RaySettings, SamplingSettings, OpticalSystem, Surface, Material, Field, RunRequest, evaluate
from sequential_optics.geometry import trace_bundle
from sequential_optics.metrics import spot_diagram, common_geometric_refocus
from sequential_optics.imaging import PixelGrid, SourceKernel, synthesize_image
from sequential_optics.io import json_value


def specimen():
    return OpticalSystem(
        (Surface(0., .02, 2., medium_after='glass'), Surface(3., -.02, 2.)),
        {'air': Material('ideal air', n_constant=1., index_temperature_invariant=True),
         'glass': Material('G-only material', wavelengths_nm=(550.,), indices=(1.5,))},
        image_z=52.)


def settings_propagation():
    cfg = EngineSettings(ray=RaySettings(direction_tol=1e-15),
        sampling=SamplingSettings(boundary_clearance_mm=1e-12, pupil_rings=2,
                                  pupil_arms=8, include_pupil_boundary=False))
    system = OpticalSystem((Surface(0., 0., 1.),),
                          {'air': Material('ideal air', n_constant=1.)}, image_z=1e-3)
    field = Field(x=np.pi/2 - 1e-13)
    bundle = trace_bundle(system, field, settings=cfg)
    direct = spot_diagram(system, bundle, settings=cfg)
    direct_focus = common_geometric_refocus(system, (bundle,), settings=cfg)
    batch = evaluate(system, RunRequest(fields=(field,), metrics=('spot', 'common_focus'),
                                       include_ray_history=False), cfg)
    wavelength = batch['states'][0]['wavelengths'][0]
    result = {'configured_direction_tol': cfg.ray.direction_tol,
              'default_direction_tol': EngineSettings().ray.direction_tol,
              'chief_dz': float(bundle.chief.last_direction[2]),
              'successful_rays': int(bundle.success_mask.sum()), 'total_rays': len(bundle.traces),
              'direct_spot_status': direct['status'],
              'batch_spot': wavelength['fields'][0]['spot'],
              'direct_common_focus_status': direct_focus['status'],
              'batch_common_focus': wavelength['common_focus']}
    assert result['successful_rays'] == result['total_rays'] == 16
    assert result['direct_spot_status'] == 'ok'
    assert result['batch_spot']['status'] == 'not_applicable'
    assert result['batch_common_focus']['status'] == 'not_applicable'
    return result


def writable_image_values():
    grid = PixelGrid([0., 1.], [0., 1.])
    kernel = SourceKernel(grid, [[1.]], 550.)
    normal = synthesize_image([kernel], [1.], [1.])
    writable = {'grid_edges': grid.x_edges_mm.flags.writeable,
                'kernel_mass': kernel.probability_mass.flags.writeable}
    kernel.probability_mass[0, 0] = 2.
    invalid_mass = synthesize_image([kernel], [1.], [1.])
    grid.x_edges_mm[1] = -1.
    invalid_grid = synthesize_image([kernel], [1.], [1.])
    result = {'caller_in_place_mutation_required': True, 'writable': writable,
              'before_mutation': {'collected': normal['collected_power'], 'detected': normal['detected_roi_power']},
              'after_mass_mutation': {'collected': invalid_mass['collected_power'],
                  'detected': invalid_mass['detected_roi_power'], 'outside': invalid_mass['outside_roi_power']},
              'after_edge_mutation_irradiance': invalid_grid['irradiance'].tolist()}
    assert invalid_mass['detected_roi_power'] > invalid_mass['collected_power']
    assert invalid_grid['irradiance'][0, 0] < 0
    return result


def malformed_thermal_indices():
    request = RunRequest(metrics=('first_order',), temperatures_c=(20.,), fail_fast=False,
        thermal_options={'lens_groups': [{'front_index': 0.0, 'back_index': 1.0}]})
    try:
        evaluate(specimen(), request)
    except Exception as exc:
        assert type(exc) is TypeError
        return {'request_construction_accepted': True, 'fail_fast': False,
                'exception_type': type(exc).__name__, 'message': str(exc),
                'classification': 'malformed-input validation defect, not missing-data defect'}
    raise AssertionError('Expected current release raw TypeError')


def optional_data_control():
    cfg = EngineSettings(sampling=SamplingSettings(pupil_rings=2, pupil_arms=8))
    out = evaluate(specimen(), RunRequest(metrics=('first_order', 'spot'), wavelengths=('G', 'B'),
                                         temperatures_c=(None, 20., 40.), include_ray_history=False), cfg)
    statuses = []
    for state in out['states']:
        item = {'requested_temperature': state['requested_temperature_c'], 'state_status': state['status'], 'wavelengths': []}
        for w in state['wavelengths']:
            item['wavelengths'].append({'nm': w['wavelength_nm'], 'first_order': w['first_order']['status'],
                'spot_or_rays': w['fields'][0].get('spot', w['fields'][0].get('rays'))['status']})
        statuses.append(item)
    assert statuses[0]['wavelengths'][0]['first_order'] == 'ok'
    assert statuses[0]['wavelengths'][0]['spot_or_rays'] == 'ok'
    assert statuses[0]['wavelengths'][1]['first_order'] == 'missing_data'
    assert statuses[1]['wavelengths'][0]['first_order'] == 'ok'
    assert statuses[2]['state_status'] == 'missing_data'
    return statuses


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'review_role': 'independent_scientific_software_developer',
              'archive_sha256': 'd854cbaba525581088fe7b890ee54e2f1fb663c8ef8c6fb562b06a66a65284ef',
              'module_loaded_from': sequential_optics.__file__,
              'runtime': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
              'SW-01': settings_propagation(), 'SW-02': writable_image_values(),
              'SW-03': malformed_thermal_indices(), 'positive_optional_data_control': optional_data_control()}
    text = json.dumps(json_value(report), ensure_ascii=False, indent=2, allow_nan=False)
    if args.output: args.output.write_text(text + '\n', encoding='utf-8')
    print(text)


if __name__ == '__main__': main()
