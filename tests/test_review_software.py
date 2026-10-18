"""Independent software-review regression cases (not implementation mirrors)."""
from dataclasses import replace
import json
import tempfile
from pathlib import Path
import unittest
import numpy as np

from sequential_optics.api import RunRequest, evaluate
from sequential_optics.errors import InputError
from sequential_optics.io import system_to_dict, system_from_dict, json_value
from sequential_optics.materials import Material
from sequential_optics.models import Surface, OpticalSystem, Field, Ray
from sequential_optics.settings import EngineSettings, RaySettings, SamplingSettings
from sequential_optics.paraxial import paraxial_properties


def specimen():
    return OpticalSystem((Surface(0., .02, 2., medium_after='glass'), Surface(5., -.02, 2.)),
                         {'air': Material('explicit ideal air', n_constant=1., index_temperature_invariant=True),
                          'glass': Material('single measured G', wavelengths_nm=(550.,), indices=(1.5,))}, image_z=55.)


class SoftwareReviewTests(unittest.TestCase):
    def test_g_only_first_order_does_not_require_loss_or_thermal(self):
        out=evaluate(specimen(), RunRequest(metrics=('first_order',)))
        self.assertEqual(out['states'][0]['wavelengths'][0]['first_order']['status'], 'ok')
        self.assertEqual(out['states'][0]['wavelengths'][0]['wavelength_nm'], 550.)

    def test_missing_blue_is_local_to_blue_analysis(self):
        out=evaluate(specimen(), RunRequest(metrics=('first_order',), wavelengths=('G','B')))
        self.assertEqual(out['states'][0]['wavelengths'][0]['first_order']['status'], 'ok')
        self.assertEqual(out['states'][0]['wavelengths'][1]['first_order']['status'], 'missing_data')

    def test_missing_thermal_is_local_to_hot_state(self):
        out=evaluate(specimen(), RunRequest(metrics=('first_order',), temperatures_c=(20.,40.)))
        self.assertEqual(out['states'][0]['wavelengths'][0]['first_order']['status'], 'ok')
        self.assertEqual(out['states'][1]['status'], 'missing_data')

    def test_non_default_green_settings_propagate(self):
        system=specimen()
        material=Material('custom G', wavelengths_nm=(532.,), indices=(1.5,))
        system=replace(system, materials={**system.materials, 'glass': material})
        out=evaluate(system, RunRequest(metrics=('first_order',)), EngineSettings(wavelengths_nm={'G':532.}))
        self.assertEqual(out['states'][0]['wavelengths'][0]['wavelength_nm'], 532.)
        self.assertEqual(out['states'][0]['wavelengths'][0]['first_order']['status'], 'ok')

    def test_noninteger_count_settings_are_rejected_before_execution(self):
        for kw in ({'ray': RaySettings(aim_max_iterations=1.5)},
                   {'ray': RaySettings(max_iterations=3.5)},
                   {'sampling': SamplingSettings(frequency_samples=4.5)},
                   {'sampling': SamplingSettings(pupil_rings=2.1)}):
            with self.subTest(kwargs=kw):
                with self.assertRaises(InputError): EngineSettings(**kw)

    def test_noninteger_stop_index_rejected_before_execution(self):
        with self.assertRaises(InputError): replace(specimen(),stop_index=.5)

    def test_ray_extreme_direction_does_not_silently_become_zero(self):
        try:
            ray=Ray(np.zeros(3),np.array([1e308,0.,1e308]))
        except InputError:
            return
        self.assertTrue(np.all(np.isfinite(ray.direction)))
        self.assertAlmostEqual(float(np.linalg.norm(ray.direction)),1.)

    def test_caller_mutation_cannot_change_settings_snapshot(self):
        mapping={'G':550.}
        settings=EngineSettings(wavelengths_nm=mapping)
        mapping['G']=532.
        self.assertEqual(settings.wavelength(),550.)

    def test_caller_mutation_cannot_replace_system_materials(self):
        system=specimen()
        original=dict(system.materials)
        system=replace(system,materials=original)
        original['glass']=Material('changed caller',n_constant=1.8)
        self.assertEqual(system.index('glass',550.),1.5)

    def test_json_thermal_lens_groups_roundtrip(self):
        system=specimen()
        system=replace(system,materials={**system.materials,'glass':Material('thermal G', wavelengths_nm=(550.,),indices=(1.5,),
                                                                         cte_per_k=1e-5,index_temperature_invariant=True)})
        options=json.loads('{"lens_groups":[{"front_index":0,"back_index":1,"anchor_fraction":0.5}]}')
        out=evaluate(system, RunRequest(metrics=('first_order',),temperatures_c=(40.,),thermal_options=options))
        self.assertEqual(out['states'][0]['wavelengths'][0]['first_order']['status'],'ok')

    def test_prescription_json_retains_single_sample_and_defaults(self):
        original=specimen()
        restored=system_from_dict(json.loads(json.dumps(system_to_dict(original))))
        self.assertEqual(original.surfaces,restored.surfaces)
        self.assertEqual(original.index('glass',550.),restored.index('glass',550.))
        a=evaluate(original,RunRequest(metrics=('first_order',)))
        b=evaluate(restored,RunRequest(metrics=('first_order',)))
        self.assertEqual(a['input_sha256'],b['input_sha256'])

    def test_result_complex_nonfinite_values_remain_explicit(self):
        converted=json_value({'otf':np.array([1+0j, .2+.3j]),'unknown':float('nan')})
        decoded=json.loads(json.dumps(converted,allow_nan=False))
        self.assertEqual(decoded['otf']['encoding'],'complex-array')
        self.assertEqual(decoded['otf']['imag'],[0.,.3])
        self.assertEqual(decoded['unknown']['encoding'],'nonfinite')

    def test_unused_material_name_cannot_change_report_efl(self):
        base=specimen()
        system=replace(base,ambient_material='vacuum',surfaces=(base.surfaces[0],replace(base.surfaces[1],medium_after='vacuum')),
                       materials={'vacuum':base.materials['air'],'glass':base.materials['glass']})
        augmented=replace(system,materials={**system.materials,'air':Material('unused air',n_constant=1.0003)})
        a=paraxial_properties(system)
        b=paraxial_properties(augmented)
        self.assertEqual(a['efl_report_mm'],b['efl_report_mm'])
        explicit=paraxial_properties(augmented,settings=EngineSettings(report_reference_material='air'))
        self.assertAlmostEqual(explicit['efl_report_mm']/a['efl_report_mm'],1.0003)
        self.assertEqual(explicit['efl_image_mm'],a['efl_image_mm'])

    def test_invalid_request_schema_rejected_at_boundary(self):
        with self.assertRaises(InputError): RunRequest(fields=({'x':0},))
        with self.assertRaises(InputError): RunRequest(thermal_options={'unrecognized_policy':'free'})


if __name__=='__main__': unittest.main()
