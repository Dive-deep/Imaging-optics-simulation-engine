"""End-to-end tests use actual public APIs, not copies of implementation equations."""
import unittest
import sys
from pathlib import Path
from dataclasses import replace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples'))
from sample_systems import synthetic_triplet, synthetic_collimator
from sequential_optics import EngineSettings, SamplingSettings, Field
from sequential_optics.api import RunRequest, evaluate
from sequential_optics.io import system_to_dict, system_from_dict, json_value
from sequential_optics.geometry import trace_bundle
from sequential_optics.paraxial import paraxial_properties
from sequential_optics.metrics import spot_diagram

class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.settings=EngineSettings(sampling=SamplingSettings(pupil_rings=5,pupil_arms=20,wave_grid_size=32,fft_padding=4))

    def test_triplet_independent_abcd(self):
        system=synthetic_triplet()
        result=paraxial_properties(system,settings=self.settings)
        self.assertAlmostEqual(result['efl_image_mm'],84.1490162331,places=8)
        self.assertAlmostEqual(result['bfl_mm'],79.9243093735,places=8)
        bundle=trace_bundle(system,Field(),settings=self.settings)
        self.assertTrue(bundle.success_mask.all())
        self.assertEqual(bundle.traces[0].positions.shape,(6,3))
        spot=spot_diagram(system,bundle)
        # Symmetry oracle is independent of geometric spot's absolute size.
        points=np.array([t.last_position for t in bundle.traces])
        self.assertLess(abs(np.average(points[:,0],weights=bundle.weights)),1e-9)
        self.assertLess(abs(np.average(points[:,1],weights=bundle.weights)),1e-9)

    def test_nominal_g_no_thermal_no_rgb_requirements(self):
        system=synthetic_triplet()
        # Strip transmission too: nominal lens-only metrics must still work.
        mats={k:replace(v,absorption_constant_per_mm=None,absorption_temperature_invariant=False)
              for k,v in system.materials.items()}
        result=evaluate(replace(system,materials=mats),RunRequest(fields=(Field(),),
              metrics=('first_order','spot','mtf'),fail_fast=True),self.settings)
        wr=result['states'][0]['wavelengths'][0]
        self.assertEqual(wr['wavelength_nm'],550.)
        for key in ('spot','mtf'):self.assertEqual(wr['fields'][0][key]['status'],'ok')

    def test_optional_data_failure_does_not_poison_other_wavelength(self):
        result=evaluate(synthetic_triplet(),RunRequest(wavelengths=('G','R'),metrics=('first_order','spot')),self.settings)
        g,r=result['states'][0]['wavelengths']
        self.assertEqual(g['first_order']['status'],'ok')
        self.assertEqual(r['first_order']['status'],'missing_data')
        self.assertEqual(r['fields'][0]['rays']['status'],'missing_data')

    def test_nominal_then_hot_failure_is_local(self):
        result=evaluate(synthetic_triplet(),RunRequest(temperatures_c=(20.,60.),metrics=('first_order',)),self.settings)
        self.assertEqual(result['states'][0]['wavelengths'][0]['first_order']['status'],'ok')
        self.assertEqual(result['states'][1]['status'],'missing_data')

    def test_prescription_json_roundtrip(self):
        system=synthetic_triplet(aspheric=True)
        restored=system_from_dict(system_to_dict(system))
        b1=trace_bundle(system,Field(.01,0.),settings=self.settings)
        b2=trace_bundle(restored,Field(.01,0.),settings=self.settings)
        self.assertTrue(np.array_equal(b1.success_mask,b2.success_mask))
        np.testing.assert_allclose(b1.chief.positions,b2.chief.positions,rtol=0,atol=0)

    def test_afocal_request_keeps_finite_efl(self):
        result=evaluate(synthetic_collimator(),RunRequest(fields=(Field(kind='height'),),output_mode='afocal',
              metrics=('first_order','spot','common_focus'),fail_fast=True),self.settings)
        wr=result['states'][0]['wavelengths'][0]
        self.assertAlmostEqual(wr['first_order']['value']['efl_image_mm'],84.1490162331,places=8)
        self.assertEqual(wr['fields'][0]['spot']['status'],'ok')
        self.assertEqual(wr['common_focus']['status'],'not_applicable')

if __name__=='__main__':unittest.main()
