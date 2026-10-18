"""Field-specific fixed-axis FFT MTF; frequency direction is not bar direction."""
from dataclasses import replace
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
from sample_systems import synthetic_triplet, synthetic_collimator
from sequential_optics import EngineSettings, SamplingSettings, Field, RunRequest, evaluate
from sequential_optics.wave import fft_psf, sample_axis_otf
from sequential_optics.source_wave import source_diffraction_kernel


class FieldMtfAxesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = EngineSettings(sampling=SamplingSettings(
            pupil_rings=5, pupil_arms=24, wave_grid_size=32, fft_padding=4,
            max_frequency_lp_mm=60., frequency_samples=7,
            max_angular_frequency_cycles_rad=10., angular_frequency_samples=5))
        cls.fields = (Field(), Field(0., .08), Field(.08, 0.), Field(.03, .03))
        output = evaluate(synthetic_triplet(), RunRequest(fields=cls.fields,
            metrics=('mtf',), include_ray_history=False, fail_fast=True), cls.cfg)
        cls.spectral = output['states'][0]['wavelengths'][0]
        cls.results = [r['mtf']['value'] for r in cls.spectral['fields']]

    def test_each_requested_field_has_both_curves_and_reference(self):
        self.assertEqual(len(self.results), len(self.fields))
        self.assertEqual(self.spectral['wavelength_nm'], self.cfg.wavelength())
        for expected, outer, wave in zip(self.fields, self.spectral['fields'], self.results):
            self.assertEqual(outer['field'], wave['field'])
            self.assertEqual(wave['field'], {'x':expected.x,'y':expected.y,'kind':expected.kind})
            for key in ('mtf_t','mtf_s','otf_t','otf_s','diffraction_mtf_t','diffraction_mtf_s'):
                self.assertEqual(wave[key].shape, wave['frequencies_lp_mm'].shape)
                self.assertAlmostEqual(wave[key][0], 1., places=12)
            self.assertIn(outer['mtf']['analysis_status'], ('conditional','underresolved'))
            self.assertEqual(wave['vendor_equivalence'], 'unverified')

    def test_frequency_axes_and_target_lines_are_perpendicular(self):
        for wave in self.results:
            self.assertEqual(wave['axis_convention_id'], 'fixed_chief_normal_xy')
            self.assertFalse(wave['field_radial_axes'])
            e1, e2, chief = wave['output_basis']
            expected_s = np.array([1.,0.,0.])-chief[0]*chief
            expected_s /= np.linalg.norm(expected_s)
            np.testing.assert_allclose(wave['mtf_frequency_axes_global']['s'], expected_s, atol=1e-12)
            np.testing.assert_allclose(wave['mtf_frequency_axes_global']['t'], np.cross(chief,expected_s), atol=1e-12)
            for label in ('s','t'):
                frequency = wave['mtf_frequency_axes_global'][label]
                bars = wave['mtf_target_line_axes_global'][label]
                self.assertAlmostEqual(float(frequency@bars), 0., places=12)
                self.assertAlmostEqual(float(frequency@chief), 0., places=12)
                self.assertAlmostEqual(float(np.linalg.norm(frequency)), 1., places=12)

    def test_rotating_y_field_to_x_exchanges_fixed_axis_mtf(self):
        # Independent symmetry control for the actual round lens prescription.
        # The inherited phase status is retained; this is not a convergence test.
        y, x = self.results[1:3]
        self.assertGreater(np.max(abs(y['mtf_s']-y['mtf_t'])), .1)
        np.testing.assert_allclose(y['mtf_t'], x['mtf_s'], atol=2e-9, rtol=0)
        np.testing.assert_allclose(y['mtf_s'], x['mtf_t'], atol=2e-9, rtol=0)
        np.testing.assert_allclose(y['diffraction_mtf_t'], x['diffraction_mtf_s'], atol=2e-9, rtol=0)

    def test_anisotropic_pupil_frequency_axes_from_exact_rectangle_overlap(self):
        # Uniform rectangle autocorrelation is triangular in each direction.
        pupil = np.zeros((32,32), complex)
        pupil[12:20,8:24] = 1.
        step = .002; vacuum_mm = .00055
        wave = fft_psf(pupil, step, vacuum_mm, padding=4)
        frequency = np.array([0.,1.,3.,5.])*step/vacuum_mm
        np.testing.assert_allclose(sample_axis_otf(wave,frequency,'x'), 1-np.array([0.,1.,3.,5.])/16, atol=1e-14)
        np.testing.assert_allclose(sample_axis_otf(wave,frequency,'y'), 1-np.array([0.,1.,3.,5.])/8, atol=1e-14)

    def test_afocal_each_field_has_separate_angular_t_s(self):
        fields = (Field(kind='height'), Field(.1,0.,'height'))
        result = evaluate(synthetic_collimator(), RunRequest(fields=fields, metrics=('mtf',),
            output_mode='afocal', include_ray_history=False, fail_fast=True), self.cfg)
        for record in result['states'][0]['wavelengths'][0]['fields']:
            wave = record['mtf']['value']
            self.assertEqual(wave['axis_convention_id'], 'fixed_chief_normal_xy')
            self.assertIn('frequencies_cycles_rad', wave)
            self.assertNotIn('frequencies_lp_mm', wave)
            for label in ('s','t'):
                self.assertEqual(wave['mtf_'+label].shape, (5,))
                self.assertAlmostEqual(wave['mtf_'+label][0], 1., places=12)

    def test_physical_source_uses_same_explicit_axis_metadata(self):
        system = replace(synthetic_triplet(), object_z=-100.)
        source = source_diffraction_kernel(system,[0.,0.],.01,settings=self.cfg,
            transmission_mode='lossless_scalar',frequencies=[0.])
        wave = source['wave']
        self.assertEqual(wave['axis_convention_id'], 'fixed_chief_normal_xy')
        np.testing.assert_allclose(wave['mtf_frequency_axes_global']['s'], [1.,0.,0.], atol=1e-12)
        self.assertAlmostEqual(wave['mtf_t'][0], 1.)
        self.assertAlmostEqual(wave['mtf_s'][0], 1.)


if __name__ == '__main__':
    unittest.main()
