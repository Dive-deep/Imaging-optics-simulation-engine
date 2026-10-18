"""Independent referee fixtures: analytic or metamorphic, not vendor truth.

These tests target prior review objections and real six-interface integration.
The fixture is deliberately synthetic and unoptimized; see reviews/sample_plan.md.
"""
import unittest
from dataclasses import replace
import numpy as np

from sequential_optics.errors import GeometryError, MissingDataError
from sequential_optics.geometry import aim_ray, intersect, trace_ray
from sequential_optics.materials import Material
from sequential_optics.models import Field, OpticalSystem, Ray, Surface
from sequential_optics.paraxial import paraxial_properties
from sequential_optics.settings import EngineSettings
from sequential_optics.thermal import scale_surface, thermal_system, validate_radial_separation


def synthetic_triplet(*, measured_wavelength=None):
    def m(name, n):
        if measured_wavelength is None:
            return Material(name, n_constant=n)
        return Material(name, wavelengths_nm=(measured_wavelength,), indices=(n,))
    materials = {'air': m('air', 1), 'positive': m('positive', 1.5),
                 'negative': m('negative', 1.62)}
    specs = [(0, 1/50, 6, 'positive'), (3, -1/50, 6, 'air'),
             (5, -1/40, 5, 'negative'), (7, 1/40, 5, 'air'),
             (9, 1/45, 6, 'positive'), (12, -1/45, 6, 'air')]
    surfaces = tuple(Surface(z=z, curvature=c, clear_radius=a, medium_after=medium)
                     for z, c, a, medium in specs)
    return OpticalSystem(surfaces, materials, image_z=91.92430937348901,
                         stop_index=2, stop_radius=4, name='referee_synthetic_triplet')


class IndependentTripletTests(unittest.TestCase):
    def test_exact_stop_rim_is_not_lost_to_looser_aim_tolerance(self):
        from sequential_optics.geometry import trace_bundle
        from sequential_optics.settings import SamplingSettings
        system = replace(synthetic_triplet(), object_z=-78.52124279501913, image_z=None)
        cfg = EngineSettings(sampling=SamplingSettings(pupil_rings=3, pupil_arms=36,
                                                      wave_grid_size=32))
        bundle = trace_bundle(system, Field(3.2, 0, 'height'), settings=cfg)
        rim = np.isclose(np.linalg.norm(bundle.pupil_xy, axis=1), 1., rtol=0, atol=1e-13)
        self.assertEqual(np.count_nonzero(rim), 36)
        self.assertTrue(all(t.success for t, selected in zip(bundle.traces, rim) if selected),
                        [(t.status, t.message) for t, selected in zip(bundle.traces, rim)
                         if selected and not t.success])
        self.assertTrue(np.all(bundle.weights[rim] == 0.))

    def test_actual_triplet_source_diffraction_global_bins_preserve_cell_power(self):
        from sequential_optics.source_wave import source_diffraction_kernel, bin_source_angular_response
        from sequential_optics.settings import SamplingSettings
        system = replace(synthetic_triplet(), object_z=-78.52124279501913, image_z=None)
        cfg = EngineSettings(sampling=SamplingSettings(pupil_rings=4, pupil_arms=16, wave_grid_size=32))
        result = source_diffraction_kernel(system, [0, 0], .001, 1., fill_factor=.64,
                                          settings=cfg, mode='afocal')
        bins = bin_source_angular_response(result, np.linspace(0, np.pi, 17),
                                          np.linspace(-np.pi, np.pi, 25))
        power = result['source_cell_collected_power']
        self.assertGreater(power, 0.)
        self.assertLess(power, np.pi*.64*.001)
        self.assertAlmostEqual(bins['bin_power'].sum()/power, 1., places=12)
        self.assertEqual(bins['registration'], 'full_unit_direction_rotation')
        # This closure verifies matched normalization/registration. The separate
        # raw-power diagnostic must exist; this test does not assert its accuracy
        # from the normalized Q*h identity.
        diagnostics = result['raw_power_diagnostics']
        self.assertGreater(diagnostics['wave_power_w_unvalidated'], 0.)
        self.assertIn(diagnostics['raw_wave_phase_sampling_status'],
                      ('necessary_check_passed', 'underresolved'))

    def test_nominal_lens_mtf_api_needs_only_requested_green_material_data(self):
        from sequential_optics.api import evaluate, RunRequest
        from sequential_optics.settings import SamplingSettings
        settings = EngineSettings(wavelengths_nm={'G': 532.}, sampling=SamplingSettings(
            pupil_rings=4, pupil_arms=16, wave_grid_size=32, frequency_samples=7))
        system = synthetic_triplet(measured_wavelength=532.)
        result = evaluate(system, RunRequest(metrics=('spot', 'mtf'), include_ray_history=False), settings)
        field = result['states'][0]['wavelengths'][0]['fields'][0]
        self.assertEqual(field['spot']['status'], 'ok')
        self.assertEqual(field['mtf']['status'], 'ok')
        value = field['mtf']['value']
        self.assertEqual(value['illumination'], 'lens_metric')
        self.assertAlmostEqual(value['wavelength_vacuum_mm'], .000532)
        self.assertEqual(len(value['mtf_s']), 7)
        # This says the requested dependency path executes; its coarse sampling
        # may correctly remain underresolved and is not a physical accuracy gate.
        self.assertIn(value['status'], ('conditional', 'underresolved'))

    def test_afocal_api_returns_angular_mtf_and_finite_vergence_reference(self):
        from sequential_optics.api import evaluate, RunRequest
        from sequential_optics.settings import SamplingSettings
        system = replace(synthetic_triplet(), object_z=-78.52124279501913, image_z=None)
        settings = EngineSettings(sampling=SamplingSettings(pupil_rings=4, pupil_arms=16,
            wave_grid_size=32, max_angular_frequency_cycles_rad=400., angular_frequency_samples=9))
        request = RunRequest(fields=(Field(.5, 0, 'height'),),
                             metrics=('mtf', 'field_curvature'), output_mode='afocal',
                             include_ray_history=False)
        field = evaluate(system, request, settings)['states'][0]['wavelengths'][0]['fields'][0]
        self.assertEqual(field['mtf']['status'], 'ok')
        self.assertEqual(len(field['mtf']['value']['mtf_s']), 9)
        self.assertEqual(field['mtf']['value']['frequencies_cycles_rad'][-1], 400.)
        self.assertEqual(field['field_curvature']['status'], 'ok')
        self.assertEqual(field['field_curvature']['value']['vergence_convention'],
                         'converging_positive_independent')

    def test_independent_first_order_constants(self):
        p = paraxial_properties(synthetic_triplet())
        self.assertAlmostEqual(p['efl_image_mm'], 84.14901623305201, places=10)
        self.assertAlmostEqual(p['bfl_mm'], 79.92430937348901, places=10)
        self.assertAlmostEqual(p['determinant'], 1, places=13)

    def test_finite_collimator_is_not_zero_power_system(self):
        system = replace(synthetic_triplet(), object_z=-78.52124279501913, image_z=None)
        p = paraxial_properties(system)
        self.assertEqual(p['status'], 'ok')
        self.assertIsNone(p['paraxial_image_z'])
        self.assertGreater(abs(p['matrix'][1, 0]), 0.01)
        self.assertAlmostEqual(p['finite_object_matrix'][1, 1], 0, places=13)
        # Actual finite-source rays still pass the six-interface trace.
        ray = aim_ray(system, Field(0.7, -0.2, 'height'), (0.2, 0.1))
        tr = trace_ray(system, ray)
        self.assertTrue(tr.success, tr.message)
        self.assertEqual(tr.positions.shape, (6, 3))

    def test_custom_green_only_data_has_no_thermal_or_rgb_dependency(self):
        settings = EngineSettings(primary_wavelength='G', wavelengths_nm={'G': 532.})
        system = synthetic_triplet(measured_wavelength=532.)
        ray = aim_ray(system, Field(0.03, 0), (0.2, -0.1), settings=settings)
        tr = trace_ray(system, ray, settings=settings)
        self.assertTrue(tr.success, tr.message)
        # This deliberately contains neither CTE, dn/dT nor absorption.
        baseline = thermal_system(system, system.reference_temperature_c)
        np.testing.assert_allclose([s.z for s in baseline.surfaces], [s.z for s in system.surfaces])
        with self.assertRaises(MissingDataError):
            trace_ray(system, ray, 650., settings=settings)

    def test_rotated_field_and_pupil_rotate_actual_six_surface_path(self):
        system = synthetic_triplet()
        field_a, pupil_a = Field(0.03, 0), np.array([0.3, -0.2])
        field_b, pupil_b = Field(0, 0.03), np.array([0.2, 0.3])
        a = trace_ray(system, aim_ray(system, field_a, pupil_a))
        b = trace_ray(system, aim_ray(system, field_b, pupil_b))
        self.assertTrue(a.success and b.success, (a.message, b.message))
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        np.testing.assert_allclose(b.positions, a.positions @ rotation.T, atol=2e-8, rtol=0)
        np.testing.assert_allclose(b.directions, a.directions @ rotation.T, atol=2e-10, rtol=0)
        np.testing.assert_allclose(b.opl, a.opl, atol=2e-8, rtol=0)

    def test_coordinate_translation_keeps_geometry_and_opd_differences(self):
        system = synthetic_triplet()
        delta = 137.4
        shifted = replace(system, surfaces=tuple(replace(s, z=s.z+delta) for s in system.surfaces),
                          image_z=system.image_z+delta)
        field = Field(0.02, -0.01)
        a, b = [], []
        for pupil in [(0, 0), (0.1, 0.2), (-0.3, 0.1)]:
            ta = trace_ray(system, aim_ray(system, field, pupil))
            tb = trace_ray(shifted, aim_ray(shifted, field, pupil))
            self.assertTrue(ta.success and tb.success, (ta.message, tb.message))
            np.testing.assert_allclose(tb.positions, ta.positions+[0, 0, delta], atol=3e-8, rtol=0)
            np.testing.assert_allclose(tb.directions, ta.directions, atol=3e-10, rtol=0)
            a.append(ta.last_opl); b.append(tb.last_opl)
        # Plane-wave eikonal may gain a common piston under origin translation.
        np.testing.assert_allclose(np.array(b)-b[0], np.array(a)-a[0], atol=3e-8, rtol=0)

    def test_full_geometry_scaling_and_wavelength_scaling_preserve_ray_directions(self):
        system = synthetic_triplet()
        factor = 2.3
        scaled = replace(system, surfaces=tuple(scale_surface(s, factor, s.z*factor) for s in system.surfaces),
                         image_z=system.image_z*factor, stop_radius=system.stop_radius*factor)
        field = Field(0.025, -0.015)
        base_settings = EngineSettings()
        scaled_settings = replace(base_settings,
                                  sampling=replace(base_settings.sampling,
                                      boundary_clearance_mm=base_settings.sampling.boundary_clearance_mm*factor))
        a = trace_ray(system, aim_ray(system, field, (0.25, 0.15)), 550.)
        b = trace_ray(scaled, aim_ray(scaled, field, (0.25, 0.15), 550.*factor, scaled_settings),
                      550.*factor, scaled_settings)
        self.assertTrue(a.success and b.success, (a.message, b.message))
        np.testing.assert_allclose(b.positions, factor*a.positions, atol=3e-8, rtol=0)
        np.testing.assert_allclose(b.directions, a.directions, atol=3e-10, rtol=0)
        np.testing.assert_allclose(b.opl/(550.*factor), a.opl/550., atol=1e-10, rtol=0)
        self.assertAlmostEqual(paraxial_properties(scaled)['efl_image_mm'],
                               factor*paraxial_properties(system)['efl_image_mm'], places=10)


class IndependentInvalidGeometryTests(unittest.TestCase):
    def test_paraboloid_two_positive_roots_chooses_first(self):
        surface = Surface(z=0, curvature=0, coefficients=(1.,), clear_radius=3,
                          surface_type='even_asphere')
        point, distance = intersect(surface, Ray(np.array([0., 0., -2.]), np.array([1., 0., 3.])))
        np.testing.assert_allclose(point, [1, 0, 1], atol=1e-9)
        self.assertAlmostEqual(distance, np.sqrt(10), places=9)

    def test_midradius_overlap_is_rejected_despite_positive_center_and_rim(self):
        # gap(r)=1−8r²+8r⁴ is positive at0 and1 but negative near1/√2.
        material = Material('synthetic', n_constant=1.5)
        system = OpticalSystem((Surface(0, 0, 1, medium_after='g'),
                                Surface(1, 0, 1, coefficients=(-8., 8.),
                                        surface_type='even_asphere', medium_after='air')),
                               {'air': Material('air', n_constant=1.), 'g': material})
        with self.assertRaises(GeometryError):
            validate_radial_separation(system)


if __name__ == '__main__':
    unittest.main()
