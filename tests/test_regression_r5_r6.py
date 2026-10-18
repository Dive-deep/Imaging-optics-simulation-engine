"""Expected-behavior regressions for independent-review findings R5 and R6.

References are a plane-interface Snell calculation, weighted ray moments, and
fixed traced-sample maxima. No continuous-pupil maximum or vendor equivalence
is inferred from these controls.
"""
from types import SimpleNamespace
import unittest

import numpy as np
from numpy.testing import assert_allclose

from sequential_optics.errors import ModelApplicabilityError, NumericalError
from sequential_optics.geometry import trace_bundle
from sequential_optics.materials import Material
from sequential_optics.metrics import (
    angular_spot, common_geometric_refocus, longitudinal_spherical_aberration,
    ray_plane_intercept, spot_diagram,
)
from sequential_optics.models import Field, OpticalSystem, RayBundle, Surface, TraceResult
from sequential_optics.paraxial import paraxial_properties


def _trace(position=(0., 0., 0.), slope=(0., 0.)):
    direction = np.array([*slope, 1.], dtype=float)
    direction /= np.linalg.norm(direction)
    return TraceResult(positions=np.array([position], dtype=float),
                       directions=np.array([direction]), opl=np.zeros(1),
                       powers=np.ones(1))


def _bundle(positions, slopes, weights):
    traces = tuple(_trace(p, d) for p, d in zip(positions, slopes))
    return RayBundle(traces, np.zeros((len(traces), 2)), np.asarray(weights, float),
                     Field(), 550., _trace())


class LSAReferenceRegressionTests(unittest.TestCase):
    def test_triplet_reference_shift_preserves_all_crossings_and_paraxial_limit(self):
        from examples.sample_systems import synthetic_triplet
        system = synthetic_triplet()
        focus = paraxial_properties(system)['paraxial_image_z']
        q = [0., 1e-3, .1, 1.]
        nominal = longitudinal_spherical_aberration(system, pupil_coordinates=q)
        shifted = longitudinal_spherical_aberration(
            system, pupil_coordinates=q, reference_focus_z=focus+10.)
        assert_allclose(shifted['axis_crossing_z_mm'], nominal['axis_crossing_z_mm'],
                        rtol=0., atol=0.)
        assert_allclose(shifted['lsa_mm'], nominal['lsa_mm']-10., rtol=0., atol=2e-14)
        self.assertEqual(shifted['axis_crossing_z_mm'][0], focus)
        self.assertEqual(shifted['paraxial_focus_z_mm'], focus)
        self.assertEqual(shifted['lsa_mm'][0], -10.)
        self.assertLess(abs(shifted['lsa_mm'][1]+10.), 1e-6)

    def test_finite_object_plane_interface_uses_actual_virtual_conjugate(self):
        air = Material('air', n_constant=1.)
        glass = Material('glass', n_constant=1.5)
        system = OpticalSystem((Surface(0., 0., 2., medium_after='glass'),),
                               {'air': air, 'glass': glass}, object_z=-10., image_z=5.)
        q = np.array([0., .01, .5, 1.])
        result = longitudinal_spherical_aberration(
            system, pupil_coordinates=q, reference_focus_z=-8.)
        # Snell at a plane: dy_out=(h/sqrt(h^2+a^2))/n, z_cross=-h*dz/dy.
        h = 2*q[1:]
        dy = h/np.sqrt(h*h+100.)/1.5
        crossings = np.r_[-15., -h*np.sqrt(1-dy*dy)/dy]
        assert_allclose(result['axis_crossing_z_mm'], crossings, rtol=0., atol=2e-13)
        assert_allclose(result['lsa_mm'], crossings+8., rtol=0., atol=2e-13)
        self.assertEqual(result['statuses'][0], 'paraxial_limit')

    def test_custom_reference_does_not_invent_finite_afocal_paraxial_limit(self):
        system = OpticalSystem((Surface(0., 0., 2.),),
                               {'air': Material('air', n_constant=1.)}, image_z=10.)
        with self.assertRaises(ModelApplicabilityError):
            longitudinal_spherical_aberration(system, reference_focus_z=10.)


class GeometricSampleMembershipTests(unittest.TestCase):
    def test_zero_weight_extremum_in_spot_keeps_quadrature_moments(self):
        positions = [[0., 0., 0.], [2., 0., 0.], [10., 0., 0.]]
        b = _bundle(positions, [[0., 0.]]*3, [1., 3., 0.])
        result = spot_diagram(SimpleNamespace(image_z=0.), b)
        assert_allclose(result['coordinates_mm'], np.asarray(positions)[:, :2])
        assert_allclose(result['centroid'], [1.5, 0.], rtol=0., atol=0.)
        assert_allclose(result['covariance'], [[.75, 0.], [0., 0.]], rtol=0., atol=0.)
        self.assertEqual(result['rms_radius'], np.sqrt(.75))
        self.assertEqual(result['geo_radius_sampled'], 8.5)
        self.assertEqual(result['transmitted_quadrature_fraction'], 1.)
        assert_allclose(result['weights'], [1., 3., 0.])
        self.assertEqual(result['quadrature_mask'].tolist(), [True, True, False])
        chief = spot_diagram(SimpleNamespace(image_z=0.), b, reference='chief')
        self.assertEqual(chief['rms_radius'], np.sqrt(3.))
        self.assertEqual(chief['geo_radius_sampled'], 10.)

    def test_angular_spot_includes_zero_weight_extremum_in_same_recorded_chart(self):
        theta = np.array([.02, .06, .30])
        b = _bundle([[0., 0., 0.]]*3, np.column_stack([np.tan(theta), np.zeros(3)]),
                    [1., 3., 0.])
        result = angular_spot(b)
        assert_allclose(result['coordinates_rad'][:, 0], theta, rtol=0., atol=1e-16)
        assert_allclose(result['centroid'], [.05, 0.], rtol=0., atol=1e-16)
        self.assertAlmostEqual(result['rms_radius'], np.sqrt(.0003), places=15)
        self.assertAlmostEqual(result['geo_radius_sampled'], .25, places=15)
        assert_allclose(result['weights'], [1., 3., 0.])

    def test_common_refocus_ignores_zero_weight_for_merit_but_keeps_plot(self):
        b = _bundle([[-1., 0., 0.], [1., 0., 0.], [100., 0., 0.]],
                    [[.125, 0.], [-.125, 0.], [0., 0.]], [1., 1., 0.])
        result = common_geometric_refocus(SimpleNamespace(image_z=4.), [b])
        self.assertAlmostEqual(result['common_image_z_mm'], 8., places=13)
        self.assertAlmostEqual(result['merit_rms_mm'], 0., places=13)
        for spot in result['fixed_plane_spots']+result['common_focus_spots']:
            self.assertEqual(len(spot['coordinates_mm']), 3)
            self.assertEqual(spot['geo_radius_sampled'], 100.)
        self.assertAlmostEqual(result['common_focus_spots'][0]['rms_radius'], 0., places=13)

    def test_physical_zero_weight_losses_excluded_but_unresolved_not_hidden(self):
        b = _bundle([[0., 0., 0.], [1., 0., 0.]], [[0., 0.]]*2, [1., 0.])
        b.traces[1].status = 'aperture_blocked'
        result = spot_diagram(SimpleNamespace(image_z=0.), b)
        self.assertEqual(result['ray_indices'].tolist(), [0])
        self.assertEqual(result['geo_radius_sampled'], 0.)
        b.traces[1].status = 'intersection_unresolved'
        with self.assertRaises(NumericalError):
            spot_diagram(SimpleNamespace(image_z=0.), b)
        with self.assertRaises(NumericalError):
            angular_spot(b)
        with self.assertRaises(NumericalError):
            common_geometric_refocus(SimpleNamespace(image_z=0.), [b])

    def test_zero_weight_success_does_not_replace_missing_integration_support(self):
        b = _bundle([[0., 0., 0.], [1., 0., 0.]], [[0., 0.]]*2, [1., 0.])
        b.traces[0].status = 'aperture_blocked'
        for metric in (lambda: spot_diagram(SimpleNamespace(image_z=0.), b),
                       lambda: angular_spot(b),
                       lambda: common_geometric_refocus(SimpleNamespace(image_z=0.), [b])):
            with self.assertRaises(ModelApplicabilityError):
                metric()

    def test_real_triplet_rim_updates_maximum_without_changing_weighted_moments(self):
        from examples.sample_systems import synthetic_triplet
        system = synthetic_triplet()
        for rings in (1, 8):
            with self.subTest(rings=rings):
                b = trace_bundle(system, Field(), rings=rings, arms=24)
                without_rim = trace_bundle(system, Field(), rings=rings, arms=24,
                                           include_boundary=False)
                result = spot_diagram(system, b)
                integration_only = spot_diagram(system, without_rim)
                ids = np.flatnonzero(b.success_mask)
                expected = np.array([ray_plane_intercept(b.traces[i], system.image_z)[0]
                                     for i in ids])
                assert_allclose(result['coordinates_mm'], expected, rtol=0., atol=0.)
                assert_allclose(result['ray_indices'], ids)
                expected_max = np.linalg.norm(expected-result['centroid'], axis=1).max()
                self.assertEqual(result['geo_radius_sampled'], expected_max)
                self.assertGreater(result['geo_radius_sampled'], integration_only['geo_radius_sampled'])
                for name in ('centroid', 'covariance', 'rms_radius'):
                    assert_allclose(result[name], integration_only[name], rtol=0., atol=0.)
                self.assertEqual(result['geometric_radius_scope'], 'all_successful_traced_samples')


if __name__ == '__main__':
    unittest.main()
