"""Coincident sequential interfaces must preserve certified real-ray paths."""
from dataclasses import replace
import unittest

import numpy as np
from numpy.testing import assert_allclose

from sequential_optics import Field, OpticalSystem, Surface
from sequential_optics.errors import GeometryError
from sequential_optics.geometry import (aim_ray, gradient, intersect, sag,
                                       trace_bundle, trace_ray, validate_geometry)
from sequential_optics.materials import Material
from sequential_optics.metrics import spot_diagram
from sequential_optics.models import Ray


def singlet():
    materials = {'air': Material('synthetic air', n_constant=1.),
                 'g': Material('synthetic glass', n_constant=1.5),
                 'intermediate': Material('synthetic intermediate', n_constant=1.2)}
    return OpticalSystem((Surface(0., .02, 10., medium_after='g'),
                          Surface(3., -.02, 10.)),
                         materials, stop_radius=5., image_z=52.4949494949)


class CoincidentInterfaceRegression(unittest.TestCase):
    def assert_equivalent_bundle(self, original, duplicate, field):
        a, b = (trace_bundle(s, field) for s in (original, duplicate))
        self.assertTrue(a.success_mask.all())
        assert_allclose(b.success_mask, a.success_mask, atol=0.)
        assert_allclose(b.pupil_xy, a.pupil_xy, atol=0.)
        assert_allclose(b.weights, a.weights, atol=0.)
        for left, right in zip(a.traces, b.traces):
            # The duplicate contributes exactly zero geometric/optical length.
            assert_allclose(right.positions[1], right.positions[0], atol=0., rtol=0.)
            self.assertEqual(right.opl[1], right.opl[0])
            assert_allclose(right.directions[1], right.directions[0], atol=2e-15)
            assert_allclose(right.positions[[0, 2]], left.positions, atol=2e-9, rtol=0.)
            assert_allclose(right.directions[[0, 2]], left.directions, atol=5e-11, rtol=0.)
            assert_allclose(right.opl[[0, 2]], left.opl, atol=2e-9, rtol=0.)
            assert_allclose(gradient(duplicate.surfaces[1], *right.positions[1, :2]),
                            gradient(original.surfaces[0], *left.positions[0, :2]),
                            atol=1e-10, rtol=0.)
        spots = [spot_diagram(s, q) for s, q in ((original, a), (duplicate, b))]
        for key in ('centroid', 'rms_radius', 'geo_radius_sampled', 'transmitted_quadrature_fraction'):
            assert_allclose(spots[0][key], spots[1][key], atol=3e-9, rtol=0.)
        return a, b

    def test_original_217_ray_regression_on_axis_and_two_off_axis_fields(self):
        base = singlet()
        duplicate = replace(base, surfaces=(base.surfaces[0], *base.surfaces))
        self.assertEqual(validate_geometry(duplicate)['status'], 'valid')
        for field in (Field(), Field(.01, .02), Field(0., .1)):
            with self.subTest(field=field):
                a, b = self.assert_equivalent_bundle(base, duplicate, field)
                self.assertEqual(len(a.traces), 217)
                self.assertEqual(sum(b.success_mask), 217)

    def test_nonzero_accepted_sag_residual_is_preserved_not_snapped(self):
        base = singlet()
        duplicate = replace(base, surfaces=(base.surfaces[0], *base.surfaces))
        ray = aim_ray(base, Field(.01, .02), [.8, 0.])
        old, new = (trace_ray(s, ray) for s in (base, duplicate))
        self.assertTrue(old.success)
        self.assertTrue(new.success, new.message)
        residual = old.positions[0, 2]-float(sag(base.surfaces[0], *old.positions[0, :2]))
        self.assertGreater(residual, 0.)  # Original confirmed failing case.
        self.assertLess(residual, 1e-10)
        assert_allclose(new.positions[1], old.positions[0], atol=0., rtol=0.)
        self.assertEqual(new.opl[1], new.opl[0])

    def test_dummy_before_internal_stop_uses_virtual_aiming(self):
        base = replace(singlet(), stop_index=1, stop_radius=3.)
        duplicate = replace(base, surfaces=(base.surfaces[0], *base.surfaces), stop_index=2)
        self.assert_equivalent_bundle(base, duplicate, Field(.04, -.03))

    def test_coincident_dummy_is_stop_and_uses_first_hit_verification(self):
        base = singlet()
        duplicate = replace(base, surfaces=(base.surfaces[0], *base.surfaces), stop_index=1)
        self.assert_equivalent_bundle(base, duplicate, Field(.04, -.03))

    def test_rear_stop_coincident_chain_aims_as_one_position_event(self):
        # Independent follow-up review: trial rays at the first rear member
        # need not already land at target z. It is the STOP event itself, not a
        # distinct prefix to reject before Newton has had a chance to converge.
        materials = {'air': Material('air', n_constant=1.),
                     'glass': Material('glass', n_constant=1.52)}
        front = Surface(0., .025, 8., conic=-.4, coefficients=(.0004, -1e-6),
                        surface_type='even_asphere', medium_after='glass')
        rear = Surface(4., -.02, 8.)
        base = OpticalSystem((front, rear), materials, stop_index=1,
                             stop_radius=3., image_z=45.)
        chain = (front, replace(front, coefficients=front.coefficients+(0.,)*4),
                 replace(front, name='observation'), rear,
                 replace(rear, surface_type='even_asphere', coefficients=(0.,)*6))
        duplicate = replace(base, surfaces=chain, stop_index=4)
        for field in (Field(.045, -.07), Field(-.08, .025)):
            for pupil in ((.2, -.4), (-.5, .3), (.7, .5), (-.8, -.2)):
                with self.subTest(field=field, pupil=pupil):
                    a, b = [trace_ray(s, aim_ray(s, field, pupil)) for s in (base, duplicate)]
                    self.assertTrue(a.success and b.success, (a.message, b.message))
                    assert_allclose(a.positions, b.positions[[0, 3]], atol=2e-9, rtol=0.)
                    assert_allclose(a.directions, b.directions[[0, 3]], atol=5e-11, rtol=0.)
                    assert_allclose(a.opl, b.opl[[0, 3]], atol=2e-9, rtol=0.)
                    assert_allclose(b.positions[3, :2], 3*np.asarray(pupil), atol=1e-9, rtol=0.)
                    assert_allclose(b.positions[3], b.positions[4], atol=0., rtol=0.)
                    self.assertEqual(b.opl[3], b.opl[4])

    def test_even_asphere_and_finite_object_coincident_prefix(self):
        original = singlet()
        front = replace(original.surfaces[0], surface_type='even_asphere', conic=-.7,
                        coefficients=(.001, -2e-6, 1e-9, -1e-12, 1e-16, -1e-20))
        base = replace(original, surfaces=(front, original.surfaces[1]),
                       stop_index=1, stop_radius=3., object_z=-100.)
        duplicate = replace(base, surfaces=(front, *base.surfaces), stop_index=2)
        self.assert_equivalent_bundle(base, duplicate, Field(.4, -.2, kind='height'))

    def test_dummy_active_aperture_is_still_applied(self):
        base = singlet()
        dummy = replace(base.surfaces[0], clear_radius=2., mechanical_radius=10.)
        system = replace(base, surfaces=(base.surfaces[0], dummy, base.surfaces[1]))
        reference = trace_bundle(base, Field(.01, .02))
        actual = trace_bundle(system, Field(.01, .02))
        expected = np.array([np.hypot(*t.positions[0, :2]) <= 2.+1e-10 for t in reference.traces])
        self.assertTrue(expected.any() and (~expected).any())
        np.testing.assert_array_equal(actual.success_mask, expected)
        self.assertTrue(all(t.status == 'aperture_blocked' for t, ok in zip(actual.traces, expected) if not ok))

    def test_dummy_mechanical_domain_is_still_applied(self):
        base = singlet()
        dummy = replace(base.surfaces[0], clear_radius=2., mechanical_radius=2.)
        system = replace(base, surfaces=(base.surfaces[0], dummy, base.surfaces[1]))
        trace = trace_ray(system, aim_ray(base, Field(.01, .02), [.8, 0.]))
        self.assertEqual(trace.status, 'escaped')
        self.assertEqual(len(trace.positions), 1)

    def test_coincident_refracting_interfaces_follow_net_snell(self):
        base = singlet()
        first = replace(base.surfaces[0], medium_after='intermediate')
        composed = replace(base, surfaces=(first, *base.surfaces))
        ray = aim_ray(base, Field(.1, -.02), [.8, .2])
        a, b = (trace_ray(s, ray) for s in (base, composed))
        self.assertTrue(b.success, b.message)
        assert_allclose(b.positions[[1, 2]], a.positions, atol=2e-12)
        assert_allclose(b.directions[[1, 2]], a.directions, atol=2e-15)
        assert_allclose(b.opl[[1, 2]], a.opl, atol=2e-12)
        self.assertEqual(b.opl[0], b.opl[1])

    def test_nearby_but_distinct_surface_keeps_nonzero_distance(self):
        base = singlet()
        nearby = replace(base.surfaces[0], z=1e-6)
        system = replace(base, surfaces=(base.surfaces[0], nearby, base.surfaces[1]))
        trace = trace_ray(system, aim_ray(base, Field(.01, .02), [.8, 0.]))
        self.assertTrue(trace.success, trace.message)
        distance = np.linalg.norm(trace.positions[1]-trace.positions[0])
        self.assertGreater(distance, 5e-7)
        self.assertAlmostEqual(trace.opl[1]-trace.opl[0], 1.5*distance, places=12)

    def test_same_vertex_distinct_asphere_keeps_nonzero_distance(self):
        base = singlet()
        distinct = replace(base.surfaces[0], surface_type='even_asphere', coefficients=(1e-4,))
        system = replace(base, surfaces=(base.surfaces[0], distinct, base.surfaces[1]))
        trace = trace_ray(system, aim_ray(base, Field(.01, .02), [.8, 0.]))
        self.assertTrue(trace.success, trace.message)
        self.assertGreater(trace.opl[1]-trace.opl[0], 1e-3)

    def test_independent_ray_small_residual_is_not_a_hit_certificate(self):
        surface = singlet().surfaces[0]
        p = np.array([4., 0., float(sag(surface, 4., 0.))+1e-13])
        with self.assertRaisesRegex(GeometryError, 'escaped'):
            intersect(surface, Ray(p, [0., 0., 1.]))

    def test_tangent_and_first_simple_root_controls_remain(self):
        surface = Surface(0., 0., 4., coefficients=(1.,), surface_type='even_asphere')
        with self.assertRaisesRegex(GeometryError, 'grazing_degenerate'):
            intersect(surface, Ray([-2., 0., 0.], [1., 0., 0.]))
        point, distance = intersect(surface, Ray([0., 0., -2.], [1., 0., 3.]))
        assert_allclose(point, [1., 0., 1.], atol=1e-9)
        self.assertAlmostEqual(distance, np.sqrt(10), places=8)


if __name__ == '__main__':
    unittest.main()
