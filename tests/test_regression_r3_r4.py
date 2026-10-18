"""Physical overlap and auxiliary-cylinder controls for review defects R3/R4."""
import math
import unittest
from dataclasses import replace

import numpy as np

from sequential_optics.errors import GeometryError, NumericalError
from sequential_optics.geometry import (
    _bound_interval, aim_ray, intersect, launch_plane_z, sag, trace_bundle, trace_ray,
)
from sequential_optics.materials import Material
from sequential_optics.metrics import spot_diagram
from sequential_optics.models import Field, OpticalSystem, Ray, Surface
from sequential_optics.settings import EngineSettings


def two_apertures(front_mechanical=None):
    return OpticalSystem(
        (Surface(0., 0., 2., mechanical_radius=front_mechanical), Surface(10., 0., 2.)),
        {'air': Material('air', n_constant=1.)}, stop_index=1, image_z=20.,
    )


class PreStopVignettingRegression(unittest.TestCase):
    def test_offset_plane_apertures_match_independent_overlap_quadrature(self):
        # STOP x -> first-plane x = STOP x - spacing*tan(field).
        # This supplies an independent per-ray oracle, not a stored engine result.
        field = Field(.1, 0.)
        baseline = two_apertures()
        wide_rim = two_apertures(4.)
        outputs = []
        for system in (baseline, wide_rim):
            bundle = trace_bundle(system, field, rings=8, arms=24, include_boundary=False)
            shift = np.array([10.*math.tan(field.x), 0.])
            expected = np.linalg.norm(2.*bundle.pupil_xy-shift, axis=1) <= 2.
            np.testing.assert_array_equal(bundle.success_mask, expected)
            self.assertEqual(int(expected.sum()), 134)
            self.assertTrue(all(t.status in ('ok', 'escaped', 'aperture_blocked')
                                for t in bundle.traces))
            weights = bundle.weights[expected]
            coords = 2.*bundle.pupil_xy[expected]+shift
            centroid = (weights/weights.sum())@coords
            rms = math.sqrt(float((weights/weights.sum())@np.sum((coords-centroid)**2, axis=1)))
            output = spot_diagram(system, bundle)
            self.assertAlmostEqual(output['transmitted_quadrature_fraction'], weights.sum(), places=14)
            self.assertAlmostEqual(output['rms_radius'], rms, places=13)
            outputs.append(output)
        # Changing only an opaque mechanical rim cannot change transmitted rays.
        self.assertAlmostEqual(outputs[0]['rms_radius'], outputs[1]['rms_radius'], places=14)

    def test_default_zero_weight_rim_no_longer_aborts_the_metric(self):
        system = two_apertures(4.)
        bundle = trace_bundle(system, Field(.1, 0.), rings=8, arms=24)
        shift = np.array([10.*math.tan(.1), 0.])
        expected = np.linalg.norm(2.*bundle.pupil_xy-shift, axis=1) <= 2.+1e-10
        np.testing.assert_array_equal(bundle.success_mask, expected)
        self.assertFalse(any(t.status == 'aim_unresolved' for t in bundle.traces))
        self.assertEqual(spot_diagram(system, bundle)['status'], 'ok')

    def test_curved_dummy_extension_does_not_enlarge_the_physical_aperture(self):
        system = replace(two_apertures(), surfaces=(Surface(0., .02, 1.), Surface(10., 0., 2.)))
        field = Field(.1, 0.)
        launch = aim_ray(system, field, [-.9, 0.])
        # No refracting-index change: the extended sag must leave the straight
        # incident direction intact while the actual finite surface rejects it.
        expected_x = -1.8-(10.-launch.origin[2])*math.tan(.1)
        self.assertAlmostEqual(launch.origin[0], expected_x, places=12)
        self.assertEqual(trace_ray(system, launch).status, 'escaped')

    def test_refracting_sphere_clipping_matches_independent_snell_threshold(self):
        # At normal object incidence, a spherical interface maps each radial
        # height monotonically to the STOP. Independently refract the clear rim
        # to obtain the limiting STOP radius; this predicts every sample's mask.
        materials = {'air': Material('air', n_constant=1.),
                     'glass': Material('glass', n_constant=1.5)}
        system = OpticalSystem((Surface(0., .02, 1., medium_after='glass'),
                                Surface(10., 0., 2., medium_after='glass')),
                               materials, stop_index=1, image_z=20.)
        r, radius, eta = 1., 50., 1./1.5
        cosine = math.sqrt(1.-(r/radius)**2)
        coefficient = eta*cosine-math.sqrt(1.-eta**2*(r/radius)**2)
        radial_direction = coefficient*r/radius
        axial_direction = eta-coefficient*cosine
        sag = radius*(1.-cosine)
        limiting_stop_radius = r+(10.-sag)*radial_direction/axial_direction
        narrow = trace_bundle(system, Field(), rings=8, arms=24, include_boundary=False)
        expected = 2.*np.linalg.norm(narrow.pupil_xy, axis=1) <= limiting_stop_radius
        self.assertTrue(np.any(expected) and np.any(~expected))
        np.testing.assert_array_equal(narrow.success_mask, expected)
        self.assertTrue(all(t.status in ('ok', 'escaped', 'aperture_blocked') for t in narrow.traces))
        widened = replace(system, surfaces=(replace(system.surfaces[0], mechanical_radius=4.),
                                             system.surfaces[1]))
        control = trace_bundle(widened, Field(), rings=8, arms=24, include_boundary=False)
        np.testing.assert_array_equal(control.success_mask, expected)
        for a, b, selected in zip(narrow.traces, control.traces, expected):
            if selected:
                np.testing.assert_allclose(a.positions, b.positions, atol=2e-10)
                np.testing.assert_allclose(a.directions, b.directions, atol=2e-12)

    def test_virtual_surface_beyond_target_is_not_a_valid_prefix(self):
        system = replace(two_apertures(), surfaces=(
            Surface(0., 0., .5, coefficients=(.5,), surface_type='even_asphere'),
            Surface(1., 0., 2.)))
        with self.assertRaisesRegex(NumericalError, 'aim_unresolved.*beyond STOP target'):
            aim_ray(system, Field(), [.9, 0.])

    def test_target_plane_residual_cannot_accept_a_later_stop_root(self):
        # z=x² and z=3x-2 cross at (x,z)=(1,1) and (2,4). The
        # desired STOP target is the second crossing, so the direction/target
        # combination is unavailable on this first-hit branch.
        system = OpticalSystem((Surface(0., 0., 4., coefficients=(1.,),
                                        surface_type='even_asphere'),),
                               {'air': Material('air', n_constant=1.)})
        with self.assertRaisesRegex(NumericalError, 'aim_unresolved'):
            aim_ray(system, Field(math.atan(1./3.), 0.), [.5, 0.])

    def test_undefined_virtual_conic_is_unresolved_not_zero_power(self):
        system = replace(two_apertures(), surfaces=(Surface(0., 1., .5), Surface(10., 0., 2.)))
        # Clear surface is valid, but x=1.8 is beyond this sphere's real sag.
        with self.assertRaisesRegex(NumericalError, 'aim_unresolved.*virtual sag continuation'):
            aim_ray(system, Field(), [.9, 0.])

    def test_invalid_real_clear_surface_domain_is_not_vignetting(self):
        system = replace(two_apertures(), surfaces=(Surface(0., 1., 2.), Surface(10., 0., 2.)))
        with self.assertRaisesRegex(GeometryError, 'sag_domain'):
            trace_bundle(system, Field(), rings=2, arms=4)


class BoundingCylinderRoundoffRegression(unittest.TestCase):
    def test_negative_roundoff_discriminant_keeps_optical_plane_hit(self):
        system = two_apertures(4.)
        pupil = np.array([math.cos(3.*math.pi/2.), math.sin(3.*math.pi/2.)])
        target = 2.*pupil
        z = launch_plane_z(system)
        d = np.array([math.tan(.1), 0., 1.]); d /= np.linalg.norm(d)
        o = np.array([*(target-(10.-z)*d[:2]/d[2]), z])
        prefix = trace_ray(system, Ray(o, d), stop_at=0, ignore_apertures=True)
        incoming = Ray(prefix.last_position, prefix.last_direction)
        a = float(incoming.direction[:2]@incoming.direction[:2])
        b = 2.*float(incoming.origin[:2]@incoming.direction[:2])
        c = float(incoming.origin[:2]@incoming.origin[:2])-4.
        self.assertLess(b*b-4*a*c, 0.)  # Original IEEE double rounding fixture.
        lo, hi = _bound_interval(system.surfaces[1], incoming, EngineSettings())
        expected_t = (10.-incoming.origin[2])/incoming.direction[2]
        self.assertLessEqual(lo, expected_t+1e-10)
        self.assertGreaterEqual(hi, expected_t-1e-10)
        hit, distance = intersect(system.surfaces[1], incoming)
        np.testing.assert_allclose(hit, [*target, 10.], atol=1e-13)
        self.assertAlmostEqual(distance, expected_t, places=13)
        self.assertGreater(incoming.direction[2], .99)  # Not optical tangency.

    def test_true_cylinder_miss_remains_miss_over_length_scales(self):
        for scale in (.001, 1., 1000.):
            with self.subTest(scale=scale):
                surface = Surface(10.*scale, 0., 2.*scale)
                ray = Ray([-1.*scale, -(2.+1e-5)*scale, 0.], [.1, 0., 1.])
                with self.assertRaisesRegex(GeometryError, 'escaped'):
                    _bound_interval(surface, ray, EngineSettings())
                with self.assertRaisesRegex(GeometryError, 'escaped'):
                    intersect(surface, ray)

    def test_auxiliary_cylinder_tangency_can_cross_curved_sag_normally(self):
        sphere = Surface(0., .1, 2.)
        surface = replace(sphere, z=10.-float(sag(sphere, 0., -2.)))
        ray = Ray([-10.*math.tan(.1), -2., 0.], [math.tan(.1), 0., 1.])
        hit, _ = intersect(surface, ray)
        np.testing.assert_allclose(hit, [0., -2., 10.], atol=2e-13)

    def test_refracting_curved_stop_rim_uses_same_physical_radial_tolerance(self):
        # Independent plane Snell fixes the prefix direction; the aimed sphere
        # rim is slightly outside r=2 in floating arithmetic, but well within
        # the declared radial tolerance. Its root must remain bracketed.
        system = OpticalSystem((Surface(0., 0., 2., medium_after='glass'),
                                Surface(10., .02, 2.)),
                               {'air': Material('air', n_constant=1.),
                                'glass': Material('glass', n_constant=1.5)},
                               stop_index=1)
        field = Field(.2, 0.)
        launch = aim_ray(system, field, [1., 0.])
        trace = trace_ray(system, launch)
        self.assertTrue(trace.success, trace.message)
        stop_z = 10.+50.*(1.-math.sqrt(1.-(2./50.)**2))
        sine = math.sin(.2)/1.5
        expected_front_x = 2.-stop_z*sine/math.sqrt(1.-sine*sine)
        self.assertAlmostEqual(trace.positions[0, 0], expected_front_x, places=9)
        np.testing.assert_allclose(trace.positions[1], [2., 0., stop_z], atol=1e-9)

    def test_actual_optical_tangency_is_still_rejected(self):
        surface = Surface(0., 0., 3., coefficients=(1.,), surface_type='even_asphere')
        with self.assertRaisesRegex(GeometryError, 'grazing_degenerate'):
            intersect(surface, Ray([-2., 0., 0.], [1., 0., 0.]))


if __name__ == '__main__':
    unittest.main()
