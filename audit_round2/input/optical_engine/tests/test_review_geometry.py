"""Independent geometry/thermal review regressions beyond module-owner tests."""
import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose

from sequential_optics.models import Surface, OpticalSystem, Field, Ray
from sequential_optics.materials import Material
from sequential_optics.settings import EngineSettings
from sequential_optics.geometry import aim_ray, trace_ray, intersect, refract, gradient, hessian
from sequential_optics.paraxial import paraxial_properties, resolve_aperture
from sequential_optics.thermal import thermal_system, LensGroup, scale_surface
from sequential_optics.errors import MissingDataError


class IndependentGeometryPhysicsReview(unittest.TestCase):
    def test_oblique_plate_output_eikonal_with_absolute_nonunit_ambient(self):
        # Analytic S_exit(x,y)=n0*d_xy·(x,y)+sqrt(ng²-n0²|d_xy|²)*thickness.
        # This checks launch phase + refracted optical paths together, not just
        # the raw sum of geometric path lengths along one ray.
        n0, ng, thickness = 1.2, 1.6, 4.
        mats = {'air': Material('synthetic ambient', n_constant=n0),
                'glass': Material('synthetic glass', n_constant=ng)}
        system = OpticalSystem((Surface(0., 0., 10., medium_after='glass'), Surface(thickness, 0., 10.)),
                               mats, stop_radius=6.)
        field = Field(.5, -.3)
        for pupil in ([0., 0.], [.4, .2], [-.5, -.3]):
            ray = aim_ray(system, field, pupil)
            result = trace_ray(system, ray)
            self.assertTrue(result.success, result.message)
            kz = np.sqrt(ng*ng-n0*n0*np.dot(ray.direction[:2], ray.direction[:2]))
            expected = n0*np.dot(ray.direction[:2], result.last_position[:2]) + kz*thickness
            self.assertAlmostEqual(result.last_opl, expected, places=10)
            assert_allclose(result.last_direction, ray.direction, atol=2e-15)

    def test_exact_reversed_asphere_path_recovers_original_ray(self):
        from examples.sample_systems import synthetic_triplet
        system = synthetic_triplet(aspheric=True)
        incident = Ray([.4, -.2, -5.], [.01, .015, 1.])
        forward = trace_ray(system, incident)
        self.assertTrue(forward.success, forward.message)
        p = forward.last_position + .25*forward.last_direction
        d = -forward.last_direction
        medium_after = [s.medium_after for s in system.surfaces]
        for j in range(len(system.surfaces)-1, -1, -1):
            surface = system.surfaces[j]
            hit, _ = intersect(surface, Ray(p, d))
            assert_allclose(hit, forward.positions[j], atol=2e-8)
            n1 = system.index(medium_after[j], 550.)
            n2 = system.index(system.ambient_material if j == 0 else medium_after[j-1], 550.)
            g = gradient(surface, *hit[:2]); normal = [-g[0], -g[1], 1.]
            d = refract(d, normal, n1, n2) if n1 != n2 else d
            p = hit
        assert_allclose(d, -incident.direction, atol=2e-10)
        distance = (incident.origin[2]-p[2])/d[2]
        assert_allclose(p+distance*d, incident.origin, atol=2e-8)

    def test_reduced_index_EFL_and_nominal_aperture_are_distinct(self):
        # Single refracting sphere: power=(1.6-1.2)/50=.008/mm,
        # image-space EFL=1.6/.008=200mm, ambient-normalized EFL=150mm.
        mats = {'air': Material('synthetic ambient', n_constant=1.2),
                'glass': Material('synthetic glass', n_constant=1.6)}
        system = OpticalSystem((Surface(0., .02, 15., medium_after='glass'),), mats)
        result = paraxial_properties(system)
        self.assertAlmostEqual(result['efl_image_mm'], 200.)
        self.assertAlmostEqual(result['efl_report_mm'], 150.)
        self.assertAlmostEqual(result['determinant'], 1.)
        ray = trace_ray(system, Ray([1e-4, 0., -1.], [0., 0., 1.]))
        real_limit = ray.last_position[2]-ray.last_position[0]*ray.last_direction[2]/ray.last_direction[0]
        self.assertAlmostEqual(real_limit, 200., places=7)
        resolved = resolve_aperture(system, f_number=10.)
        self.assertAlmostEqual(resolved.stop_semi_diameter, 7.5)


class IndependentThermalPhysicsReview(unittest.TestCase):
    def test_expansion_obeys_generic_model_validity_range(self):
        material = Material('bounded synthetic material', n_constant=1.5,
                            cte_per_k=1e-5, temperature_range_c=(0., 85.))
        with self.assertRaises(MissingDataError):
            material.expansion_scale(200., 20.)
        with self.assertRaises(MissingDataError):
            material.expansion_scale(20., 200.)

    def test_scaled_asphere_gradient_and_curvature_transform(self):
        original = Surface(0., .02, 3., conic=-.5, coefficients=(.001, 1e-5, -1e-8, 1e-11, -1e-14, 1e-17),
                           surface_type='even_asphere')
        s = 1.02
        hot = scale_surface(original, s)
        for xy in ([0., 0.], [.2, -.3], [1.7, .8]):
            xy = np.asarray(xy)
            assert_allclose(gradient(hot, *(s*xy)), gradient(original, *xy), rtol=5e-15, atol=1e-16)
            assert_allclose(hessian(hot, *(s*xy)), hessian(original, *xy)/s, rtol=5e-15, atol=1e-16)

    def test_thermal_whole_system_similarity_under_explicit_anchor_motion(self):
        # Analytic similarity is a stronger oracle than expecting only a plausible
        # trend in focal length. All paths scale, all Snell directions stay fixed.
        # Anchor motions here are explicitly prescribed mathematical test data.
        from examples.sample_systems import synthetic_triplet
        base = synthetic_triplet(aspheric=True)
        alpha, dt = 1e-4, 100.
        mats = {k: replace(v, index_temperature_invariant=True, cte_per_k=alpha)
                for k, v in base.materials.items()}
        base = replace(base, materials=mats)
        groups = tuple(LensGroup(i, i+1, anchor_shift_per_k_mm=alpha*(base.surfaces[i].z+base.surfaces[i+1].z)/2)
                       for i in (0, 2, 4))
        hot = thermal_system(base, 120., lens_groups=groups)
        s = 1+alpha*dt
        assert_allclose([v.z for v in hot.surfaces], s*np.array([v.z for v in base.surfaces]), atol=2e-15)
        original = Ray([.3, -.2, -4.], [.01, .02, 1.])
        a = trace_ray(base, original)
        b = trace_ray(hot, Ray(s*original.origin, original.direction))
        self.assertTrue(a.success, a.message); self.assertTrue(b.success, b.message)
        assert_allclose(b.positions, s*a.positions, atol=1e-9)
        assert_allclose(b.directions, a.directions, atol=2e-12)
        assert_allclose(b.opl, s*a.opl, atol=1e-9)
        self.assertAlmostEqual(paraxial_properties(hot)['efl_image_mm']/paraxial_properties(base)['efl_image_mm'], s)


if __name__ == '__main__': unittest.main()
