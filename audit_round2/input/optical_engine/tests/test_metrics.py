"""Independent analytic or symmetry references; no Zemax equivalence claim.

Sample choices were challenged by the sample referee: use unequal quadratures,
chief/centroid separation, clipped-pupil rank tests, rotated fields, and retain
nominal scale drift even when refitted distortion is zero.
"""
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from numpy.testing import assert_allclose

from sequential_optics.models import Field, RayBundle, TraceResult
from sequential_optics.errors import InputError, NumericalError, ModelApplicabilityError
from sequential_optics.metrics import (
    spot_diagram, angular_spot, ray_plane_intercept, tangential_sagittal_basis,
    focus_from_ray_lines, common_geometric_refocus, zernike_terms, zernike_basis,
    zernike_fit, map_distortion, coma_from_fan, differential_field_curvature,
    chief_mapping, distortion, ray_fans,
    longitudinal_spherical_aberration, grid_mapping,
)


def trace(position, slope=(0., 0.)):
    d = np.array([*slope, 1.]); d /= np.linalg.norm(d)
    return TraceResult(positions=np.array([position], float), directions=np.array([d]),
                       opl=np.array([0.]), powers=np.array([1.]))


def bundle(positions, slopes=None, weights=None, field=None, chief=None):
    positions = np.array(positions, float)
    if slopes is None:
        slopes = np.zeros((len(positions), 2))
    traces = tuple(trace(p, s) for p, s in zip(positions, slopes))
    return RayBundle(traces, np.zeros((len(positions), 2)),
                     np.ones(len(traces)) if weights is None else np.array(weights),
                     field or Field(), 550., chief or trace([0., 0., 0.]))


class SpotTests(unittest.TestCase):
    def test_weighted_spot_chief_centroid_distinction(self):
        b = bundle([[0, 0, 0], [2, 0, 0], [0, 3, 0]], weights=[1, 2, 1])
        s = SimpleNamespace(image_z=0.)
        centroid = spot_diagram(s, b)
        chief = spot_diagram(s, b, reference='chief')
        assert_allclose(centroid['centroid'], [1., .75])
        self.assertAlmostEqual(centroid['rms_radius']**2, 2.6875)
        self.assertAlmostEqual(chief['rms_radius']**2, 4.25)
        assert_allclose(centroid['chief_xy_mm'], [0., 0.])
        assert_allclose(centroid['basis_T'], [0, 1])
        self.assertAlmostEqual(centroid['sigma_T_mm']**2, 1.6875)

    def test_rotating_field_and_rays_preserves_TS(self):
        angle = .63
        r = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        xy = np.array([[0., 0.], [1., .2], [-.3, 2.], [.5, -.1]])
        pos = np.column_stack([xy, np.zeros(4)])
        rotated = np.column_stack([xy @ r.T, np.zeros(4)])
        f = np.array([.3, .7]); fr = r @ f
        b1 = bundle(pos, field=Field(*f, kind='height'))
        b2 = bundle(rotated, field=Field(*fr, kind='height'))
        a = spot_diagram(SimpleNamespace(image_z=0.), b1)
        b = spot_diagram(SimpleNamespace(image_z=0.), b2)
        assert_allclose(a['coordinates_TS_mm'], b['coordinates_TS_mm'], atol=1e-14)
        self.assertAlmostEqual(a['rms_radius'], b['rms_radius'])

    def test_physical_blocking_differs_from_unresolved(self):
        b = bundle([[0, 0, 0], [2, 0, 0]])
        b.traces[1].status = 'aperture_blocked'
        out = spot_diagram(SimpleNamespace(image_z=1.), b)
        self.assertEqual(out['transmitted_quadrature_fraction'], .5)
        b.traces[1].status = 'intersection_unresolved'
        with self.assertRaises(NumericalError):
            spot_diagram(SimpleNamespace(image_z=1.), b)

    def test_virtual_and_grazing_intercepts_explicit(self):
        b = bundle([[0, 0, 2]], slopes=[[.2, .1]])
        out = spot_diagram(SimpleNamespace(image_z=0.), b)
        self.assertTrue(out['virtual_intercepts'][0])
        assert_allclose(out['coordinates_mm'][0], [-.4, -.2])
        with self.assertRaises(ModelApplicabilityError):
            spot_diagram(SimpleNamespace(image_z=0.), b, allow_virtual=False)
        tr = trace([0, 0, 0]); tr.directions[0] = [1, 0, 0]
        with self.assertRaises(ModelApplicabilityError):
            ray_plane_intercept(tr, 2.)

    def test_angular_spot_uses_directions_not_output_position(self):
        b = bundle([[1, 2, 0], [-1, 3, 0]], slopes=[[.1, 0], [-.1, 0]])
        out = angular_spot(b)
        self.assertAlmostEqual(out['rms_radius'], np.arctan(.1))
        shifted = bundle([[40, 90, 200], [-6, 7, 200]], slopes=[[.1, 0], [-.1, 0]])
        assert_allclose(out['coordinates_rad'], angular_spot(shifted)['coordinates_rad'])


class FocusTests(unittest.TestCase):
    def test_one_common_focus_not_individual_envelope(self):
        a = np.array([[-1., 0.], [1., 0.]])
        result = focus_from_ray_lines([a, 2*a], [-a/10, -a/10])
        self.assertAlmostEqual(result['common_image_z_mm'], 15.)
        assert_allclose(result['local_best_focus_diagnostic_mm'], [10, 20])
        assert_allclose(result['per_field_rms_mm'], [.5, .5])
        self.assertGreater(result['merit_rms_mm'], 0.)

    def test_unequal_pupil_and_field_weights_match_direct_minimum(self):
        a = [np.array([[-1., 0.], [2., 0.], [0., 1.]]), np.array([[0., 2.], [2., 3.], [-1., 1.]])]
        b = [-a[0]/8, -a[1]/20]
        weights = [np.array([1., 4., 2.]), np.array([3., 2., 1.])]
        fieldweights = np.array([1., 3.])
        result = focus_from_ray_lines(a, b, pupil_weights=weights, field_weights=fieldweights)
        z = result['common_image_z_mm']
        def objective(z):
            total = 0.
            for ah, bh, wh, fh in zip(a, b, weights, fieldweights/fieldweights.sum()):
                x = ah + z*bh; wn = wh/wh.sum(); center = wn@x
                total += fh*np.sum(wn*np.sum((x-center)**2, axis=1))
            return total
        self.assertLess(objective(z), objective(z-.001))
        self.assertLess(objective(z), objective(z+.001))
        self.assertAlmostEqual(result['merit_rms_mm']**2, objective(z))
        # Multiplying one field's pupil weights must NOT change field weight.
        scaled = focus_from_ray_lines(a, b, pupil_weights=[100*weights[0], weights[1]], field_weights=fieldweights)
        self.assertAlmostEqual(z, scaled['common_image_z_mm'])

    def test_common_api_evaluates_every_field_same_plane(self):
        aa = np.array([[-1., 0., 0.], [1., 0., 0.]])
        b1 = bundle(aa, -aa[:, :2]/10)
        b2 = bundle(2*aa, -aa[:, :2]/10)
        result = common_geometric_refocus(SimpleNamespace(image_z=5.), [b1, b2])
        self.assertAlmostEqual(result['common_image_z_mm'], 15.)
        assert_allclose([s['image_z_mm'] for s in result['common_focus_spots']], [15., 15.])
        assert_allclose([s['rms_radius'] for s in result['common_focus_spots']], [.5, .5])

    def test_bounds_and_collimated_degenerate_case(self):
        a = np.array([[-1., 0.], [1., 0.]])
        bounded = focus_from_ray_lines([a], [-a/10], bounds=(2, 7))
        self.assertEqual(bounded['common_image_z_mm'], 7.)
        degenerate = focus_from_ray_lines([a], [np.zeros_like(a)])
        self.assertEqual(degenerate['status'], 'no_unique_finite_focus')
        self.assertIsNone(degenerate['common_image_z_mm'])


class ZernikeTests(unittest.TestCase):
    @staticmethod
    def disk():
        # Uniform-area disk quadrature in t=rho^2, independent of implementation.
        t, w = np.polynomial.legendre.leggauss(16)
        radius = np.sqrt((t+1)/2)
        phi = np.linspace(0, 2*np.pi, 64, endpoint=False)
        x = (radius[:, None]*np.cos(phi)).ravel()
        y = (radius[:, None]*np.sin(phi)).ravel()
        return np.column_stack([x, y]), np.repeat(w/2/len(phi), len(phi))

    def test_full_disk_normalization_and_known_polynomials(self):
        xy, w = self.disk(); terms = zernike_terms(6)
        basis = zernike_basis(xy, terms)
        assert_allclose(basis.T @ (w[:, None]*basis), np.eye(len(terms)), atol=2e-14)
        x, y = xy.T; r2 = x*x+y*y
        assert_allclose(basis[:, terms.index((3, 1, 'c'))], np.sqrt(8)*(3*r2-2)*x, atol=3e-15)
        assert_allclose(basis[:, terms.index((4, 0, '0'))], np.sqrt(5)*(6*r2*r2-6*r2+1), atol=3e-15)

    def test_clipped_nonuniform_pupil_exact_fit_immutable_OPD(self):
        xy, w = self.disk(); x, y = xy.T; r2 = x*x+y*y
        opd = .0007+2e-4*x+3e-4*np.sqrt(3)*(2*r2-1)+1e-4*np.sqrt(8)*(3*r2-2)*y
        original = opd.copy()
        mask = (x > -.4) & (y < .8)
        result = zernike_fit(xy, opd, n_max=3, weights=w*(1+.2*x), mask=mask,
                             removed_terms=[(2, 0, '0')])
        coeff = dict(zip(result['terms'], result['coefficients_mm']))
        self.assertEqual(result['status'], 'ok')
        self.assertAlmostEqual(coeff[(0, 0, '0')], .0007, places=14)
        self.assertAlmostEqual(coeff[(1, 1, 'c')], .0001, places=14)
        self.assertAlmostEqual(coeff[(2, 0, '0')], .0003, places=14)
        self.assertAlmostEqual(coeff[(3, 1, 's')], .0001, places=14)
        self.assertLess(result['residual_rms_mm'], 1e-17)
        assert_allclose(opd, original, atol=0, rtol=0)
        assert_allclose(result['diagnostic_opd_mm'], opd[mask]-.0003*np.sqrt(3)*(2*r2[mask]-1), atol=1e-17)
        self.assertFalse(np.shares_memory(result['diagnostic_opd_mm'], opd))

    def test_rank_deficiency_explicit_and_wave_units_configurable(self):
        xy = np.zeros((8, 2))
        result = zernike_fit(xy, np.ones(8)*.0005, n_max=2, wavelength_nm=500)
        self.assertEqual(result['status'], 'fit_rank_deficient')
        self.assertEqual(result['rank'], 1)
        self.assertTrue(np.isinf(result['condition_number']))
        assert_allclose(result['coefficients_waves'], result['coefficients_mm']/.0005)

    def test_invalid_chart_not_silently_rescaled(self):
        with self.assertRaises(InputError):
            zernike_fit([[1.1, 0]], [0.], n_max=0)


class MappingFanTests(unittest.TestCase):
    def test_radial_cubic_distortion_reference_and_axis_limit(self):
        q = np.array([[0., 0.], [1., 0.], [1., 2.]])
        g = 3*q*(1+.02*np.sum(q*q, axis=1))[:, None]
        result = map_distortion(q, g, 3*np.eye(2))
        assert_allclose(result['radial_percent'], [0, 2, 10])
        self.assertTrue(result['origin_limit_used'][0])

    def test_refitted_distortion_cannot_hide_nominal_scale_drift(self):
        q = np.array([[0., 0.], [1., 2.], [-3., .5]])
        result = map_distortion(q, 1.01*q, 1.01*np.eye(2), nominal_mapping=q)
        assert_allclose(result['radial_percent'], 0, atol=1e-14)
        assert_allclose(result['nominal_change'], .01*q, atol=1e-15)
        assert_allclose(result['actual_mapping'], 1.01*q)

    def test_even_odd_coma_descriptor_preserves_cross_component(self):
        q = np.linspace(-1, 1, 11)
        v = np.column_stack([2*q*q+3*q+4*q**3, 7*q*q])
        result = coma_from_fan(q, v)
        assert_allclose(result['even'], np.column_stack([2*q*q, 7*q*q]), atol=2e-15)
        assert_allclose(result['odd'], np.column_stack([3*q+4*q**3, np.zeros_like(q)]), atol=2e-15)
        with self.assertRaises(InputError):
            coma_from_fan([0, .5, 1], [0, 1, 2])

    def test_parabasal_curvature_rotated_analytic_family(self):
        field = Field(.3, .7, 'height')
        t, s = tangential_sagittal_basis(field)
        basis = np.column_stack([t, s])
        system = SimpleNamespace(image_z=12., surfaces=[SimpleNamespace(medium_after='air', z=0.)],
                                 stop_semi_diameter=1., index=lambda m, w: 1.)
        def family(system, f, pupil, wl, cfg):
            a = basis.T @ pupil
            x = basis @ (np.array([2., 3.])*a)
            slope = basis @ (-np.array([2/10, 3/20])*a)
            return trace([*x, 0.], slope)
        with patch('sequential_optics.metrics._trace_at', side_effect=family):
            result = differential_field_curvature(system, field)
        self.assertEqual(result['status'], 'ok')
        assert_allclose(result['focus_z_TS_mm'], [10, 20], atol=2e-13)
        assert_allclose(result['field_curvature_TS_mm'], [-2, 8], atol=2e-13)
        assert_allclose(result['vergence_TS_per_m'], [-500, 125], atol=1e-10)

    def test_chief_mapping_small_field_jacobian_real_intercepts(self):
        system = SimpleNamespace(image_z=0.)
        def family(system, f, pupil, wl, cfg):
            q = np.tan([f.x, f.y]); pos = 3*q*(1+.02*np.sum(q*q))
            return trace([*pos, 0.])
        fields = [Field(), Field(np.arctan(1), 0), Field(np.arctan(1), np.arctan(2))]
        with patch('sequential_optics.metrics._trace_at', side_effect=family):
            result = distortion(system, fields)
        assert_allclose(result['radial_percent'], [0, 2, 10], atol=3e-10)
        assert_allclose(result['reference_jacobian'], 3*np.eye(2), atol=1e-11)


class RealLensRegressionTests(unittest.TestCase):
    @staticmethod
    def system():
        # Already debated sample: moderate aperture, unoptimized three lenses,
        # one supplied G index only and no thermal data. Not catalog truth.
        from examples.sample_systems import synthetic_triplet
        return synthetic_triplet()

    def test_axis_field_curvature_at_exact_paraxial_focus_is_well_conditioned(self):
        system = self.system()
        result = differential_field_curvature(system, Field())
        self.assertEqual(result['status'], 'ok')
        assert_allclose(result['focus_z_TS_mm'], [system.image_z, system.image_z], atol=2e-7)
        self.assertLess(max(result['derivative_relative_change']), 1e-7)

    def test_real_offaxis_curvature_rotates_with_field_meridian(self):
        system = self.system()
        tangent = np.array([.03, .02]); phi = .63
        r = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        a = differential_field_curvature(system, Field(*np.arctan(tangent)))
        b = differential_field_curvature(system, Field(*np.arctan(r @ tangent)))
        self.assertEqual(a['status'], 'ok'); self.assertEqual(b['status'], 'ok')
        assert_allclose(a['focus_z_TS_mm'], b['focus_z_TS_mm'], atol=2e-7)
        self.assertGreater(abs(a['astigmatic_focus_separation_mm']), .01)

    def test_real_axis_LSA_not_conflated_with_parabasal_focus(self):
        system = self.system()
        result = longitudinal_spherical_aberration(system, pupil_coordinates=[0, .001, .5, 1])
        self.assertEqual(result['statuses'][0], 'paraxial_limit')
        self.assertLess(abs(result['lsa_mm'][1]), 1e-6)
        self.assertGreater(abs(result['lsa_mm'][-1]), .01)
        self.assertTrue(np.all(np.isfinite(result['axis_crossing_z_mm'])))

    def test_real_axis_fan_odd_and_offaxis_coma_even(self):
        system = self.system()
        axial = ray_fans(system, Field(), points=9)
        paired = coma_from_fan(axial['pupil_q'], axial['T_fan_TS_mm'])
        assert_allclose(paired['even'], 0., atol=1e-10)
        assert_allclose(axial['T_fan_TS_mm'][:, 1], 0., atol=1e-10)
        offaxis = ray_fans(system, Field(.03, .02), points=9)
        coma = coma_from_fan(offaxis['pupil_q'], offaxis['T_fan_TS_mm'])
        self.assertGreater(np.nanmax(abs(coma['even'][:, 0])), 1e-5)

    def test_grid_uses_original_fields_and_same_image_plane(self):
        system = self.system()
        result = grid_mapping(system, [-.02, 0, .02], [-.01, 0, .01], field_kind='angle')
        self.assertEqual(result['mapped_grid'].shape, (3, 3, 2))
        self.assertEqual(result['image_z_mm'], system.image_z)
        assert_allclose(result['mapped_grid'][1, 1], [0, 0], atol=1e-12)
        assert_allclose(result['mapped_grid'][0, 0], -result['mapped_grid'][-1, -1], atol=1e-10)
        self.assertEqual(result['fields'][0], Field(-.02, -.01))

    def test_finite_object_collimator_angular_fan_needs_no_image_plane(self):
        from dataclasses import replace
        from examples.sample_systems import synthetic_collimator
        from sequential_optics.geometry import aim_ray, trace_ray
        system = synthetic_collimator()
        field = Field(.3, .1, 'height')
        fan = ray_fans(system, field, points=5, angular=True)
        self.assertIsNone(fan['image_z_mm'])
        self.assertEqual(fan['units'], 'rad')
        self.assertTrue(np.all(np.isfinite(fan['T_fan_TS_rad'])))
        shifted = ray_fans(replace(system, image_z=1000.), field, points=5, angular=True)
        assert_allclose(fan['T_fan_TS_rad'], shifted['T_fan_TS_rad'], atol=0, rtol=0)
        # Same one ray gives identical chart coordinates in angular spot and fan.
        pupil = fan['basis_T']
        chief = trace_ray(system, aim_ray(system, field, [0., 0.]))
        tr = trace_ray(system, aim_ray(system, field, pupil))
        rb = RayBundle((tr,), pupil[None, :], np.ones(1), field, 550., chief)
        expected = angular_spot(rb, reference='chief')['coordinates_rad'][0]
        assert_allclose(fan['T_fan_XY_rad'][-1], expected, atol=1e-15)
        d = tr.last_direction; dc = chief.last_direction
        expected_ts = np.arctan2([d@fan['angular_basis_T'], d@fan['angular_basis_S']], d@dc)
        assert_allclose(fan['T_fan_TS_rad'][-1], expected_ts, atol=1e-15)


if __name__ == '__main__':
    unittest.main()
