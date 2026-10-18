"""Adversarial supplied-data and uniform-thermal invariants.

All glasses here are explicitly synthetic analytic specimens, not grades.
These tests verify mathematical invariants, not experimental material accuracy.
"""
import json
import math
import unittest
from dataclasses import replace

import numpy as np

from sequential_optics.errors import InputError, MissingDataError, GeometryError, ModelApplicabilityError
from sequential_optics.materials import Material
from sequential_optics.models import Surface, OpticalSystem
from sequential_optics.thermal import LensGroup, thermal_system, scale_surface, validate_radial_separation
from sequential_optics.geometry import sag
from sequential_optics.paraxial import paraxial_properties


def lens_system(material=None, two=False):
    air = Material('ideal vacuum-index air', n_constant=1., index_temperature_invariant=True,
                   absorption_constant_per_mm=0., absorption_temperature_invariant=True)
    glass = material or Material('synthetic G-only glass', wavelengths_nm=(550.,), indices=(1.5,))
    surfaces = [Surface(0., .02, 3., medium_after='g'), Surface(5., -.02, 3., medium_after='air')]
    materials = {'air': air, 'g': glass}
    if two:
        materials['h'] = Material('synthetic second glass', n_constant=1.6, cte_per_k=2e-4, index_temperature_invariant=True)
        surfaces += [Surface(8., .015, 3., medium_after='h'), Surface(12., -.015, 3., medium_after='air')]
    return OpticalSystem(tuple(surfaces), materials, object_z=-100., image_z=55., stop_radius=2.)


class MaterialDataTests(unittest.TestCase):
    def test_single_green_index_needs_no_other_data(self):
        m = Material('measured only G', wavelengths_nm=(550.,), indices=(1.52,))
        self.assertEqual(m.index(550.), 1.52)
        self.assertEqual(m.index(550., 20.), 1.52)
        with self.assertRaises(MissingDataError): m.index(650.)
        with self.assertRaises(MissingDataError): m.index(550., 21.)
        with self.assertRaises(MissingDataError): m.absorption(550.)
        self.assertEqual(m.expansion_scale(20., 20.), 1.)

    def test_constant_index_is_explicit_achromatic_not_thermal(self):
        m = Material('explicit achromatic', n_constant=1.5)
        self.assertEqual(m.index(450.), m.index(650.))
        with self.assertRaises(MissingDataError): m.index(550., 40.)
        m = replace(m, index_temperature_invariant=True)
        self.assertEqual(m.index(450., 40.), 1.5)

    def test_index_and_geometry_require_only_relevant_thermal_data(self):
        thermooptic = Material('index only', n_constant=1.5, dn_dT_per_k=2e-5)
        self.assertAlmostEqual(thermooptic.index(550., 70.), 1.501)
        with self.assertRaises(MissingDataError): thermooptic.expansion_scale(70.)
        expansion = Material('geometry only', wavelengths_nm=(550.,), indices=(1.5,), cte_per_k=1e-5)
        self.assertAlmostEqual(expansion.expansion_scale(70.), 1.0005)
        with self.assertRaises(MissingDataError): expansion.index(550., 70.)

    def test_temperature_wavelength_bilinear_table(self):
        m = Material('two dimensional data', wavelengths_nm=(500., 600.), temperatures_c=(0., 100.),
                     temperature_indices=((1.6, 1.5), (1.7, 1.6)))
        self.assertAlmostEqual(m.index(550., 25.), 1.575)
        with self.assertRaises(MissingDataError): m.index(550., 101.)
        with self.assertRaises(MissingDataError): m.index(499., 25.)

    def test_vacuum_wavelength_and_relative_index_independent(self):
        a = 1.0003
        m = Material('relative with vacuum axis', wavelengths_nm=(500., 600.), indices=(1.6, 1.5),
                     index_basis='reference_air', wavelength_basis='vacuum', reference_air_index=a)
        self.assertAlmostEqual(m.index(550.), 1.55 * a)
        both = replace(m, wavelength_basis='reference_air')
        self.assertAlmostEqual(both.index(550. * a), 1.55 * a)
        abs_air_axis = replace(both, index_basis='absolute')
        self.assertAlmostEqual(abs_air_axis.index(550. * a), 1.55)
        self.assertNotAlmostEqual(m.index(550.), both.index(550.), places=8)

    def test_missing_reference_air_and_relative_thermal_table_rejected(self):
        with self.assertRaises(MissingDataError): Material('ambiguous', n_constant=1.5, index_basis='reference_air')
        with self.assertRaises(InputError):
            Material('ambiguous T table', wavelengths_nm=(550.,), temperatures_c=(20.,),
                     temperature_indices=((1.5,),), index_basis='reference_air', reference_air_index=1.0003)

    def test_absolute_dn_dT_not_multiplied_by_reference_air_again(self):
        m = Material('relative reference, absolute derivative', n_constant=1.5, index_basis='reference_air',
                     reference_air_index=1.0003, dn_dT_per_k=1e-4)
        self.assertAlmostEqual(m.index(550., 40.), 1.5 * 1.0003 + .002)

    def test_spectral_dn_dT_table(self):
        m = Material('spectral derivative', wavelengths_nm=(500., 600.), indices=(1.6, 1.5),
                     dn_dT_values_per_k=(1e-5, 3e-5))
        self.assertAlmostEqual(m.index(550., 70.), 1.551)

    def test_sellmeier_known_synthetic_identity_and_pole(self):
        m = Material('synthetic analytic', sellmeier_b=(1.25,), sellmeier_c_um2=(0.,))
        self.assertAlmostEqual(m.index(550.), 1.5)
        pole = replace(m, sellmeier_c_um2=(.25,))
        with self.assertRaises(InputError): pole.index(500.)

    def test_model_ranges_and_bad_inputs(self):
        m = Material('bounded', n_constant=1.5, wavelength_range_nm=(450., 650.), temperature_range_c=(0., 85.))
        with self.assertRaises(MissingDataError): m.index(700.)
        with self.assertRaises(MissingDataError): m.index(550., 100.)
        for kwargs in ({'n_constant': 1.5, 'indices': (1.5,), 'wavelengths_nm': (550.,)},
                       {'n_constant': float('nan')}, {'n_constant': -1},
                       {'wavelengths_nm': (600., 500.), 'indices': (1.5, 1.6)}):
            with self.assertRaises(InputError): Material('invalid', **kwargs)

    def test_transmission_uses_power_optical_path_length(self):
        m = Material.from_internal_transmission(name='synthetic absorbing glass', wavelength_nm=(500., 600.),
                                                transmission=(.81, 1.), thickness_mm=2., n_constant=1.5)
        mu = m.absorption(500.)
        self.assertAlmostEqual(math.exp(-mu * 4.), .81 ** 2)
        self.assertAlmostEqual(m.absorption(550.), mu / 2)
        with self.assertRaises(MissingDataError): m.absorption(500., 40.)
        self.assertEqual(replace(m, absorption_temperature_invariant=True).absorption(600., 40.), 0.)
        for tau in (0., -0.1, 1.01):
            with self.assertRaises(InputError):
                Material.from_internal_transmission(name='bad tau', wavelength_nm=(550.,), transmission=(tau,), thickness_mm=1., n_constant=1.5)

    def test_serialization_retains_conventions(self):
        m = Material('supplied', wavelengths_nm=(500., 600.), indices=(1.6, 1.5), dn_dT_values_per_k=(1e-5, 2e-5),
                     cte_per_k=1e-5, index_basis='reference_air', reference_air_index=1.0003, source='synthetic test')
        restored = Material.from_dict(json.loads(json.dumps(m.to_dict())))
        self.assertEqual(m, restored)
        self.assertEqual(m.index(550., 40.), restored.index(550., 40.))


class ThermalGeometryTests(unittest.TestCase):
    def test_reference_noop_needs_no_cte_or_dn_dT(self):
        base = lens_system()
        result = thermal_system(base, 20.)
        self.assertEqual(base.surfaces, result.surfaces)
        self.assertEqual(result.index('g', 550.), 1.5)
        self.assertAlmostEqual(paraxial_properties(base)['efl_image_mm'], paraxial_properties(result)['efl_image_mm'])

    def test_geometry_only_thermal_does_not_require_spectrum_or_dn_dT(self):
        base = lens_system(Material('single G expands', wavelengths_nm=(550.,), indices=(1.5,), cte_per_k=1e-5))
        result = thermal_system(base, 70.)
        self.assertAlmostEqual(result.surfaces[1].z - result.surfaces[0].z, 5.0025)
        with self.assertRaises(MissingDataError): result.index('g', 550.)

    def test_all_asphere_terms_obey_exact_similarity(self):
        for curvature in (0., .02):
            surface = Surface(0., curvature, 3., conic=-.8,
                              coefficients=(1e-3, 1e-5, -1e-7, 1e-9, 1e-11, -1e-13),
                              surface_type='even_asphere', mechanical_radius=4.)
            scaled = scale_surface(surface, 1.04)
            r = np.linspace(0., 3.12, 101)
            np.testing.assert_allclose(sag(scaled, r, 0.), 1.04 * sag(surface, r / 1.04, 0.), rtol=3e-15, atol=1e-16)
            self.assertEqual(scaled.conic, surface.conic)
            self.assertEqual(scaled.mechanical_radius, 4.16)

    def test_different_materials_and_anchor_gap_drift(self):
        m = Material('synthetic first glass', n_constant=1.5, cte_per_k=1e-4, index_temperature_invariant=True)
        base = lens_system(m, two=True)
        hot = thermal_system(base, 120.)
        self.assertAlmostEqual(hot.surfaces[0].z, -.025)
        self.assertAlmostEqual(hot.surfaces[1].z, 5.025)
        self.assertAlmostEqual(hot.surfaces[2].z, 7.96)
        self.assertAlmostEqual(hot.surfaces[3].z, 12.04)
        self.assertAlmostEqual(hot.surfaces[2].z - hot.surfaces[1].z, 2.935)
        front = thermal_system(base, 120., anchor_fraction=0.)
        self.assertEqual(front.surfaces[0].z, 0.)
        self.assertEqual(front.surfaces[2].z, 8.)
        self.assertNotEqual(front.surfaces[2].z - front.surfaces[1].z, hot.surfaces[2].z - hot.surfaces[1].z)

    def test_fixed_gaps_not_fixed_independent_anchors(self):
        base = lens_system(Material('expanding', n_constant=1.5, cte_per_k=1e-4), two=True)
        hot = thermal_system(base, 120., layout='fixed_gaps')
        self.assertAlmostEqual(hot.surfaces[2].z - hot.surfaces[1].z, 3.)
        groups = (LensGroup(0, 1), LensGroup(2, 3, anchor_shift_per_k_mm=.01))
        with self.assertRaises(InputError): thermal_system(base, 120., layout='fixed_gaps', lens_groups=groups)

    def test_stop_object_and_image_policies(self):
        base = lens_system(Material('expanding', n_constant=1.5, cte_per_k=1e-4))
        hot = thermal_system(base, 120.)
        self.assertAlmostEqual(hot.stop_radius, 2.02)
        self.assertEqual(hot.object_z, base.object_z)
        self.assertEqual(hot.image_z, base.image_z)
        followed = thermal_system(base, 120., stop_policy='fixed_radius', object_policy='first_surface_offset', image_policy='last_surface_offset')
        self.assertEqual(followed.stop_radius, 2.)
        self.assertAlmostEqual(followed.object_z, -100.025)
        self.assertAlmostEqual(followed.image_z, 55.025)

    def test_baseline_sweep_order_and_no_cumulative_drift(self):
        base = lens_system(Material('expanding', n_constant=1.5, cte_per_k=1e-4))
        high1 = thermal_system(base, 120.)
        thermal_system(base, -20.)
        high2 = thermal_system(base, 120.)
        self.assertEqual(high1.surfaces, high2.surfaces)
        self.assertEqual(thermal_system(base, 20.).surfaces, base.surfaces)
        with self.assertRaises(InputError): thermal_system(high1, 20.)

    def test_cte_reference_distinct_from_geometry_reference(self):
        m = Material('length origin at zero', n_constant=1.5, reference_temperature_c=10.,
                     cte_reference_temperature_c=0., cte_per_k=1e-3)
        self.assertAlmostEqual(m.expansion_scale(100., 20.), 1.1 / 1.02)
        self.assertAlmostEqual(m.expansion_scale(20., 100.) * m.expansion_scale(100., 20.), 1.)

    def test_cte_honors_declared_model_temperature_validity(self):
        m = Material('limited model', n_constant=1.5, cte_per_k=1e-5, temperature_range_c=(0., 85.))
        with self.assertRaises(MissingDataError): m.expansion_scale(200., 20.)
        with self.assertRaises(MissingDataError): m.expansion_scale(20., 200.)
        # Geometry identity does not ask the material to predict a change.
        self.assertEqual(m.expansion_scale(200., 200.), 1.)

    def test_cte_definitions_and_piecewise_integral(self):
        m = Material('varying instantaneous alpha', n_constant=1.5, cte_model='instantaneous',
                     cte_temperatures_c=(0., 100.), cte_values_per_k=(1e-4, 3e-4))
        self.assertAlmostEqual(m.expansion_scale(100., 0.), math.exp(.02))
        self.assertAlmostEqual(m.expansion_scale(0., 100.), math.exp(-.02))
        refderiv = replace(m, cte_model='reference_derivative', cte_reference_temperature_c=0.)
        self.assertAlmostEqual(refderiv.expansion_scale(100., 0.), 1.02)
        with self.assertRaises(MissingDataError): m.expansion_scale(101., 0.)
        with self.assertRaises(InputError): replace(m, cte_model='linear_reference')

    def test_thermal_focal_similarity_with_invariant_index(self):
        base = replace(lens_system(Material('constant-index expansion', n_constant=1.5, cte_per_k=1e-4,
                                            index_temperature_invariant=True)), object_z=None)
        hot = thermal_system(base, 120.)
        f0 = paraxial_properties(base)['efl_image_mm']
        f1 = paraxial_properties(hot)['efl_image_mm']
        self.assertAlmostEqual(f1 / f0, 1.01)

    def test_internal_collision_missed_by_center_and_edge_is_rejected(self):
        # Gap(r)=1 - 5r² + 5r⁴: positive at r=0 and 1, negative at r²=.5.
        a = Surface(0., 0., 1., medium_after='g')
        b = Surface(1., 0., 1., coefficients=(-5., 5.), surface_type='even_asphere', medium_after='air')
        base = replace(lens_system(), surfaces=(a, b), stop_radius=.5)
        self.assertGreater(b.z + sag(b, 0., 0.), 0)
        self.assertGreater(b.z + sag(b, 1., 0.), 0)
        with self.assertRaises(GeometryError): validate_radial_separation(base)

    def test_cte_nonpositive_length_and_bad_grouping(self):
        m = Material('invalid at extreme T', n_constant=1.5, cte_per_k=.1)
        with self.assertRaises(InputError): m.expansion_scale(-20., 20.)
        base = lens_system(Material('g', n_constant=1.5, cte_per_k=0.))
        with self.assertRaises(InputError): thermal_system(base, 50., lens_groups=(LensGroup(0, 1, material_id='air'),))
        with self.assertRaises(InputError): LensGroup(0, 2)


if __name__ == '__main__': unittest.main()
