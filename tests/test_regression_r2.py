"""R2: trace-time provenance must precede construction of every wave cache.

The immutable round-2 audit remains untouched. These assertions expect the
correct behavior of the working root source, rather than replaying old bugs.
"""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import sys
import unittest
import numpy as np

from sequential_optics.errors import InputError
from sequential_optics.geometry import trace_bundle
from sequential_optics.materials import Material
from sequential_optics.models import Field, OpticalSystem, Surface
from sequential_optics.provenance import require_bundle_match
from sequential_optics.settings import EngineSettings, SamplingSettings
from sequential_optics.source_wave import source_diffraction_kernel
from sequential_optics.wave import analyze_wave, make_boundary_state

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
from sample_systems import synthetic_triplet


class TraceProvenanceRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.settings = EngineSettings(sampling=SamplingSettings(
            pupil_rings=4, pupil_arms=16, wave_grid_size=32, fft_padding=4))
        air = Material('monochromatic synthetic air', wavelengths_nm=(550.,), indices=(1.,),
                       index_temperature_invariant=True)
        cls.system = OpticalSystem((Surface(0., 0., 2.), Surface(1., 0., 2.)),
                                   {'air': air}, stop_index=1, image_z=40.)
        cls.bundle = trace_bundle(cls.system, Field(), settings=cls.settings)

    def assertStaleRejected(self, system, bundle=None, settings=None, **kwargs):
        with self.assertRaisesRegex(InputError, 'provenance mismatch'):
            analyze_wave(system, self.bundle if bundle is None else bundle,
                         self.settings if settings is None else settings,
                         frequencies=[0., 1.], **kwargs)

    def test_round2_last_surface_change_requires_retrace(self):
        system = synthetic_triplet()
        bundle = trace_bundle(system, Field(), settings=self.settings)
        baseline = analyze_wave(system, bundle, self.settings, frequencies=[0., 60.])
        surfaces = (*system.surfaces[:-1], replace(system.surfaces[-1],
                    curvature=system.surfaces[-1].curvature * .995))
        changed = replace(system, surfaces=surfaces)
        self.assertStaleRejected(changed, bundle)
        self.assertStaleRejected(changed, bundle, boundary_state=baseline['boundary_state'])
        fresh = analyze_wave(changed, trace_bundle(changed, Field(), settings=self.settings),
                             self.settings, frequencies=[0., 60.])
        self.assertEqual(fresh['trace_provenance_status'], 'verified_engine_trace')
        # Independent optical evidence that this is a materially changed input.
        self.assertGreater(abs(fresh['mtf_s'][1] - baseline['mtf_s'][1]), .05)

    def test_stop_radius_and_index_are_trace_inputs_even_with_cached_boundary(self):
        state = make_boundary_state(self.system, self.bundle, self.settings)
        for changed in (replace(self.system, stop_radius=1.), replace(self.system, stop_index=0)):
            with self.subTest(stop_index=changed.stop_index, radius=changed.stop_radius):
                self.assertStaleRejected(changed)
                self.assertStaleRejected(changed, boundary_state=state)
                with self.assertRaisesRegex(InputError, 'provenance mismatch'):
                    make_boundary_state(changed, self.bundle, self.settings)

    def test_temperature_material_conjugate_and_apertures_are_trace_inputs(self):
        modified_surfaces = (*self.system.surfaces[:-1], replace(self.system.surfaces[-1],
                             mechanical_radius=2.1))
        changed_index = replace(self.system, materials={'air': Material('changed G',
                                wavelengths_nm=(550.,), indices=(1.01,))})
        systems = (replace(self.system, temperature_c=21.), changed_index,
                   replace(self.system, object_z=-100.),
                   replace(self.system, surfaces=modified_surfaces))
        for changed in systems:
            with self.subTest(changed=changed):
                self.assertStaleRejected(changed)

    def test_field_wavelength_weights_and_intermediate_ray_mutation_are_rejected(self):
        mutated = deepcopy(self.bundle)
        mutated.traces[0].positions[0, 0] += .01  # final hit and direction unchanged
        changed_weight = self.bundle.weights.copy(); changed_weight[0] *= 1.1
        bundles = (replace(self.bundle, field=Field(.01, 0.)),
                   replace(self.bundle, wavelength_nm=650.),
                   replace(self.bundle, weights=changed_weight), mutated)
        for changed in bundles:
            with self.subTest(kind=repr(changed.field), wavelength=changed.wavelength_nm):
                self.assertStaleRejected(self.system, changed)

    def test_ray_solver_settings_are_preserved_but_fft_sampling_may_change(self):
        changed_ray = replace(self.settings, ray=replace(self.settings.ray, aim_tolerance_mm=1e-8))
        self.assertStaleRejected(self.system, settings=changed_ray)
        changed_launch = replace(self.settings, sampling=replace(self.settings.sampling,
                                 boundary_clearance_mm=.002))
        self.assertStaleRejected(self.system, settings=changed_launch)
        # Explicit sample counts are represented by actual bundle arrays, so
        # requested FFT/reported frequencies or future ray defaults do not
        # falsely invalidate already-traced samples.
        wave_only = replace(self.settings, sampling=replace(self.settings.sampling,
                            pupil_rings=9, wave_grid_size=48, fft_padding=2))
        result = analyze_wave(self.system, self.bundle, wave_only, frequencies=[0., 1.])
        self.assertEqual(result['trace_provenance_status'], 'verified_engine_trace')

    def test_ima_only_motion_explicitly_reuses_fixed_boundary(self):
        state = make_boundary_state(self.system, self.bundle, self.settings)
        moved = replace(self.system, image_z=self.system.image_z + .5)
        a = analyze_wave(moved, self.bundle, self.settings, boundary_state=state,
                         frequencies=[0., 1.])
        b = analyze_wave(self.system, self.bundle, self.settings, image_z=moved.image_z,
                         boundary_state=state, frequencies=[0., 1.])
        self.assertIs(a['boundary_state'], state)
        np.testing.assert_allclose(a['mtf_s'], b['mtf_s'], rtol=0, atol=0)

    def test_source_bridge_rejects_changed_later_interface(self):
        system = replace(self.system, object_z=-100.)
        bundle = trace_bundle(system, Field(kind='height'), settings=self.settings)
        changed = replace(system, surfaces=(system.surfaces[0],
                          replace(system.surfaces[1], curvature=.001)))
        with self.assertRaisesRegex(InputError, 'provenance mismatch'):
            source_diffraction_kernel(changed, [0., 0.], .001, 1., bundle=bundle,
                settings=self.settings, mode='afocal', transmission_mode='lossless_scalar')

    def test_manual_boundaries_require_opt_in_and_cannot_override_engine_mismatch(self):
        manual = replace(self.bundle, provenance=None)
        with self.assertRaisesRegex(InputError, 'no trace provenance'):
            make_boundary_state(self.system, manual, self.settings)
        result = analyze_wave(self.system, manual, self.settings, mode='afocal',
                              frequencies=[0.], allow_unverified_bundle=True)
        self.assertEqual(result['trace_provenance_status'], 'unverified_external_bundle')
        self.assertEqual(result['boundary_state'].trace_provenance_status, 'unverified_external_bundle')
        with self.assertRaisesRegex(InputError, 'no trace provenance'):
            analyze_wave(self.system, manual, self.settings, mode='afocal',
                         boundary_state=result['boundary_state'])
        self.assertStaleRejected(replace(self.system, stop_radius=1.),
                                 allow_unverified_bundle=True)

    def test_source_amplitude_cache_owns_snapshot_and_detects_replacement(self):
        amplitudes = np.ones(len(self.bundle.traces))
        state = make_boundary_state(self.system, self.bundle, self.settings,
                                    source_boundary_amplitudes=amplitudes)
        amplitudes *= 2
        np.testing.assert_array_equal(state.source_amplitudes, 1.)
        self.assertFalse(state.source_amplitudes.flags.writeable)
        with self.assertRaisesRegex(InputError, 'Boundary snapshot'):
            analyze_wave(self.system, self.bundle, self.settings,
                         boundary_state=replace(state, source_amplitudes=amplitudes),
                         frequencies=[0., 1.])

    def test_new_unused_cte_absorption_does_not_require_retracing_geometric_rays(self):
        # Trace identity reads only data used by this monochromatic calculation.
        air = replace(self.system.materials['air'], cte_per_k=0., absorption_constant_per_mm=.01)
        updated = replace(self.system, materials={'air': air})
        self.assertEqual(require_bundle_match(updated, self.bundle, self.settings),
                         'verified_engine_trace')

    def test_boundary_geometry_is_immutable_and_replace_does_not_reseal(self):
        state = make_boundary_state(self.system, self.bundle, self.settings)
        for name in ('origin', 'basis', 'uv_mm', 'opl_mm', 'triangles', 'valid_ray_mask',
                     'pupil_xy', 'sampled_jacobians_mm2', 'chief_position', 'chief_direction'):
            with self.subTest(array=name):
                with self.assertRaises(ValueError):
                    getattr(state, name).setflags(write=True)
        opl = state.opl_mm.copy(); opl[0] += .0001
        changed_states = (replace(state, opl_mm=opl),
                          replace(state, reference_distance_mm=state.reference_distance_mm + 1.),
                          replace(state, output_index=1.1))
        for changed in changed_states:
            with self.subTest(state=changed.content_fingerprint):
                with self.assertRaisesRegex(InputError, 'Boundary snapshot contents changed'):
                    analyze_wave(self.system, self.bundle, self.settings,
                                 boundary_state=changed, frequencies=[0., 1.])

    def test_physical_source_recomputes_loss_amplitudes_on_valid_ray_reuse(self):
        air = Material('air', n_constant=1., absorption_constant_per_mm=0.)
        plate = Material('index-matched plate', n_constant=1., absorption_constant_per_mm=0.)
        system = OpticalSystem((Surface(0., 0., .1, medium_after='plate'),
                                Surface(1., 0., .1)), {'air': air, 'plate': plate},
                               object_z=-100., image_z=40., stop_index=1)
        bundle = trace_bundle(system, Field(kind='height'), settings=self.settings)
        changed = replace(system, materials={'air': air,
                          'plate': replace(plate, absorption_constant_per_mm=np.log(2.))})
        args = dict(source_xy_mm=[0., 0.], source_area_mm2=.001, source_radiance=1.,
                    bundle=bundle, settings=self.settings, mode='afocal')
        clear = source_diffraction_kernel(system, **args)
        absorbing = source_diffraction_kernel(changed, **args)
        self.assertAlmostEqual(absorbing['collected_power_density'] / clear['collected_power_density'],
                               .5, delta=1e-6)


if __name__ == '__main__':
    unittest.main()
