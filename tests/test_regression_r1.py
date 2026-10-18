"""R1: requested angular frequency must be resolved by angular cell quadrature."""
import unittest
from dataclasses import replace
import numpy as np
from sequential_optics.wave import afocal_fft, angular_otf, analyze_wave
from sequential_optics.settings import EngineSettings, SamplingSettings
from sequential_optics.models import OpticalSystem, Surface, Field
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_bundle
from sequential_optics.source_wave import source_diffraction_kernel
from sequential_optics.errors import ModelApplicabilityError, InputError


class AngularSamplingRegression(unittest.TestCase):
    def gaussian(self, padding):
        n=64; spacing=.01; u=(np.arange(n)-n//2)*spacing
        x,y=np.meshgrid(u,u)
        return afocal_fft(np.exp(-(x*x+y*y)/.06**2),spacing,.00055,padding=padding)

    def test_alias_is_rejected_and_dc_preserved(self):
        result=self.gaussian(4); alias=64*4*.01/.00055
        with self.assertRaisesRegex(ModelApplicabilityError,'underresolved'):
            angular_otf(result,[alias])
        checked=angular_otf(result,[0.,100.,alias,-alias],return_diagnostics=True)
        np.testing.assert_array_equal(checked['resolved_mask'],[True,True,False,False])
        self.assertAlmostEqual(checked['otf'][0].real,1.,places=12)
        self.assertTrue(np.isnan(checked['otf'][2:]).all())
        self.assertGreater(checked['unresolved_cell_power_fraction'][2],.99)

    def test_padding_refinement_removes_false_alias_and_converges(self):
        alias=64*4*.01/.00055
        fine=angular_otf(self.gaussian(16),[100.,alias],return_diagnostics=True)
        finer=angular_otf(self.gaussian(24),[100.,alias],return_diagnostics=True)
        self.assertTrue(fine['resolved_mask'].all())
        np.testing.assert_allclose(fine['otf'],finer['otf'],atol=1e-8)
        self.assertLess(abs(fine['otf'][1]),1e-7)
        self.assertFalse(fine['convergence_verified'])

    def test_exact_angle_span_not_direction_cosine_cutoff(self):
        r=self.gaussian(4)
        # Deliberately high angle: same cell width in l is a larger angle width.
        n=8; l=np.linspace(.78,.92,n); m=np.linspace(-.03,.03,n)
        ll,mm=np.meshgrid(l,m);g=np.sqrt(1-ll*ll-mm*mm)
        r.update(l=l,m=m,dl=l[1]-l[0],dm=m[1]-m[0],
                 alpha_rad=np.arctan2(ll,g),beta_rad=np.arctan2(mm,g),
                 cell_power=np.full((n,n),1/n**2))
        f=.45/r['dl'] # passes the direction-cosine Nyquist test
        d=angular_otf(r,[f],return_diagnostics=True)
        self.assertFalse(d['resolved_mask'][0])
        self.assertGreater(d['unresolved_cell_power_fraction'][0],.1)

    def test_high_level_lens_and_reference_report_underresolved(self):
        system=OpticalSystem((Surface(0.,0.,.3),),{'air':Material('air',n_constant=1.)})
        cfg=EngineSettings(sampling=SamplingSettings(pupil_rings=5,pupil_arms=24,wave_grid_size=32))
        bundle=trace_bundle(system,Field(),settings=cfg)
        amp=np.exp(-np.sum((.3*bundle.pupil_xy)**2,axis=1)/.06**2)
        r=analyze_wave(system,bundle,cfg,mode='afocal',source_boundary_amplitudes=amp,
                       frequencies=[0.,100.,4*(2*.3*1.05)/.00055])
        self.assertEqual(r['status'],'underresolved')
        self.assertTrue(r['sampling_diagnostics']['phase_sampling_necessary_check_passed'])
        self.assertFalse(r['sampling_diagnostics']['angular_frequency_sampling_necessary_check_passed'])
        self.assertTrue(np.isnan(r['mtf_s'][-1]))
        self.assertTrue(np.isnan(r['diffraction_mtf_s'][-1]))
        self.assertAlmostEqual(r['mtf_s'][0],1.,places=12)

    def test_physical_source_path_uses_same_frequency_guard(self):
        system=OpticalSystem((Surface(0.,0.,.1),),{'air':Material('air',n_constant=1.)},object_z=-10.)
        cfg=EngineSettings(sampling=SamplingSettings(pupil_rings=3,pupil_arms=16,wave_grid_size=16,fft_padding=2))
        r=source_diffraction_kernel(system,[0,0],.001,1.,settings=cfg,mode='afocal',
                                   transmission_mode='lossless_scalar',frequencies=[0.,1e8])['wave']
        self.assertEqual(r['status'],'underresolved')
        self.assertFalse(r['angular_frequency_sampling']['necessary_check_passed'])
        self.assertTrue(np.isnan(r['mtf_s'][-1]))

    def test_sampling_limits_live_in_settings(self):
        with self.assertRaises(InputError):
            EngineSettings(sampling=SamplingSettings(angular_otf_phase_step_limit_rad=4.))
        cfg=EngineSettings(sampling=SamplingSettings(angular_otf_phase_step_limit_rad=.01))
        d=angular_otf(self.gaussian(4),[100.],settings=cfg,return_diagnostics=True)
        self.assertFalse(d['resolved_mask'][0])


if __name__=='__main__': unittest.main()
