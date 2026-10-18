"""Matched physical-source radiometry/diffraction integration fixtures."""
import unittest
import numpy as np
from dataclasses import replace
from scipy.integrate import dblquad
from sequential_optics.models import Surface,OpticalSystem,Field
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_bundle
from sequential_optics.settings import EngineSettings,SamplingSettings
from sequential_optics.source_wave import (source_diffraction_kernel,spherical_triangle_measures,
                                         wave_to_source_kernel,bin_source_angular_response)
from sequential_optics.imaging import PixelGrid,synthesize_image
from sequential_optics.errors import InputError,ModelApplicabilityError


def aperture(radius=.1):
    return OpticalSystem((Surface(0.,0.,radius),),
             {'air':Material('G only',wavelengths_nm=(550.,),indices=(1.,))},
             object_z=-10.,image_z=10.)


class SourceWaveTests(unittest.TestCase):
    def setUp(self):
        self.settings=EngineSettings(sampling=SamplingSettings(pupil_rings=5,pupil_arms=32,wave_grid_size=32,fft_padding=2))

    def test_spherical_projected_measure_independent_surface_integral(self):
        directions=np.array([[0.,0.,1.],[1.,0.,1.],[0.,1.,1.]])
        omega,projected=spherical_triangle_measures(directions)
        expected=dblquad(lambda y,x:1/(1+x*x+y*y)**2,0,1,lambda x:0,lambda x:1-x)[0]
        self.assertAlmostEqual(projected,expected,places=12)
        self.assertGreater(omega,projected)
        tiny=np.array([[0.,0.,1.],[1e-6,0.,1.],[0.,1e-6,1.]])
        self.assertAlmostEqual(spherical_triangle_measures(tiny)[1]/5e-13,1.,places=10)
        self.assertAlmostEqual(spherical_triangle_measures(directions[::-1])[1],-projected,places=12)

    def test_lambertian_aperture_Q_converges_with_boundary_polygon(self):
        system=aperture();expected=np.pi*.1**2/(10**2+.1**2)
        errors=[]
        for arms in (12,48):
            bundle=trace_bundle(system,Field(kind='height'),settings=self.settings,rings=5,arms=arms)
            r=source_diffraction_kernel(system,[0,0],.001,1.,settings=self.settings,bundle=bundle,
                                       transmission_mode='lossless_scalar',mode='afocal')
            errors.append(abs(r['collected_power_density']-expected))
            self.assertLess(abs(r['raw_power_diagnostics']['boundary_grid_to_ray_ratio']-1),.1)
        self.assertLess(errors[1],errors[0]/10)
        self.assertLess(errors[1]/expected,.003)

    def test_radiance_fill_and_source_area_do_not_change_normalized_psf(self):
        system=aperture();bundle=trace_bundle(system,Field(kind='height'),settings=self.settings)
        args=dict(system=system,source_xy_mm=[0,0],settings=self.settings,bundle=bundle,transmission_mode='lossless_scalar')
        a=source_diffraction_kernel(**args,source_area_mm2=.001,source_radiance=1.,fill_factor=.5)
        b=source_diffraction_kernel(**args,source_area_mm2=.003,source_radiance=4.,fill_factor=.5)
        self.assertAlmostEqual(b['collected_power_density']/a['collected_power_density'],4.,places=12)
        self.assertAlmostEqual(b['source_cell_collected_power']/a['source_cell_collected_power'],12.,places=12)
        self.assertTrue(np.allclose(b['boundary_triangle_amplitude_v_m']/a['boundary_triangle_amplitude_v_m'],np.sqrt(12)))
        self.assertTrue(np.allclose(a['wave']['psf'],b['wave']['psf'],rtol=1e-12,atol=1e-9))
        self.assertAlmostEqual(b['raw_power_diagnostics']['wave_power_w_unvalidated']/a['raw_power_diagnostics']['wave_power_w_unvalidated'],12.,places=10)

    def test_matched_Q_kernel_image_power_and_crop(self):
        system=aperture();grid=PixelGrid(np.linspace(-1,1,17),np.linspace(-1,1,17))
        result=source_diffraction_kernel(system,[0,0],.001,1.,settings=self.settings,
                         transmission_mode='lossless_scalar',image_grid=grid,direct_check_points=[[0.,0.]])
        self.assertIn('local_expansion_diagnostics',result['wave'])
        direct=result['wave']['direct_scalar_check']
        self.assertTrue(np.allclose(direct['exact_intensity'],direct['local_fft_integral_intensity'],rtol=1e-7))
        kernel=result['kernel'];self.assertGreater(kernel.captured_fraction,.99)
        image=synthesize_image([kernel],[result['collected_power_density']],[.001])
        self.assertAlmostEqual(image['collected_power'],result['source_cell_collected_power'],places=15)
        self.assertAlmostEqual(image['detected_roi_power'],kernel.captured_fraction*result['source_cell_collected_power'],places=15)
        crop=PixelGrid(np.linspace(-.01,.01,5),np.linspace(-.01,.01,5))
        small,_=wave_to_source_kernel(result,crop)
        self.assertLess(small.captured_fraction,kernel.captured_fraction*.5)
        self.assertGreater(small.captured_fraction,0.)

    def test_angular_bin_rotation_retains_cell_power(self):
        system=aperture();result=source_diffraction_kernel(system,[.04,.01],.001,1.,settings=self.settings,
                                      mode='afocal',transmission_mode='lossless_scalar')
        b=bin_source_angular_response(result,np.linspace(0,np.pi,20),np.linspace(-np.pi,np.pi,24))
        self.assertAlmostEqual(b['bin_power'].sum(),result['source_cell_collected_power'],places=15)
        self.assertAlmostEqual(np.sum(b['intensity']*b['solid_angle_sr']),result['source_cell_collected_power'],places=15)
        self.assertEqual(b['registration'],'full_unit_direction_rotation')

    def test_offaxis_projection_requires_explicit_model_choice(self):
        result=source_diffraction_kernel(aperture(),[.2,0],.001,1.,settings=self.settings,transmission_mode='lossless_scalar')
        grid=PixelGrid(np.linspace(-1,1,17),np.linspace(-1,1,17))
        with self.assertRaises(ModelApplicabilityError): wave_to_source_kernel(result,grid)
        kernel,note=wave_to_source_kernel(result,grid,allow_local_detector_projection=True)
        self.assertTrue(note['approximation']);self.assertGreater(kernel.captured_fraction,.99)

    def test_source_bundle_mismatch_and_absorbing_output_rejected(self):
        system=aperture();bundle=trace_bundle(system,Field(.1,0,'height'),settings=self.settings)
        with self.assertRaises(InputError):
            source_diffraction_kernel(system,[0,0],.001,1.,settings=self.settings,bundle=bundle)
        absorbing=replace(system,materials={'air':Material('absorbing',n_constant=1.,absorption_constant_per_mm=.01)})
        with self.assertRaises(ModelApplicabilityError):
            source_diffraction_kernel(absorbing,[0,0],.001,1.,settings=self.settings)

if __name__=='__main__': unittest.main()
