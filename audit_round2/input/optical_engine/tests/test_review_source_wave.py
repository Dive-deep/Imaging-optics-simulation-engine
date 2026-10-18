"""Independent review of the physical source-to-wave amplitude bridge."""
import unittest
from dataclasses import replace
import numpy as np
from sequential_optics.models import Surface, OpticalSystem, Field
from sequential_optics.materials import Material
from sequential_optics.settings import EngineSettings, SamplingSettings
from sequential_optics.geometry import trace_bundle
from sequential_optics.source_wave import source_diffraction_kernel, wave_to_source_kernel
from sequential_optics.imaging import PixelGrid, synthesize_image
from sequential_optics.errors import InputError


SETTINGS=EngineSettings(sampling=SamplingSettings(pupil_rings=5,pupil_arms=48,wave_grid_size=64,fft_padding=2))


def transparent_aperture(index=1.):
    return OpticalSystem((Surface(0.,0.,.1),),{'air':Material('homogeneous',n_constant=index)},
                         object_z=-100.,image_z=20.)


class IndependentSourceWaveReviewTests(unittest.TestCase):
    def test_absolute_n_changes_electric_amplitude_not_collected_packet_power(self):
        args=dict(source_xy_mm=[0,0],source_area_mm2=.003,source_radiance=2.,
                  fill_factor=.64,settings=SETTINGS,transmission_mode='lossless_scalar',mode='afocal')
        a=source_diffraction_kernel(transparent_aperture(1.),**args)
        b=source_diffraction_kernel(transparent_aperture(1.44),**args)
        self.assertAlmostEqual(a['collected_power_density'],b['collected_power_density'],places=15)
        np.testing.assert_allclose(b['boundary_triangle_amplitude_v_m']/a['boundary_triangle_amplitude_v_m'],
                                   1/np.sqrt(1.44),rtol=2e-12)
        # This is compared BEFORE PSF normalization to an independent analytic
        # Lambertian circular-aperture flux, not normalized Q*h closure.
        expected=2*.64*.003*np.pi*.1**2/(100**2+.1**2)
        for result in (a,b):
            diagnostic=result['raw_power_diagnostics']
            self.assertLess(abs(diagnostic['boundary_grid_power_w']/expected-1),.025)
            self.assertLess(abs(diagnostic['wave_power_w_unvalidated']/expected-1),.025)
            self.assertEqual(diagnostic['raw_wave_phase_sampling_status'],'necessary_check_passed')

    def test_absorbing_index_matched_plate_is_applied_once(self):
        surfaces=(Surface(0.,0.,.2,medium_after='plate'),Surface(1.,0.,.1))
        def prescription(mu):
            return OpticalSystem(surfaces,{'air':Material('air',n_constant=1.,absorption_constant_per_mm=0.),
                  'plate':Material('matched absorbing medium',n_constant=1.,absorption_constant_per_mm=mu)},
                  object_z=-10.,image_z=12.,stop_index=1,stop_radius=.01)
        grid=PixelGrid(np.linspace(-1,1,9),np.linspace(-1,1,9))
        args=dict(source_xy_mm=[0,0],source_area_mm2=.001,source_radiance=1.,settings=SETTINGS,image_grid=grid)
        clear=source_diffraction_kernel(prescription(0.),**args)
        half=source_diffraction_kernel(prescription(np.log(2.)),**args)
        self.assertAlmostEqual(half['collected_power_density']/clear['collected_power_density'],.5,delta=3e-7)
        # No additional attenuation is permitted when the matched kernel is
        # placed on the detector and multiplied by its already attenuated Q.
        images=[]
        for result in (clear,half):
            images.append(synthesize_image([result['kernel']],[result['collected_power_density']],[.001]))
        self.assertAlmostEqual(images[1]['detected_roi_power']/images[0]['detected_roi_power'],.5,delta=2e-6)

    def test_chief_parallel_projection_registers_global_detector_position(self):
        theta=.4;phi=.7
        e3=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        e1=np.array([1.,0.,0.])-e3[0]*e3;e1/=np.linalg.norm(e1)
        e2=np.cross(e3,e1);basis=np.array([e1,e2,e3])
        q0=np.array([1.2,-.7,20.]);u=.2;v=.1
        virtual=q0+u*e1+v*e2
        expected=(virtual+(q0[2]-virtual[2])/e3[2]*e3)[:2]
        psf=np.zeros((3,3));psf[1,1]=100.
        fixture={'wavelength_nm':550.,'wave':{'mode':'focal','illumination':'physical_source_triangle_boundary_amplitude',
            'output_basis':basis,'image_reference_point_mm':q0,'dx_mm':.1,'dy_mm':.1,
            'x_mm':np.array([.1,.2,.3]),'y_mm':np.array([0.,.1,.2]),'psf':psf}}
        grid=PixelGrid(np.linspace(expected[0]-.2,expected[0]+.2,81),np.linspace(expected[1]-.2,expected[1]+.2,81))
        kernel,note=wave_to_source_kernel(fixture,grid,allow_local_detector_projection=True)
        xx,yy=np.meshgrid(grid.x_centers_mm,grid.y_centers_mm)
        centroid=np.array([np.sum(xx*kernel.probability_mass),np.sum(yy*kernel.probability_mass)])
        self.assertAlmostEqual(kernel.captured_fraction,1.,places=11)
        np.testing.assert_allclose(centroid,expected,atol=2e-5)
        self.assertTrue(note['approximation'])

    def test_same_field_old_bundle_cannot_change_physical_source_z(self):
        original=transparent_aperture()
        bundle=trace_bundle(original,Field(kind='height'),settings=SETTINGS)
        changed=replace(original,object_z=-200.)
        with self.assertRaisesRegex(InputError,'physical source point'):
            source_diffraction_kernel(changed,[0,0],.001,1.,bundle=bundle,settings=SETTINGS,
                                      transmission_mode='lossless_scalar',mode='afocal')


if __name__=='__main__': unittest.main()
