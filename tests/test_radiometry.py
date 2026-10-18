"""Independent physical fixtures: measures, Fresnel, source/Q kernels.

Tests check analytic identities and conservation, not a Zemax reference or a
measured prescription. Finite image is geometric; source diffraction remains
an explicit bridge requiring source-conditioned coherent amplitudes.
"""
import unittest
from dataclasses import replace
import numpy as np
from scipy.integrate import dblquad
from sequential_optics.models import Surface,OpticalSystem,Ray
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_ray
from sequential_optics.errors import InputError,MissingDataError,ModelApplicabilityError
from sequential_optics.radiometry import (RectangularSource,fresnel_power,trace_transmission,
    far_field,uniform_scene_ri,spherical_bin_solid_angles,slope_rectangle_solid_angle,
    angle_chart_jacobian,projected_solid_angle_aperture,source_direction_samples)
from sequential_optics.imaging import (PixelGrid,SourceKernel,synthesize_image,
    geometric_source_kernel,geometric_image,registered_spectral_sum,combine_complex_otfs)


def plate(absorption=0.,air_absorption=0.):
    air=Material('air',n_constant=1.,absorption_constant_per_mm=air_absorption)
    glass=Material('glass',n_constant=1.5,absorption_constant_per_mm=absorption)
    return OpticalSystem((Surface(0,0,20,medium_after='glass'),Surface(2,0,20,medium_after='air')),
                         {'air':air,'glass':glass},object_z=-10,image_z=10)


def aperture(radius=2.,object_z=-10.,image_z=10.):
    air=Material('air',n_constant=1.,absorption_constant_per_mm=0.)
    return OpticalSystem((Surface(0,0,radius),),{'air':air},object_z=object_z,image_z=image_z)


class TestRadiometry(unittest.TestCase):
    def test_fresnel_normal_and_brewster_and_tir(self):
        f=fresnel_power(1,1.5,1)
        self.assertAlmostEqual(f['Ts'],.96);self.assertAlmostEqual(f['Tp'],.96)
        b=fresnel_power(1,1.5,np.cos(np.arctan(1.5)))
        self.assertAlmostEqual(b['Rp'],0.,places=14)
        self.assertTrue(fresnel_power(1.5,1,np.cos(np.deg2rad(60)))['tir'])
        self.assertEqual(fresnel_power(1,1,0)['Ts'],1.)

    def test_normal_plate_direct_power_no_etalon(self):
        system=plate();ray=Ray([0,0,-10],[0,0,1]);t=trace_ray(system,ray)
        for mode in ('unpolarized_jones_power','scalar_interface_average'):
            result=trace_transmission(system,ray,t,mode=mode)
            self.assertAlmostEqual(result['transmission'],.96**2,places=12)
            self.assertAlmostEqual(result['ledger_residual'],0.,places=12)

    def test_oblique_jones_not_repeated_unpolarized_average(self):
        system=plate();angle=np.deg2rad(45);ray=Ray([0,0,-10],[np.sin(angle),0,np.cos(angle)])
        trace=trace_ray(system,ray);self.assertTrue(trace.success)
        f=fresnel_power(1,1.5,np.cos(angle));g=fresnel_power(1.5,1,f['cos_transmission'])
        exact=.5*(f['Ts']*g['Ts']+f['Tp']*g['Tp'])
        j=trace_transmission(system,ray,trace);s=trace_transmission(system,ray,trace,mode='scalar_interface_average')
        self.assertAlmostEqual(j['transmission'],exact,places=12)
        self.assertGreater(abs(j['transmission']-s['transmission']),1e-4)

    def test_absorption_uses_oblique_length(self):
        system=plate(absorption=.1);a=np.deg2rad(30);ray=Ray([0,0,-10],[np.sin(a),0,np.cos(a)])
        trace=trace_ray(system,ray);result=trace_transmission(system,ray,trace)
        bare=trace_transmission(plate(),ray,trace)
        ct=fresnel_power(1,1.5,np.cos(a))['cos_transmission']
        self.assertAlmostEqual(result['transmission']/bare['transmission'],np.exp(-.1*2/ct),places=12)
        self.assertAlmostEqual(result['ledger_residual'],0.,places=12)

    def test_missing_absorption_is_lazy_and_explicit(self):
        system=plate(absorption=None);ray=Ray([0,0,-10],[0,0,1]);trace=trace_ray(system,ray)
        self.assertTrue(trace.success)
        result=trace_transmission(system,ray,trace)
        self.assertIn('absorption_missing_assumed_zero:glass',result['assumptions'])
        with self.assertRaises(MissingDataError):trace_transmission(system,ray,trace,missing_absorption='error')
        self.assertEqual(trace_transmission(system,ray,trace,mode='lossless_scalar',missing_absorption='error')['transmission'],1.)

    def test_panel_flux_and_bin_average_cosine(self):
        source=RectangularSource(0,radiance=2.)
        theta=np.linspace(0,np.pi/2,9);phi=np.linspace(-np.pi,np.pi,9)
        result=far_field(None,source,theta,phi,sample_count=16384)
        self.assertAlmostEqual(result['exit_power'],np.pi*source.area_mm2*2*.64,places=10)
        self.assertAlmostEqual(np.sum(result['bin_power']),result['exit_power'],places=10)
        # Mean over each theta band is A L (cos(theta0)+cos(theta1))/2.
        observed=np.sum(result['bin_power'],axis=1)/np.sum(result['solid_angle_sr'],axis=1)
        expected=source.area_mm2*source.panel_radiance*(np.cos(theta[:-1])+np.cos(theta[1:]))/2
        np.testing.assert_allclose(observed,expected,rtol=.004)

    def test_fill_applied_once_and_normalization_cancels(self):
        source=RectangularSource(0,radiance=2.,fill_factor=.5)
        theta=np.linspace(0,np.pi/2,6);phi=np.linspace(-np.pi,np.pi,7)
        a=far_field(None,source,theta,phi,sample_count=2048)
        b=far_field(None,replace(source,fill_factor=1.),theta,phi,sample_count=2048)
        np.testing.assert_allclose(a['intensity']*2,b['intensity'])
        np.testing.assert_allclose(a['relative_far_field_intensity'],b['relative_far_field_intensity'])
        c=replace(source,radiance_basis='panel_average')
        self.assertAlmostEqual(c.hemisphere_power,2*source.hemisphere_power)

    def test_projected_solid_angle_and_limited_cos4(self):
        exact=np.pi*.5**2/(10**2+.5**2)
        axial=projected_solid_angle_aperture([0,0,10],0,.5)
        self.assertAlmostEqual(axial,exact,places=12)
        # Only a small flat aperture subtending a slow cone is a cos^4 fixture.
        tiny=projected_solid_angle_aperture([0,0,10],0,.02)
        off=projected_solid_angle_aperture([4,0,10],0,.02)
        self.assertAlmostEqual(off/tiny,(10/np.sqrt(116))**4,delta=3e-6)

    def test_reverse_actual_ray_ri_raw_integral(self):
        system=aperture()
        result=uniform_scene_ri(system,[[0,0],[2,0]],source_radiance=3.,direction_grid_size=80,
                               direction_cosine_bounds=(-.4,.4,-.25,.25),mode='lossless_scalar',axis_index=0)
        expected=np.array([projected_solid_angle_aperture([x,0,10],0,2)*3 for x in (0,2)])
        np.testing.assert_allclose(result['irradiance'],expected,rtol=.018)
        self.assertAlmostEqual(result['relative_illumination'][0],1.)
        self.assertGreater(result['irradiance'][0],result['irradiance'][1])
        self.assertIn('truncated_image_direction_domain_requires_coverage_check',result['assumptions'])

    def test_solid_angle_jacobians(self):
        self.assertAlmostEqual(np.sum(spherical_bin_solid_angles([0,np.pi/2],[-np.pi,np.pi])),2*np.pi)
        p0,p1,q0,q1=-.4,.8,-.3,.6
        exact=slope_rectangle_solid_angle(p0,p1,q0,q1)
        numerical=dblquad(lambda q,p:1/(1+p*p+q*q)**1.5,p0,p1,lambda p:q0,lambda p:q1)[0]
        self.assertAlmostEqual(exact,numerical,places=11)
        a=.3;b=-.2;eps=1e-5
        rectangle=slope_rectangle_solid_angle(np.tan(a-eps/2),np.tan(a+eps/2),np.tan(b-eps/2),np.tan(b+eps/2))
        self.assertAlmostEqual(rectangle/eps**2,angle_chart_jacobian(a,b),delta=3e-7)

    def test_entrance_sampling_measure_and_coverage_guard(self):
        system=aperture(radius=2.)
        d,w,_=source_direction_samples((0,0),-10,2048,system=system,launch_plane_z_mm=-1,launch_radius_mm=2)
        expected=np.pi*4/(81+4)
        self.assertAlmostEqual(np.sum(w),expected,delta=expected*1e-4)
        with self.assertRaises(InputError):source_direction_samples((3,0),-10,16,system=system,launch_plane_z_mm=-1,launch_radius_mm=2)
        with self.assertRaises(InputError):source_direction_samples((0,0),-10,16,system=system,launch_plane_z_mm=1,launch_radius_mm=2)

    def test_absorbing_output_medium_rejects_infinite_far_field(self):
        system=plate(air_absorption=.01)
        source=RectangularSource(-10)
        with self.assertRaises(ModelApplicabilityError):
            far_field(system,source,[0,np.pi/2],[-np.pi,np.pi],sample_count=8)

    def test_source_packet_no_index_square_energy_gain(self):
        system=OpticalSystem((Surface(0,0,50,medium_after='glass'),),
            {'air':Material('air',n_constant=1.,absorption_constant_per_mm=0.),
             'glass':Material('glass',n_constant=1.5,absorption_constant_per_mm=0.)},object_z=-1,image_z=1)
        source=RectangularSource(-1,width_mm=.1,height_mm=.1,radiance=1.)
        result=far_field(system,source,[0,np.pi/2],[-np.pi,np.pi],sample_count=128)
        self.assertGreater(result['direct_efficiency'],0.)
        self.assertLess(result['direct_efficiency'],1.)
        self.assertAlmostEqual(result['source_hemisphere_power'],result['exit_power']+result['nondirect_power'])


class TestImaging(unittest.TestCase):
    def test_flux_q_kernel_and_crop_no_extra_jacobian(self):
        grid=PixelGrid([0,2,4],[0,3])
        a=SourceKernel(grid,[[.5,.25]],550)
        b=SourceKernel(grid,[[0.,1.]],550)
        result=synthesize_image([a,b],[2,3],[4,5],bitmap_weights=[1,.5])
        np.testing.assert_allclose(result['cell_power'],[[4,9.5]])
        self.assertAlmostEqual(result['collected_power'],15.5)
        self.assertAlmostEqual(result['outside_roi_power'],2.)
        self.assertAlmostEqual(np.sum(result['irradiance']*grid.cell_areas_mm2),13.5)

    def test_reject_lens_metric_and_bad_mass(self):
        grid=PixelGrid([0,1],[0,1])
        with self.assertRaises(InputError):SourceKernel(grid,[[1]],550,provenance='lens_metric')
        with self.assertRaises(InputError):SourceKernel(grid,[[1.1]],550)

    def test_spectral_registration_and_complex_otf(self):
        grid=PixelGrid([0,1,2],[0,1])
        a=synthesize_image([SourceKernel(grid,[[1,0]],450)],[1],[1])
        b=synthesize_image([SourceKernel(grid,[[0,1]],650)],[1],[1])
        np.testing.assert_allclose(registered_spectral_sum([a,b])['cell_power'],[[1,1]])
        result=combine_complex_otfs([[1,1j],[1,-1j]],[1,1])
        np.testing.assert_allclose(result['mtf'],[1,0],atol=1e-15)

    def test_end_to_end_geometric_kernel_independent_reference(self):
        system=aperture(radius=2.)
        grid=PixelGrid(np.linspace(-5,5,21),np.linspace(-5,5,21))
        result=geometric_source_kernel(system,[0,0],grid,source_radiance=2.,fill_factor=.5,
                sample_count=2048,launch_plane_z_mm=-1,launch_radius_mm=2,mode='lossless_scalar')
        # Collected cone reaches R=2 at a distance 10, not the launch plane's cone.
        expected=np.pi*4/104
        self.assertAlmostEqual(result['collected_power_density'],expected,delta=expected*.003)
        self.assertAlmostEqual(result['captured_fraction'],1.,places=12)

    def test_end_to_end_bitmap_scales_and_keeps_roi_loss(self):
        system=aperture(radius=2.)
        grid=PixelGrid(np.linspace(-1,1,9),np.linspace(-1,1,9))
        source=RectangularSource(-10,width_mm=.2,height_mm=.2,radiance=1.)
        a=geometric_image(system,[[1,1]],source,grid,samples_per_source=128,mode='lossless_scalar',launch_plane_z_mm=-1,launch_radius_mm=2)
        b=geometric_image(system,[[2,2]],source,grid,samples_per_source=128,mode='lossless_scalar',launch_plane_z_mm=-1,launch_radius_mm=2)
        np.testing.assert_allclose(b['cell_power'],2*a['cell_power'])
        self.assertGreater(a['outside_roi_power'],0.)
        self.assertEqual(a['source_pixel_metadata'],(320,80))


if __name__=='__main__':unittest.main()
