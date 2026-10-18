"""Independent physics review regressions: reciprocity, etendue, and crop power."""
import unittest
from dataclasses import replace
import numpy as np
from sequential_optics.models import Surface,OpticalSystem,Ray
from sequential_optics.materials import Material
from sequential_optics.geometry import trace_ray
from sequential_optics.radiometry import (uniform_scene_ri,trace_transmission,_reverse_system,
                                         far_field,RectangularSource,fresnel_power)
from sequential_optics.imaging import PixelGrid,SourceKernel,synthesize_image
from sequential_optics.errors import ModelApplicabilityError


class IndependentRadiometryPhysicsReview(unittest.TestCase):
    def test_reverse_ri_high_index_escape_cone_and_index_factor(self):
        # An infinite planar lossless interface illuminated by uniform source
        # radiance has an image-side escape cone sin(theta)<=ns/ni.
        # Integrating ni^2/ns^2 * dl dm over that cone gives pi exactly.
        ns,ni=1.,1.5
        system=OpticalSystem((Surface(0.,0.,50.,medium_after='image'),),
             {'air':Material('object',n_constant=ns),'image':Material('image',n_constant=ni)},
             object_z=-1.,image_z=1.)
        result=uniform_scene_ri(system,[[0.,0.]],source_radiance=1.,mode='lossless_scalar',
                  direction_grid_size=80,direction_cosine_bounds=(-.7,.7,-.7,.7))
        self.assertAlmostEqual(result['irradiance'][0],np.pi,delta=.012*np.pi)
        self.assertGreater(result['trace_status_counts'].get('tir',0),0)

    def test_nonmeridional_jones_direct_path_reciprocity(self):
        air=Material('air',n_constant=1.,absorption_constant_per_mm=0.)
        glass=Material('glass',n_constant=1.6,absorption_constant_per_mm=0.)
        system=OpticalSystem((Surface(0.,1/15.,5.,medium_after='glass'),
                              Surface(2.,-1/22.,5.,medium_after='air')),
                             {'air':air,'glass':glass},object_z=-3.,image_z=10.)
        ray=Ray([.8,.3,-3.],[.04,.07,1.])
        forward=trace_ray(system,ray)
        self.assertTrue(forward.success,forward.message)
        a=trace_transmission(system,ray,forward)
        reverse=_reverse_system(system)
        reverse_origin=forward.last_position+forward.last_direction
        reverse_origin=reverse_origin*np.array([1.,1.,-1.])
        reverse_direction=-forward.last_direction*np.array([1.,1.,-1.])
        rray=Ray(reverse_origin,reverse_direction)
        backward=trace_ray(reverse,rray)
        self.assertTrue(backward.success,backward.message)
        b=trace_transmission(reverse,rray,backward)
        self.assertTrue(np.allclose(backward.positions*np.array([1.,1.,-1.]),forward.positions[::-1],atol=1e-8))
        self.assertAlmostEqual(a['transmission'],b['transmission'],places=10)
        self.assertAlmostEqual(a['ledger_residual'],0.,places=13)
        self.assertAlmostEqual(b['ledger_residual'],0.,places=13)

    def test_crop_and_nonuniform_pixel_measure_preserve_raw_power(self):
        grid=PixelGrid([-2.,-.5,0.,2.],[-1.,0.,3.])
        mass=np.array([[.1,.05,.2],[.04,.01,.1]])
        kernel=SourceKernel(grid,mass,550.)
        result=synthesize_image([kernel],[3.],[2.],bitmap_weights=[.5])
        self.assertAlmostEqual(result['collected_power'],3.)
        self.assertAlmostEqual(result['detected_roi_power'],3.*mass.sum())
        self.assertAlmostEqual(result['outside_roi_power'],3.*(1-mass.sum()))
        self.assertAlmostEqual(np.sum(result['irradiance']*grid.cell_areas_mm2),result['detected_roi_power'])
        # A display crop cannot rescale its mass to unity.
        self.assertLess(np.sum(result['cell_power']),result['collected_power'])

    def test_grazing_equal_media_are_not_tir(self):
        f=fresnel_power(1.5,1.5,0.)
        self.assertFalse(f['tir']);self.assertEqual(f['Ts'],1.);self.assertEqual(f['Tp'],1.)

    def test_absorbing_infinite_output_needs_explicit_model(self):
        system=OpticalSystem((Surface(0.,0.,10.,medium_after='absorbing'),),
            {'air':Material('air',n_constant=1.,absorption_constant_per_mm=0.),
             'absorbing':Material('out',n_constant=1.,absorption_constant_per_mm=.1)},object_z=-1.)
        source=RectangularSource(-1.,width_mm=.1,height_mm=.1,radiance=1.)
        with self.assertRaises(ModelApplicabilityError):
            far_field(system,source,[0.,np.pi/2],[-np.pi,np.pi],sample_count=8)
        # Explicit loss-free ray-distribution metric remains available.
        value=far_field(system,source,[0.,np.pi/2],[-np.pi,np.pi],sample_count=8,mode='lossless_scalar')
        self.assertGreater(value['exit_power'],0.)

if __name__=='__main__': unittest.main()
