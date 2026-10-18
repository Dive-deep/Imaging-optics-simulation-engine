"""Independent camera-design review: physically derived plane-interface foci.

These tests use public computational APIs. There is no Zemax golden file.
"""
import unittest
import numpy as np
from sequential_optics.models import OpticalSystem,Surface,Field
from sequential_optics.materials import Material
from sequential_optics.metrics import differential_field_curvature,ray_fans
from sequential_optics.api import evaluate,RunRequest


def plane_system():
    return OpticalSystem((Surface(0,0,20,medium_after='glass'),),
        {'air':Material('air',n_constant=1.),'glass':Material('glass',n_constant=1.5)},
        object_z=-10,image_z=5)


class IndependentCameraReviewTests(unittest.TestCase):
    def test_coddington_plane_refraction_foci_and_signs(self):
        # Exact derivative of vector Snell at a plane: oblique finite source.
        # ci=a/sqrt(a²+h²), ct=sqrt(1-eta²(1-ci²));
        # delta_T=-a ct³/(eta ci³); delta_S=-a ct/(eta ci).
        # Both are virtual and split at finite incidence, even for plane sag.
        a=10.;h=5.;eta=1/1.5;ci=a/np.hypot(a,h);ct=np.sqrt(1-eta**2*(1-ci**2))
        expected=np.array([-a*ct**3/(eta*ci**3),-a*ct/(eta*ci)])
        result=differential_field_curvature(plane_system(),Field(h,0,'height'),reference_z=0)
        self.assertEqual(result['status'],'ok')
        np.testing.assert_allclose(result['focus_z_TS_mm'],expected,rtol=2e-8,atol=1e-7)
        self.assertTrue(np.all(result['focus_z_TS_mm']<0))
        self.assertLess(result['astigmatic_focus_separation_mm'],0)

    def test_plane_refraction_foci_rotate_covariantly(self):
        h=5.;phi=.71
        axis=differential_field_curvature(plane_system(),Field(h,0,'height'),reference_z=0)
        rotated=differential_field_curvature(plane_system(),Field(h*np.cos(phi),h*np.sin(phi),'height'),reference_z=0)
        np.testing.assert_allclose(axis['focus_z_TS_mm'],rotated['focus_z_TS_mm'],rtol=3e-8)

    def test_finite_angular_fans_use_three_dimensional_meridians(self):
        field=Field(5,7,'height');system=plane_system();samples=np.array([-.3,0,.3])
        result=ray_fans(system,field,angular=True,pupil_coordinates=samples)
        radial=np.array([5.,7.]);t=radial/np.linalg.norm(radial);s=np.array([t[1],-t[0]])
        def outgoing(target):
            incident=np.array([target[0]-5,target[1]-7,10.]);incident/=np.linalg.norm(incident)
            xy=incident[:2]/1.5
            return np.array([*xy,np.sqrt(1-xy@xy)])
        dc=outgoing([0,0]);T=np.array([*t,0.]);T-=np.dot(T,dc)*dc;T/=np.linalg.norm(T)
        S=np.cross(T,dc)
        for name,pupil_axis in [('T',t),('S',s)]:
            expected=[]
            for q in samples:
                d=outgoing(q*system.stop_semi_diameter*pupil_axis)
                expected.append([np.arctan2(d@T,d@dc),np.arctan2(d@S,d@dc)])
            np.testing.assert_allclose(result[name+'_fan_TS_rad'],expected,atol=1e-12)
        np.testing.assert_allclose(result['T_fan_TS_rad'][:,1],0,atol=1e-12)

    def test_batch_api_focal_and_angular_ray_fans_supported(self):
        for output in ('focal','afocal'):
            result=evaluate(plane_system(),RunRequest(fields=(Field(1,2,'height'),),
                    metrics=('ray_fans',),output_mode=output,fail_fast=True))
            self.assertEqual(result['states'][0]['wavelengths'][0]['fields'][0]['ray_fans']['status'],'ok')


if __name__=='__main__':unittest.main()
