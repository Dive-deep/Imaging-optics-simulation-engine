"""Independent laws and adversarial prescriptions, not output snapshots."""
import unittest
from dataclasses import replace
import numpy as np
from sequential_optics.models import Surface, OpticalSystem, Ray, Field
from sequential_optics.settings import EngineSettings
from sequential_optics.errors import GeometryError, NumericalError, TotalInternalReflection, InputError
from sequential_optics.geometry import sag, gradient, hessian, intersect, refract, trace_ray, aim_ray, trace_bundle, validate_geometry
from sequential_optics.paraxial import paraxial_properties, resolve_aperture


class ConstantIndex:
    def __init__(self,n): self.n=n
    def index(self,wavelength_nm,temperature_c=None): return self.n


def system(surfaces,**kwargs):
    return OpticalSystem(tuple(surfaces),{'air':ConstantIndex(1.),'glass':ConstantIndex(1.5)},**kwargs)


class SagTests(unittest.TestCase):
    def test_zero_curvature_a2_is_not_plane(self):
        s=Surface(0,0,2,coefficients=(.1,),surface_type='even_asphere')
        self.assertAlmostEqual(float(sag(s,1,1)),.2)
        np.testing.assert_allclose(gradient(s,1,1),[.2,.2])
        np.testing.assert_allclose(hessian(s,0,0),np.eye(2)*.2)
        p=paraxial_properties(system([replace(s,medium_after='glass')]))
        self.assertAlmostEqual(p['power_abs_per_mm'],.1)

    def test_gradient_hessian_independent_differences(self):
        s=Surface(0,.03,4,conic=-.7,coefficients=(.001,-2e-5,1e-7,-1e-10,2e-14,-1e-18),surface_type='even_asphere')
        x=np.array([1.1,-.8]);h=1e-4
        g=np.array([(float(sag(s,*(x+np.eye(2)[i]*h)))-float(sag(s,*(x-np.eye(2)[i]*h))))/(2*h) for i in range(2)])
        np.testing.assert_allclose(gradient(s,*x),g,rtol=1e-8,atol=1e-10)
        hh=np.column_stack([(gradient(s,*(x+np.eye(2)[i]*h))-gradient(s,*(x-np.eye(2)[i]*h)))/(2*h) for i in range(2)])
        np.testing.assert_allclose(hessian(s,*x),hh,rtol=1e-8,atol=1e-10)

    def test_domain_fault_and_vectorized_shapes(self):
        s=Surface(0,1,2)
        with self.assertRaises(GeometryError): sag(s,2,0)
        self.assertEqual(hessian(s,np.array([0,.1]),0).shape,(2,2,2))

    def test_interior_overlap_not_axis_or_rim(self):
        a=Surface(0,0,1,medium_after='glass')
        b=Surface(1,0,1,coefficients=(-8,8),surface_type='even_asphere')
        self.assertGreater(1+float(sag(b,0,0)),0)
        self.assertGreater(1+float(sag(b,1,0)),0)
        with self.assertRaisesRegex(GeometryError,'surface_overlap'): validate_geometry(system([a,b]))


class IntersectionTests(unittest.TestCase):
    def test_two_positive_roots_first_hit(self):
        s=Surface(0,0,4,coefficients=(1,),surface_type='even_asphere')
        hit,t=intersect(s,Ray([0,0,-2],[1,0,3]))
        np.testing.assert_allclose(hit,[1,0,1],atol=1e-9)
        self.assertAlmostEqual(t,np.sqrt(10),places=8)

    def test_tangent_not_silently_skipped(self):
        s=Surface(0,0,3,coefficients=(1,),surface_type='even_asphere')
        with self.assertRaisesRegex(GeometryError,'grazing_degenerate'):
            intersect(s,Ray([-2,0,0],[1,0,0]))

    def test_first_opaque_hit_not_later_active_hit(self):
        s=Surface(0,0,.8,mechanical_radius=4,coefficients=(1,),surface_type='even_asphere')
        tr=trace_ray(system([s]),Ray([0,0,-2],[1,0,3]))
        self.assertEqual(tr.status,'aperture_blocked')
        self.assertFalse(tr.success)

    def test_no_hit_and_sphere_original_branch(self):
        with self.assertRaisesRegex(GeometryError,'escaped'):
            intersect(Surface(0,0,1),Ray([2,0,-1],[0,0,1]))
        s=Surface(0,.1,4)
        hit,t=intersect(s,Ray([1,0,-1],[0,0,1]))
        self.assertAlmostEqual(hit[2],10-np.sqrt(99),places=9)

    def test_zero_thickness_dummy_surface(self):
        tr=trace_ray(system([Surface(0,0,3),Surface(0,0,3)]),Ray([.2,0,-1],[0,0,1]))
        self.assertTrue(tr.success,tr.message)
        np.testing.assert_allclose(tr.opl,[1,1],atol=1e-12)


class SnellTraceTests(unittest.TestCase):
    def test_vector_snell_invariants_and_orientation(self):
        theta=.5;d=np.array([np.sin(theta),0,np.cos(theta)])
        a=refract(d,[0,0,-1],1,1.5)
        np.testing.assert_allclose(a,refract(d,[0,0,1],1,1.5))
        self.assertAlmostEqual(np.linalg.norm(a),1)
        self.assertAlmostEqual(1.5*a[0],np.sin(theta))
        self.assertGreater(a[2],0)
        with self.assertRaises(TotalInternalReflection): refract([.9,0,np.sqrt(.19)],[0,0,1],1.5,1)

    def test_vector_snell_normalization_accepts_large_finite_scales(self):
        expected=refract([1,0,1],[0,0,-1],1.,1.5)
        actual=refract([1e308,0,1e308],[0,0,-1e308],1.,1.5)
        self.assertTrue(np.all(np.isfinite(actual)))
        np.testing.assert_allclose(actual,expected,atol=1e-15)

    def test_plate_restores_angle_and_correct_lateral_shift_opl(self):
        theta=.3;thickness=4;theta_t=np.arcsin(np.sin(theta)/1.5)
        sys=system([Surface(0,0,10,medium_after='glass'),Surface(thickness,0,10)])
        ray=Ray([-np.tan(theta),0,-1],[np.sin(theta),0,np.cos(theta)])
        tr=trace_ray(sys,ray)
        self.assertTrue(tr.success,tr.message)
        np.testing.assert_allclose(tr.last_direction,ray.direction,atol=1e-12)
        self.assertAlmostEqual(tr.last_position[0],thickness*np.tan(theta_t),places=9)
        self.assertAlmostEqual(tr.last_opl,1/np.cos(theta)+1.5*thickness/np.cos(theta_t),places=9)
        self.assertEqual(tr.last_power,1.)

    def test_aperture_tolerance(self):
        sys=system([Surface(0,0,1,mechanical_radius=2)])
        self.assertTrue(trace_ray(sys,Ray([1,0,-1],[0,0,1])).success)
        self.assertEqual(trace_ray(sys,Ray([1.001,0,-1],[0,0,1])).status,'aperture_blocked')


class AimFirstOrderTests(unittest.TestCase):
    def make_lens(self,**kwargs):
        return system([Surface(0,.02,6,medium_after='glass'),Surface(3,-.02,6)],stop_index=1,stop_radius=3,**kwargs)

    def test_stop_actual_sag_finite_and_infinite(self):
        for field,sys in [(Field(.03,-.01),self.make_lens()),(Field(.5,.1,'height'),self.make_lens(object_z=-30))]:
            for target in [[0,0],[.3,-.4]]:
                ray=aim_ray(sys,field,target)
                tr=trace_ray(sys,ray,stop_at=sys.stop_index)
                self.assertTrue(tr.success,tr.message)
                np.testing.assert_allclose(tr.last_position[:2],np.asarray(target)*3,atol=2e-9)
                self.assertAlmostEqual(tr.last_position[2],3+float(sag(sys.surfaces[1],*tr.last_position[:2])),places=9)

    def test_stop_inside_glass_dummy(self):
        sys=system([Surface(0,.02,6,medium_after='glass'),Surface(1.5,0,6,medium_after='glass'),Surface(3,-.02,6)],stop_index=1,stop_radius=2)
        ray=aim_ray(sys,Field(.04,0),[.6,0]);tr=trace_ray(sys,ray)
        self.assertTrue(tr.success,tr.message)
        self.assertAlmostEqual(tr.positions[1,0],1.2,places=8)

    def test_plane_wave_launch_phase_and_bundle_weights(self):
        sys=system([Surface(0,0,10)],stop_radius=3)
        field=Field(.04,.01)
        ray=aim_ray(sys,field,[.5,.2]);tr=trace_ray(sys,ray)
        self.assertAlmostEqual(tr.last_opl,float(ray.direction@tr.last_position),places=11)
        bundle=trace_bundle(sys,field,rings=4,arms=12)
        self.assertTrue(bundle.success_mask.all())
        self.assertAlmostEqual(bundle.weights.sum(),1.)
        self.assertAlmostEqual(np.sum(bundle.weights*np.sum(bundle.pupil_xy**2,axis=1)),.5)
        self.assertEqual(len(bundle.traces),4*12+12+1)
        self.assertEqual(np.sum(bundle.weights==0),13)
        self.assertEqual(np.sum(np.isclose(np.linalg.norm(bundle.pupil_xy,axis=1),1.)),12)
        interior=trace_bundle(sys,field,rings=4,arms=12,include_boundary=False)
        self.assertEqual(len(interior.traces),48)
        self.assertAlmostEqual(interior.weights.sum(),1.)

    def test_first_order_real_ray_limit_and_epd(self):
        sys=self.make_lens();p=paraxial_properties(sys)
        self.assertAlmostEqual(p['determinant'],1.,places=12)
        tr=trace_ray(sys,Ray([1e-4,0,-1],[0,0,1]))
        z=tr.last_position[2]-tr.last_position[0]*tr.last_direction[2]/tr.last_direction[0]
        self.assertAlmostEqual(z,p['paraxial_image_z'],places=7)
        resolved=resolve_aperture(sys,epd_mm=2)
        self.assertGreater(resolved.stop_semi_diameter,0)
        ray=Ray([1,0,-1],[0,0,1]);tr=trace_ray(resolved,ray,ignore_apertures=True)
        self.assertAlmostEqual(resolved.stop_semi_diameter,tr.positions[1,0],places=3)
        with self.assertRaises(InputError): resolve_aperture(sys,epd_mm=2,f_number=4)

    def test_primary_wavelength_is_configured_and_no_thermal_required(self):
        settings=EngineSettings(primary_wavelength='green',wavelengths_nm={'green':532})
        b=trace_bundle(self.make_lens(),Field(),settings=settings,rings=2,arms=4)
        self.assertEqual(b.wavelength_nm,532)
        self.assertTrue(b.success_mask.all())

    def test_stop_TIR_is_physical_not_aim_solver_failure(self):
        sys=system([Surface(0,0,3)],ambient_material='glass')
        ray=aim_ray(sys,Field(1.,0),[.2,0])
        trace=trace_ray(sys,ray)
        self.assertEqual(trace.status,'tir')

    def test_finite_object_angle_is_not_silently_infinite_illumination(self):
        with self.assertRaisesRegex(InputError,'height field'):
            aim_ray(self.make_lens(object_z=-30),Field(.02,0),[0,0])

    def test_finite_collimator_exact_STOP_rim_not_false_clipped_by_aim_tolerance(self):
        # Independent review exposed five false STOP clips with default
        # aim_tol=1e-9 > aperture/intersection_tol=1e-10 on this finite field.
        surfaces=(Surface(0,1/50,6,medium_after='crown'),Surface(3,-1/50,6),
                  Surface(5,-1/40,5,medium_after='flint'),Surface(7,1/40,5),
                  Surface(9,1/45,6,medium_after='crown'),Surface(12,-1/45,6))
        materials={'air':ConstantIndex(1.),'crown':ConstantIndex(1.5),'flint':ConstantIndex(1.62)}
        sys=OpticalSystem(surfaces,materials,object_z=-78.52124279501913,stop_index=2,stop_radius=4.)
        bundle=trace_bundle(sys,Field(3.2,0,'height'),rings=10,arms=36)
        rim=np.isclose(np.linalg.norm(bundle.pupil_xy,axis=1),1.)
        self.assertEqual(rim.sum(),36)
        for trace in np.asarray(bundle.traces,dtype=object)[rim]:
            self.assertTrue(trace.success,trace.message)
            self.assertLessEqual(np.hypot(*trace.positions[2,:2]),4+EngineSettings().ray.intersection_tol_mm)

    def test_report_reference_explicit_unused_material_independence(self):
        sys=self.make_lens()
        original=paraxial_properties(sys)
        modified=replace(sys,materials={**sys.materials,'unused':ConstantIndex(1.0003)})
        same=paraxial_properties(modified)
        self.assertEqual(original['efl_report_mm'],same['efl_report_mm'])
        chosen=paraxial_properties(modified,settings=EngineSettings(report_reference_material='unused'))
        self.assertAlmostEqual(chosen['efl_report_mm']/same['efl_report_mm'],1.0003)
        self.assertEqual(chosen['efl_image_mm'],same['efl_image_mm'])


if __name__=='__main__': unittest.main()
