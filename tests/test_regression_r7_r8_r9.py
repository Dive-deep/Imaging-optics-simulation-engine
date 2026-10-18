"""Expected-behavior regressions for independent software review R7–R9.

The numerical boundary specimen, flux ledger and malformed group inputs come
from audit_round2/software; these import the working engine, never that snapshot.
"""
import unittest

import numpy as np

from sequential_optics import (EngineSettings, RaySettings, SamplingSettings,
    OpticalSystem, Surface, Material, Field, RunRequest, evaluate)
from sequential_optics.errors import InputError
from sequential_optics.geometry import trace_bundle
from sequential_optics.metrics import spot_diagram, common_geometric_refocus
from sequential_optics.imaging import PixelGrid, SourceKernel, synthesize_image
from sequential_optics.thermal import LensGroup, thermal_system


def lens():
    return OpticalSystem(
        (Surface(0., .02, 2., medium_after='glass'), Surface(3., -.02, 2.)),
        {'air': Material('ideal air', n_constant=1., index_temperature_invariant=True),
         'glass': Material('G-only material', wavelengths_nm=(550.,), indices=(1.5,))},
        image_z=52.)


class TestR7SettingsPropagation(unittest.TestCase):
    def test_batch_uses_requested_tolerance_for_spot_and_common_focus(self):
        settings=EngineSettings(ray=RaySettings(direction_tol=1e-15),
            sampling=SamplingSettings(boundary_clearance_mm=1e-12, pupil_rings=2,
                                      pupil_arms=8, include_pupil_boundary=False))
        system=OpticalSystem((Surface(0., 0., 1.),),
            {'air': Material('ideal air', n_constant=1.)}, image_z=1e-3)
        field=Field(x=np.pi/2-1e-13)
        bundle=trace_bundle(system,field,settings=settings)
        direct=spot_diagram(system,bundle,settings=settings)
        focus=common_geometric_refocus(system,(bundle,),settings=settings)
        result=evaluate(system,RunRequest(fields=(field,),metrics=('spot','common_focus'),
            include_ray_history=False,fail_fast=True),settings)
        w=result['states'][0]['wavelengths'][0]
        self.assertEqual(int(bundle.success_mask.sum()),len(bundle.traces))
        self.assertEqual(w['fields'][0]['spot']['status'],'ok')
        self.assertEqual(w['common_focus']['status'],'ok')
        self.assertEqual(w['common_focus']['value']['status'],focus['status'])
        np.testing.assert_allclose(w['fields'][0]['spot']['value']['coordinates_mm'],direct['coordinates_mm'])
        self.assertEqual(w['fields'][0]['spot']['value']['rms_radius'],direct['rms_radius'])
        self.assertEqual(result['settings']['ray']['direction_tol'],settings.ray.direction_tol)


class TestR8ImmutableImageValues(unittest.TestCase):
    def test_inputs_are_copied_and_array_write_flags_cannot_be_unlocked(self):
        edges=np.array([0.,1.]);mass=np.array([[1.]])
        grid=PixelGrid(edges,edges);kernel=SourceKernel(grid,mass,550.)
        edges[1]=-1.;mass[0,0]=2.
        for array in (grid.x_edges_mm,grid.y_edges_mm,kernel.probability_mass):
            self.assertFalse(array.flags.writeable)
            with self.assertRaises(ValueError): array.flat[0]=-1.
            with self.assertRaises(ValueError): array.setflags(write=True)
            # Nor can exposing the ndarray base unlock the immutable buffer.
            if isinstance(array.base,np.ndarray):
                with self.assertRaises(ValueError): array.base.setflags(write=True)
        result=synthesize_image([kernel],[1.],[1.])
        self.assertEqual(result['collected_power'],1.)
        self.assertEqual(result['detected_roi_power'],1.)
        self.assertEqual(result['power_balance_residual'],0.)
        self.assertEqual(result['irradiance'][0,0],1.)

    def test_synthesis_revalidates_forged_kernel_mass(self):
        # frozen dataclasses are not a security boundary against reflection.
        # A deserializer/foreign extension can bypass constructors; reject it.
        for mass in ([[2.]], [[-1.]], [[np.nan]], [[1.,0.]]):
            with self.subTest(mass=mass):
                kernel=SourceKernel(PixelGrid([0,1],[0,1]),[[1.]],550.)
                object.__setattr__(kernel,'probability_mass',np.array(mass))
                with self.assertRaises(InputError): synthesize_image([kernel],[1.],[1.])

    def test_synthesis_revalidates_forged_grid_and_provenance(self):
        grid=PixelGrid([0,1],[0,1]);kernel=SourceKernel(grid,[[1.]],550.)
        object.__setattr__(grid,'x_edges_mm',np.array([0.,-1.]))
        with self.assertRaises(InputError): synthesize_image([kernel],[1.],[1.])
        kernel=SourceKernel(PixelGrid([0,1],[0,1]),[[1.]],550.)
        object.__setattr__(kernel,'provenance','lens_metric')
        with self.assertRaises(InputError): synthesize_image([kernel],[1.],[1.])

    def test_cropped_kernel_preserves_power_and_nonuniform_pixel_measure(self):
        grid=PixelGrid([0.,1.,3.],[0.,2.])
        kernel=SourceKernel(grid,[[.2,.3]],550.)
        result=synthesize_image([kernel],[2.],[3.],bitmap_weights=[.5])
        np.testing.assert_allclose(result['cell_power'],[[.6,.9]])
        np.testing.assert_allclose(result['irradiance'],[[.3,.225]])
        self.assertAlmostEqual(result['detected_roi_power'],1.5)
        self.assertAlmostEqual(result['outside_roi_power'],1.5)
        self.assertAlmostEqual(result['collected_power'],3.)
        self.assertAlmostEqual(result['power_balance_residual'],0.)

    def test_roundoff_mass_excess_is_exposed_in_power_ledger(self):
        kernel=SourceKernel(PixelGrid([0.,1.],[0.,1.]),[[1.+5e-13]],550.)
        result=synthesize_image([kernel],[1.],[1.])
        self.assertGreater(result['power_balance_residual'],0.)
        self.assertLess(result['power_balance_residual'],1e-12)
        self.assertEqual(result['outside_roi_power'],0.)

    def test_nonfinite_power_and_unrepresentable_cell_areas_are_rejected(self):
        kernel=SourceKernel(PixelGrid([0.,1.],[0.,1.]),[[1.]],550.)
        with self.assertRaises(InputError): synthesize_image([kernel],[1e300],[1e300])
        with self.assertRaises(InputError): PixelGrid([0.,1e200],[0.,1e200])
        with self.assertRaises(InputError): PixelGrid([0.,1e-200],[0.,1e-200])
        tiny=SourceKernel(PixelGrid([0.,1e-155],[0.,1e-155]),[[1.]],550.)
        with self.assertRaises(InputError): synthesize_image([tiny],[1.],[1.])


class TestR9ThermalIndices(unittest.TestCase):
    def test_float_bool_and_nonscalar_indices_are_rejected(self):
        for front,back in ((0.,1.),(False,True),(np.bool_(False),np.bool_(True)),
                           ('0','1'),(np.array(0),np.array(1)),(None,1)):
            with self.subTest(front=front,back=back):
                with self.assertRaises(InputError): LensGroup(front,back)

    def test_numpy_integer_scalars_are_normalized_and_work(self):
        group=LensGroup(np.int64(0),np.int32(1))
        self.assertIs(type(group.front_index),int)
        self.assertIs(type(group.back_index),int)
        result=thermal_system(lens(),20.,lens_groups=(group,))
        self.assertEqual(result.temperature_c,20.)
        self.assertEqual(result.surfaces,lens().surfaces)

    def test_bad_dictionary_rejected_when_request_is_constructed(self):
        with self.assertRaises(InputError):
            RunRequest(metrics=('first_order',),temperatures_c=(20.,),fail_fast=False,
                thermal_options={'lens_groups':[{'front_index':0.,'back_index':1.}]})
        with self.assertRaises(InputError): RunRequest(thermal_options={'lens_groups':[object()]})

    def test_runtime_group_revalidation_stays_local_to_thermal_state(self):
        group=LensGroup(0,1)
        request=RunRequest(metrics=('first_order',),temperatures_c=(None,20.),fail_fast=False,
            thermal_options={'lens_groups':[group]})
        object.__setattr__(group,'front_index',0.)
        result=evaluate(lens(),request)
        self.assertEqual(result['states'][0]['status'],'ok')
        self.assertEqual(result['states'][0]['wavelengths'][0]['first_order']['status'],'ok')
        self.assertEqual(result['states'][1]['status'],'failed')
        self.assertEqual(result['states'][1]['error_type'],'InputError')

    def test_optional_nominal_green_needs_no_unused_thermal_or_spectral_data(self):
        settings=EngineSettings(sampling=SamplingSettings(pupil_rings=2,pupil_arms=8))
        result=evaluate(lens(),RunRequest(metrics=('first_order','spot'),wavelengths=('G','B'),
            temperatures_c=(None,20.,40.),include_ray_history=False),settings)
        nominal,reference,hot=result['states']
        for state in (nominal,reference):
            self.assertEqual(state['wavelengths'][0]['first_order']['status'],'ok')
            self.assertEqual(state['wavelengths'][0]['fields'][0]['spot']['status'],'ok')
            self.assertEqual(state['wavelengths'][1]['first_order']['status'],'missing_data')
        self.assertEqual(hot['status'],'missing_data')


if __name__=='__main__': unittest.main()
