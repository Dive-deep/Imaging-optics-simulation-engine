"""Independent review C3/C4: physical spectral measure and portable inputs."""
from dataclasses import replace
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from sequential_optics import EngineSettings, RunRequest, evaluate
from sequential_optics.errors import InputError
from sequential_optics.imaging import PixelGrid, SourceKernel, synthesize_image, registered_spectral_sum
from sequential_optics.io import load_system, save_system, system_to_dict
from sequential_optics.materials import Material
from sequential_optics.models import OpticalSystem, Surface
from sequential_optics.settings import RaySettings, SamplingSettings


def spectral_inputs():
    # Unequal physical cell areas expose accidental broadcasting or wrong measure.
    grid=PixelGrid([0.,1.,3.],[0.,2.])
    first=synthesize_image([SourceKernel(grid,[[.75,.25]],450.)],[2.],[3.])
    second=synthesize_image([SourceKernel(grid,[[.1,.4]],650.)],[5.],[2.])
    return grid,first,second


class SpectralCompositionBoundaryTests(unittest.TestCase):
    def test_actual_images_preserve_registered_cell_power_and_area_measure(self):
        grid,a,b=spectral_inputs()
        result=registered_spectral_sum([a,b])
        np.testing.assert_allclose(result['cell_power'],[[5.5,5.5]])
        np.testing.assert_allclose(result['irradiance'],[[2.75,1.375]])
        self.assertAlmostEqual(float(np.sum(result['irradiance']*grid.cell_areas_mm2)),11.)
        self.assertEqual(result['wavelengths_nm'],[450.,650.])
        # Cropped tails stay lost; spectral composition must not normalize them.
        self.assertGreater(a['collected_power']+b['collected_power'],11.)

    def test_response_maps_wavelengths_to_channels_once(self):
        grid,a,b=spectral_inputs()
        result=registered_spectral_sum([a,b],response=[[1.,0.],[0.,2.],[.25,.75]])
        expected=np.array([[[4.5,1.5]],[[2.,8.]],[[1.875,3.375]]])
        np.testing.assert_allclose(result['cell_power'],expected)
        np.testing.assert_allclose(result['irradiance']*grid.cell_areas_mm2,expected)

    def test_zero_power_and_zero_response_are_valid(self):
        grid,a,b=spectral_inputs()
        zero={**a,'cell_power':np.zeros(grid.shape)}
        np.testing.assert_array_equal(registered_spectral_sum([zero])['cell_power'],np.zeros(grid.shape))
        np.testing.assert_array_equal(registered_spectral_sum([a,b],[[0.,0.]])['cell_power'],np.zeros((1,*grid.shape)))

    def test_wrong_shape_cannot_broadcast_into_valid_irradiance(self):
        grid=PixelGrid([0.,1.,2.],[0.,1.,2.])
        for power in ([[1.]], [1.,1.], [[[1.,1.],[1.,1.]]]):
            with self.subTest(shape=np.shape(power)):
                with self.assertRaises(InputError):
                    registered_spectral_sum([dict(grid=grid,cell_power=power,wavelength_nm=550.)])

    def test_negative_nonfinite_and_complex_power_are_rejected(self):
        grid,a,_=spectral_inputs()
        for power in ([[1.,-1.]], [[1.,np.nan]], [[1.,np.inf]], [[1.+2.j,0.]], [[1.],[1.,2.]]):
            with self.subTest(power=power):
                with self.assertRaises(InputError): registered_spectral_sum([{**a,'cell_power':power}])

    def test_invalid_wavelengths_are_rejected(self):
        _,a,_=spectral_inputs()
        for wavelength in (0.,-1.,np.nan,np.inf,'unknown',None,
                           np.complex64(550.+20.j),np.complex128(550.+20.j)):
            with self.subTest(wavelength=wavelength):
                with self.assertRaises(InputError): registered_spectral_sum([{**a,'wavelength_nm':wavelength}])

    def test_missing_required_payload_and_invalid_grid_are_input_errors(self):
        _,a,_=spectral_inputs()
        for bad in ({}, {**a,'grid':None}, None):
            with self.subTest(payload=type(bad).__name__):
                with self.assertRaises(InputError): registered_spectral_sum([bad])

    def test_different_physical_grid_is_rejected_even_if_shape_matches(self):
        _,a,b=spectral_inputs()
        b={**b,'grid':PixelGrid([1.,2.,4.],[0.,2.])}
        with self.assertRaises(InputError): registered_spectral_sum([a,b])

    def test_invalid_response_is_rejected(self):
        _,a,b=spectral_inputs()
        for response in ([1.,1.], [[1.]], [[1.,-1.]], [[1.,np.nan]],
                         [[1.,np.inf]], [[1.+1.j,0.]], [[1.],[1.,2.]], np.empty((0,2))):
            with self.subTest(response=response):
                with self.assertRaises(InputError): registered_spectral_sum([a,b],response)

    def test_sum_overflow_from_finite_engine_generated_images_is_rejected(self):
        grid=PixelGrid([0.,1.],[0.,1.])
        results=[synthesize_image([SourceKernel(grid,[[1.]],w)],[6e307],[1.])
                 for w in (450.,550.,650.)]
        with self.assertRaises(InputError): registered_spectral_sum(results)

    def test_response_product_overflow_is_rejected(self):
        grid=PixelGrid([0.,1.],[0.,1.])
        result=synthesize_image([SourceKernel(grid,[[1.]],550.)],[1e308],[1.])
        with self.assertRaises(InputError): registered_spectral_sum([result],[[2.]])

    def test_total_channel_power_overflow_is_rejected(self):
        grid=PixelGrid([0.,1.,2.],[0.,1.])
        results=[synthesize_image([SourceKernel(grid,k,w)],[1e308],[1.])
                 for k,w in (([[1.,0.]],450.),([[0.,1.]],650.))]
        with self.assertRaises(InputError): registered_spectral_sum(results)

    def test_irradiance_overflow_is_rejected_even_with_finite_power(self):
        grid=PixelGrid([0.,1e-200],[0.,1e-100])
        with self.assertRaises(InputError):
            registered_spectral_sum([dict(grid=grid,cell_power=[[1e10]],wavelength_nm=550.)])


class PortableNumpyInputTests(unittest.TestCase):
    def specimen(self):
        return OpticalSystem((Surface(np.float32(0.),np.float32(0.),np.int64(1)),),
            {'air':Material('explicit ideal air',n_constant=1.)},
            image_z=np.float32(1.),stop_index=np.int64(0))

    def test_all_accepted_numpy_setting_counts_are_normalized_and_roundtrip(self):
        ray_keys=('root_samples','max_iterations','aim_max_iterations')
        sample_keys=('pupil_rings','pupil_arms','wave_grid_size','fft_padding',
                     'frequency_samples','angular_frequency_samples')
        ray=RaySettings(**{key:np.int64(getattr(RaySettings(),key)) for key in ray_keys})
        sampling=SamplingSettings(**{key:np.int32(getattr(SamplingSettings(),key)) for key in sample_keys})
        settings=EngineSettings(ray=ray,sampling=sampling,wavelengths_nm={'G':np.int64(550)})
        for key in ray_keys: self.assertIs(type(getattr(settings.ray,key)),int)
        for key in sample_keys: self.assertIs(type(getattr(settings.sampling,key)),int)
        # Replacing a frozen configuration must not mutate the caller's object.
        self.assertIs(type(ray.root_samples),np.int64)
        self.assertIs(type(sampling.pupil_rings),np.int32)
        with TemporaryDirectory() as d:
            path=Path(d)/'settings.json'
            settings.to_json(path)
            payload=json.loads(path.read_text())
            self.assertIs(type(payload['sampling']['pupil_rings']),int)
            self.assertEqual(EngineSettings.from_json(path),settings)

    def test_numpy_stop_and_other_numeric_prescription_values_roundtrip(self):
        system=self.specimen()
        self.assertIs(type(system.stop_index),int)
        payload=system_to_dict(system)
        json.dumps(payload,allow_nan=False)  # public to_dict must already be portable
        with TemporaryDirectory() as d:
            path=Path(d)/'system.json'
            save_system(system,path)
            restored=load_system(path)
        self.assertEqual(restored.surfaces,system.surfaces)
        self.assertEqual(restored.stop_index,system.stop_index)
        request=RunRequest(metrics=('first_order',))
        original=evaluate(system,request)
        roundtrip=evaluate(restored,request)
        self.assertEqual(original['states'][0]['wavelengths'][0]['first_order']['status'],'ok')
        self.assertEqual(original['input_sha256'],roundtrip['input_sha256'])

    def test_bool_and_float_counts_and_stop_remain_rejected(self):
        for invalid in (True,np.bool_(True),1.,np.float64(1.)):
            with self.subTest(value=invalid,type=type(invalid).__name__):
                with self.assertRaises(InputError): EngineSettings(ray=RaySettings(max_iterations=invalid))
                with self.assertRaises(InputError): EngineSettings(sampling=SamplingSettings(pupil_rings=invalid))
                with self.assertRaises(InputError): replace(self.specimen(),stop_index=invalid)


if __name__=='__main__': unittest.main()
