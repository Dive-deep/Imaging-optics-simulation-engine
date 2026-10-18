import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
import numpy as np
from sequential_optics.settings import EngineSettings, SamplingSettings
from sequential_optics.models import Surface, Field, Ray
from sequential_optics.errors import InputError
from sequential_optics.io import json_value

class ConfigurationTests(unittest.TestCase):
    def test_default_and_user_green(self):
        s=EngineSettings();self.assertEqual(s.wavelength(),550.)
        edited=EngineSettings(wavelengths_nm={'G':532.},sampling=SamplingSettings(wave_grid_size=32))
        self.assertEqual(edited.wavelength(),532.)
        with self.assertRaises(InputError): edited.wavelength('B')
        with TemporaryDirectory() as d:
            p=Path(d)/'settings.json';edited.to_json(p)
            self.assertEqual(EngineSettings.from_json(p),edited)

    def test_invalid_physics_inputs(self):
        for kwargs in ({'curvature':float('nan')},{'clear_radius':-1},{'surface_type':'fake'},
                       {'coefficients':(1,)},{'mechanical_radius':.1}):
            args=dict(z=0,curvature=0,clear_radius=1);args.update(kwargs)
            with self.assertRaises(InputError):Surface(**args)
        with self.assertRaises(InputError):Field(x=np.pi/2)
        with self.assertRaises(InputError):Ray(np.zeros(3),np.zeros(3))

    def test_complex_otf_not_discarded(self):
        encoded=json_value(np.array([1+0j,-.4+.3j]))
        self.assertEqual(encoded['imag'],[0.,.3])
        self.assertEqual(encoded['real'],[1.,-.4])

    def test_sequential_rows_radius_thickness_and_diameter(self):
        from sequential_optics.io import from_sequential_rows
        from sequential_optics.materials import Material
        m={'air':Material('ideal',n_constant=1.),'glass':Material('synthetic',n_constant=1.5)}
        rows=[{'radius_mm':50.,'thickness_mm':3.,'active_diameter_mm':10.,'medium_after':'glass'},
              {'radius_mm':-50.,'thickness_mm':48.,'active_diameter_mm':10.,'medium_after':'air'}]
        s=from_sequential_rows(rows,m)
        self.assertEqual(s.surfaces[1].z,3.)
        self.assertEqual(s.surfaces[0].clear_radius,5.)
        self.assertEqual(s.image_z,51.)
        with self.assertRaises(InputError):from_sequential_rows(rows,m,image_z=51.)

    def test_inactive_height_preserved_and_scaled_as_metadata(self):
        from sequential_optics.thermal import scale_surface
        s=Surface(0,0,1,mechanical_radius=2,inactive_height_mm=.3)
        self.assertAlmostEqual(scale_surface(s,2).inactive_height_mm,.6)

if __name__=='__main__':unittest.main()
