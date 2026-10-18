"""Independent area controls for detector overlap; no lens geometry is involved."""
import unittest
import numpy as np

from sequential_optics.errors import InputError, ModelApplicabilityError
from sequential_optics.imaging import PixelGrid
from sequential_optics.materials import Material
from sequential_optics.models import OpticalSystem, Surface
from sequential_optics.settings import EngineSettings, SamplingSettings
from sequential_optics.source_wave import _polygon_area, source_diffraction_kernel, wave_to_source_kernel


def unit_cell(origin, dx=.01, dy=None, basis=None):
    """A convex rectangle with exactly one unit of probability mass."""
    dy = dx if dy is None else dy
    return {'wavelength_nm': 550., 'wave': {
        'mode': 'focal',
        'illumination': 'physical_source_triangle_boundary_amplitude',
        'output_basis': np.eye(3) if basis is None else basis,
        'image_reference_point_mm': np.array([*origin, 20.]),
        'dx_mm': dx, 'dy_mm': dy,
        'x_mm': np.array([0.]), 'y_mm': np.array([0.]),
        'psf': np.array([[1. / (dx * dy)]]),
    }}


def independent_rectangle_overlap(x_edges, y_edges, dx, dy):
    """Exact separable interval-overlap formula, independent of polygon code."""
    x = np.maximum(0., np.minimum(x_edges[1:], dx / 2) -
                   np.maximum(x_edges[:-1], -dx / 2)) / dx
    y = np.maximum(0., np.minimum(y_edges[1:], dy / 2) -
                   np.maximum(y_edges[:-1], -dy / 2)) / dy
    return y[:, None] * x[None, :]


class DetectorOverlapTranslationTests(unittest.TestCase):
    def test_convex_unit_cell_full_capture_is_translation_invariant(self):
        for shift, dx in ((0., .01), (5., .01), (1000., .0001)):
            with self.subTest(shift=shift, dx=dx):
                edges = shift + np.linspace(-dx, dx, 21)
                self.assertTrue(np.all(np.diff(edges) > 0))
                kernel, _ = wave_to_source_kernel(unit_cell((shift, shift), dx),
                                                  PixelGrid(edges, edges))
                self.assertAlmostEqual(kernel.captured_fraction, 1., delta=2e-14)
                expected = independent_rectangle_overlap(edges-shift, edges-shift, dx, dx)
                np.testing.assert_allclose(kernel.probability_mass, expected,
                                           atol=2e-15, rtol=1e-12)

    def test_partial_crop_matches_interval_integral_without_renormalizing(self):
        dx, dy = .01, .003
        for ox, oy in ((0., 0.), (5., -5.), (1000., 2000.)):
            with self.subTest(origin=(ox, oy)):
                xe = ox + np.linspace(-dx/4, dx, 19)
                ye = oy + np.linspace(-dy/2, dy/4, 13)
                kernel, _ = wave_to_source_kernel(unit_cell((ox, oy), dx, dy), PixelGrid(xe, ye))
                expected = independent_rectangle_overlap(xe-ox, ye-oy, dx, dy)
                np.testing.assert_allclose(kernel.probability_mass, expected,
                                           atol=2e-15, rtol=1e-12)
                self.assertAlmostEqual(kernel.captured_fraction, .75**2, delta=3e-11)

    def test_affine_tilted_anisotropic_cell_half_plane_and_full_capture(self):
        theta, phi = .45, .7
        chief = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
        first = np.array([1., 0., 0.]) - chief[0]*chief
        first /= np.linalg.norm(first)
        basis = np.array([first, np.cross(chief, first), chief])
        np.testing.assert_allclose(basis@basis.T, np.eye(3), atol=3e-16)
        self.assertGreater(np.linalg.det(basis), .999999999999)
        dx, dy = .007, .013
        for ox, oy in ((0., 0.), (5., -7.), (1000., -1200.)):
            with self.subTest(origin=(ox, oy)):
                source = unit_cell((ox, oy), dx, dy, basis)
                full = PixelGrid(ox+np.linspace(-.04, .04, 31),
                                 oy+np.linspace(-.04, .04, 35))
                cropped = PixelGrid(ox+np.linspace(0., .04, 16),
                                    oy+np.linspace(-.04, .04, 35))
                kfull, _ = wave_to_source_kernel(source, full, allow_local_detector_projection=True)
                khalf, _ = wave_to_source_kernel(source, cropped, allow_local_detector_projection=True)
                # Every nonsingular linear image of a centered rectangle is
                # centrosymmetric. Any half-plane through its center has 1/2 area.
                self.assertAlmostEqual(kfull.captured_fraction, 1., delta=2e-14)
                self.assertAlmostEqual(khalf.captured_fraction, .5, delta=2e-14)

    def test_multiple_fft_cells_keep_independent_cell_weights(self):
        dx, dy = .01, .004
        source = unit_cell((5., -4.), dx, dy)
        wave = source['wave']
        wave['x_mm'] = np.array([-dx, 0., dx])
        wave['y_mm'] = np.array([-dy, 0., dy])
        weights = np.arange(1., 10.).reshape(3, 3) / 45.
        wave['psf'] = weights/(dx*dy)
        xe = 5. + np.linspace(-dx*1.1, dx*1.5, 23)
        ye = -4. + np.linspace(-dy*1.5, dy*.6, 29)
        kernel, _ = wave_to_source_kernel(source, PixelGrid(xe, ye))
        expected = np.zeros(kernel.probability_mass.shape)
        for j, y in enumerate(wave['y_mm']):
            for i, x in enumerate(wave['x_mm']):
                expected += weights[j, i]*independent_rectangle_overlap(xe-5.-x, ye+4.-y, dx, dy)
        np.testing.assert_allclose(kernel.probability_mass, expected, atol=2e-15, rtol=1e-12)
        self.assertLess(kernel.captured_fraction, 1.)

    def test_area_does_not_subtract_large_global_products(self):
        original = np.array([[-.003, -.001], [.003, -.001], [.003, .001], [-.003, .001]])
        for offset in (0., 5., 1000.):
            polygon = original+offset
            # Use the represented side lengths as oracle; translation cannot
            # restore bits already lost when the caller forms global vertices.
            area = (polygon[1, 0]-polygon[0, 0])*(polygon[2, 1]-polygon[1, 1])
            self.assertAlmostEqual(_polygon_area(polygon), area, delta=1e-20)
            self.assertAlmostEqual(_polygon_area(polygon[::-1]), area, delta=1e-20)

    def test_invalid_detector_edges_and_singular_projection_remain_rejected(self):
        for edges in ([0., 1., .5], [0., 0., 1.], [0., np.nan, 1.]):
            with self.assertRaises(InputError):
                PixelGrid(edges, [0., 1.])
        basis = np.array([[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]])
        with self.assertRaises(ModelApplicabilityError):
            wave_to_source_kernel(unit_cell((0., 0.), basis=basis),
                                  PixelGrid([-1., 1.], [-1., 1.]),
                                  allow_local_detector_projection=True)

    def test_actual_producer_translated_overlap_preserves_entire_fft_power(self):
        system = OpticalSystem((Surface(0., 0., 5.),),
                               {'air': Material('air', n_constant=1.)},
                               object_z=-100., image_z=20.)
        settings = EngineSettings(sampling=SamplingSettings(
            pupil_rings=3, pupil_arms=16, wave_grid_size=16, fft_padding=2))
        result = source_diffraction_kernel(system, [25., 25.], .01, 1.,
                                          settings=settings, frequencies=[0.],
                                          transmission_mode='lossless_scalar')
        wave = result['wave']
        original_origin = wave['image_reference_point_mm'][:2]
        corners = np.array([[x, y]
            for x in (wave['x_mm'][0]-wave['dx_mm']/2, wave['x_mm'][-1]+wave['dx_mm']/2)
            for y in (wave['y_mm'][0]-wave['dy_mm']/2, wave['y_mm'][-1]+wave['dy_mm']/2)])
        projected = corners@np.linalg.inv(wave['output_basis'][:2, :2]).T
        self.assertTrue(np.all(np.abs(projected) < .02))
        local_edges = np.linspace(-.02, .02, 41)
        reference_mass = None
        for origin in (original_origin, np.zeros(2)):
            moved = dict(result, wave=dict(wave, image_reference_point_mm=np.r_[origin, 20.]))
            grid = PixelGrid(origin[0]+local_edges, origin[1]+local_edges)
            kernel, _ = wave_to_source_kernel(moved, grid, allow_local_detector_projection=True)
            self.assertAlmostEqual(kernel.captured_fraction, 1., delta=3e-14)
            if reference_mass is None:
                reference_mass = kernel.probability_mass
            else:
                np.testing.assert_allclose(kernel.probability_mass, reference_mass,
                                           atol=1e-13, rtol=1e-10)
        # This is an area/registration regression. It does not promote the
        # producer's phase sampling status to validated optical accuracy.


if __name__ == '__main__':
    unittest.main()
