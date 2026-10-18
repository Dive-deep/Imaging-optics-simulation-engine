"""Referee Fourier counterexamples with independent finite-grid analytic answers.

They test the declared scalar/discrete profile, not vendor equivalence or vector
diffraction. The alias specimen deliberately demonstrates a missing-sampling
ambiguity: identical aliased samples cannot be detected from wrapped data alone.
"""
import unittest
import numpy as np
from sequential_optics.wave import fft_psf, afocal_fft, angular_otf


class IndependentWaveTests(unittest.TestCase):
    def test_rectangular_pupil_exact_discrete_triangular_otf_both_axes(self):
        pupil = np.zeros((48, 64), complex)
        pupil[16:32, 22:42] = 1
        dx, dy, wavelength, index = .002, .003, .00055, 1.2
        result = fft_psf(pupil, (dx, dy), wavelength, index=index, padding=2)
        lm = wavelength/index
        tx = np.maximum(0, 1-np.abs(result['frequency_x'])*lm/(20*dx))
        ty = np.maximum(0, 1-np.abs(result['frequency_y'])*lm/(16*dy))
        expected = ty[:, None]*tx[None, :]
        np.testing.assert_allclose(result['otf'].real, expected, atol=2e-13, rtol=0)
        np.testing.assert_allclose(result['otf'].imag, 0, atol=2e-13, rtol=0)

    def test_raw_parseval_scale_and_amplitude_squared_before_normalization(self):
        y, x = np.mgrid[-16:16, -20:20]
        pupil = np.exp(-(x*x/85+y*y/49))*np.exp(0.1j*x+0.005j*y*y)
        spacing, wavelength, index = (.0015, .0021), .000532, 1.4
        base = fft_psf(pupil, spacing, wavelength, index=index, padding=3)
        twice = fft_psf(2*pupil, spacing, wavelength, index=index, padding=3)
        expected_raw = (wavelength/index)**2*np.sum(np.abs(pupil)**2)*np.prod(spacing)
        self.assertAlmostEqual(base['raw_integral']/expected_raw, 1, places=12)
        self.assertAlmostEqual(twice['raw_integral']/base['raw_integral'], 4, places=12)
        np.testing.assert_allclose(twice['psf'], base['psf'], rtol=1e-13, atol=0)
        np.testing.assert_allclose(twice['otf'], base['otf'], rtol=1e-13, atol=1e-14)
        far = afocal_fft(pupil, (.004, .007), wavelength, index=index, padding=3)
        far2 = afocal_fft(2*pupil, (.004, .007), wavelength, index=index, padding=3)
        self.assertAlmostEqual(far2['electric_field_wave_power_w']/far['electric_field_wave_power_w'], 4, places=12)
        np.testing.assert_allclose(far2['cell_power'], far['cell_power'], atol=1e-14, rtol=1e-13)

    def test_single_dft_direction_uses_exact_angle_not_direction_cosine(self):
        # Periodic single DFT mode is an algebra fixture, not a finite physical
        # aperture with padding=1 claimed to have a perfect delta far field.
        n, wavelength, dl = 64, .00055, .025
        du = wavelength/(n*dl)
        coordinate = (np.arange(n)-n//2)*du
        u, v = np.meshgrid(coordinate, coordinate)
        l0, m0 = .5, .25
        field = np.exp(2j*np.pi*(l0*u+m0*v)/wavelength)
        result = afocal_fft(field, du, wavelength, padding=1)
        gamma0 = np.sqrt(1-l0*l0-m0*m0)
        alpha0, beta0 = np.arctan2(l0, gamma0), np.arctan2(m0, gamma0)
        frequencies = np.array([0., .5, 1.3])
        np.testing.assert_allclose(angular_otf(result, frequencies, axis='x'),
                                   np.exp(-2j*np.pi*frequencies*alpha0), atol=5e-13)
        np.testing.assert_allclose(angular_otf(result, frequencies, axis='y'),
                                   np.exp(-2j*np.pi*frequencies*beta0), atol=5e-13)
        valid = result['propagating_mask']
        solid_angle_integral = np.sum(result['h_solid_angle'][valid]
                                      /result['gamma'][valid])*result['dl']*result['dm']
        self.assertAlmostEqual(solid_angle_integral, 1, places=12)
        self.assertGreater(abs(alpha0-l0), .04)

    def test_commensurate_refinement_can_hide_phase_alias(self):
        wavelength, width, phase_cycles = .00055, .2, 64
        peaks = []
        for n in (32, 64, 70, 160):
            coordinate = (np.arange(n)-n//2)*width/n
            phase = np.exp(2j*np.pi*phase_cycles*coordinate/width)
            pupil = np.broadcast_to(phase, (8, n))
            result = fft_psf(pupil, (width/n, .005), wavelength, padding=2)
            peak = np.unravel_index(np.argmax(result['psf']), result['psf'].shape)
            peaks.append(result['x_mm'][peak[1]])
        # Two apparently stable refinements report the wrong zero shift.
        self.assertAlmostEqual(peaks[0], 0, places=12)
        self.assertAlmostEqual(peaks[1], 0, places=12)
        # A noncommensurate grid exposes instability but does not yet resolve it.
        self.assertGreater(abs(peaks[2]), .005)
        # Only the adequate unwrapped phase bandwidth gives the true shift.
        self.assertAlmostEqual(peaks[3], wavelength*phase_cycles/width, places=12)


if __name__ == '__main__':
    unittest.main()
