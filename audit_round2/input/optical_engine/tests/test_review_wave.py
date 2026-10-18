"""Independent geometry-author review of diffraction: analytic physical oracles."""
import unittest
import numpy as np
from sequential_optics.models import Surface, OpticalSystem, Field
from sequential_optics.materials import Material
from sequential_optics.geometry import aim_ray, trace_ray, trace_bundle
from sequential_optics.wave import (direct_rs, afocal_fft, fft_psf, sample_axis_otf,
                                    make_boundary_state, VACUUM_IMPEDANCE_OHM)


class IndependentDiffractionReviewTests(unittest.TestCase):
    def test_RS_uniform_disk_axis_closed_form_complex_field(self):
        # Integrating the RS-I kernel over an exactly uniform disk gives
        # exp(ikL) - L/sqrt(L²+a²)*exp(ik sqrt(L²+a²)). This independent
        # closed form checks phase, absolute amplitude and obliquity together.
        radius=1.3;distance=12.;wavelength=.05;index=1.23
        nodes,weights=np.polynomial.legendre.leggauss(100)
        q=(nodes+1)*radius**2/2
        rho=np.sqrt(q)
        # On-axis cylindrical symmetry makes one azimuth sample sufficient;
        # the exact annulus-area weight accounts for the whole 2pi rotation.
        area=np.pi*radius**2/2*weights
        actual=direct_rs(np.ones(len(q)),rho,np.zeros(len(q)),[0,0,distance],
                         wavelength,index=index,area_weights_mm2=area)[0]
        k=2*np.pi*index/wavelength;r=np.sqrt(distance**2+radius**2)
        expected=np.exp(1j*k*distance)-distance/r*np.exp(1j*k*r)
        self.assertLess(abs(actual-expected),2e-11)

    def test_afocal_raw_power_has_actual_oblique_cosine(self):
        # A discrete periodic plane wave occupies one exact spectral mode.
        # Physical flux through the plane is n/(2Z0)|E|² area cos(theta).
        n=1.33;wavelength=.00055;medium=wavelength/n
        count=64;du=medium;mode=20
        coords=(np.arange(count)-count//2)*du
        xx,yy=np.meshgrid(coords,coords)
        l=mode/count;gamma=np.sqrt(1-l*l)
        field=np.exp(2j*np.pi*l*xx/medium)
        output=afocal_fft(field,du,wavelength,index=n,padding=1)
        expected=n/(2*VACUUM_IMPEDANCE_OHM)*(count*du)**2*1e-6*gamma
        self.assertAlmostEqual(output['electric_field_wave_power_w']/expected,1.,places=12)
        flat=afocal_fft(np.ones_like(field),du,wavelength,index=n,padding=1)
        self.assertAlmostEqual(output['electric_field_wave_power_w']/flat['electric_field_wave_power_w'],gamma,places=12)

    def test_focal_frequency_scale_uses_medium_wavelength(self):
        count=96;d=.0025;x=(np.arange(count)-count//2)*d
        xx,yy=np.meshgrid(x,x);pupil=(xx*xx+yy*yy<.1**2)
        vacuum=fft_psf(pupil,d,.00055,index=1.,padding=4)
        immersed=fft_psf(pupil,d,.00055,index=1.5,padding=4)
        self.assertAlmostEqual(vacuum['dx_mm']/immersed['dx_mm'],1.5,places=12)
        f=np.array([30.,60.,120.])
        np.testing.assert_allclose(sample_axis_otf(immersed,1.5*f),sample_axis_otf(vacuum,f),atol=2e-12)

    def test_refracted_OPL_gradient_matches_ray_direction_on_boundary(self):
        materials={'air':Material('air',n_constant=1.),'glass':Material('glass',n_constant=1.5)}
        surfaces=(Surface(0,.02,4,medium_after='glass'),Surface(3,-.02,4))
        for field,obj in [(Field(.04,.02),None),(Field(.7,.2,'height'),-40.)]:
            system=OpticalSystem(surfaces,materials,object_z=obj,stop_index=1,stop_radius=2.)
            bundle=trace_bundle(system,field,rings=3,arms=12)
            state=make_boundary_state(system,bundle,mode='afocal')
            def boundary(pupil):
                tr=trace_ray(system,aim_ray(system,field,pupil))
                self.assertTrue(tr.success,tr.message)
                d=tr.last_direction
                t=np.dot(state.origin-tr.last_position,state.basis[2])/np.dot(d,state.basis[2])
                p=tr.last_position+t*d
                return state.basis[:2]@(p-state.origin),tr.last_opl+state.output_index*t,d
            target=np.array([.3,-.2]);_,_,direction=boundary(target)
            expected=state.output_index*(state.basis[:2]@direction)
            errors=[]
            for h in [2e-3,1e-3]:
                j=np.empty((2,2));ds=np.empty(2)
                for axis in range(2):
                    offset=np.eye(2)[axis]*h
                    up,sp,_=boundary(target+offset);um,sm,_=boundary(target-offset)
                    j[:,axis]=(up-um)/(2*h);ds[axis]=(sp-sm)/(2*h)
                actual=np.linalg.solve(j.T,ds)
                errors.append(np.linalg.norm(actual-expected))
            self.assertLess(errors[-1],2e-7)
            self.assertLess(errors[-1],errors[0]+1e-9)


if __name__=='__main__': unittest.main()
