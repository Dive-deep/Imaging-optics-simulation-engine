"""Analytic and counterexample-oriented scalar wave tests (no vendor oracle)."""
import unittest
from dataclasses import replace
import numpy as np
from sequential_optics.wave import (fft_psf,otf_from_psf,sample_axis_otf,circular_mtf,
    afocal_fft,angular_otf,angular_chart_jacobian,direct_rs,focal_direct_reference,
    analyze_wave,make_boundary_state,ModelApplicabilityError,InputError)
from sequential_optics.models import Surface,OpticalSystem,Field,RayBundle,TraceResult
from sequential_optics.settings import EngineSettings,SamplingSettings


class ConstantIndex:
    def index(self,wavelength_nm,temperature_c=None): return 1.


def grid(n=128,extent=.12):
    d=2*extent/n
    x=(np.arange(n)-n//2)*d
    xx,yy=np.meshgrid(x,x)
    return x,xx,yy,d


def perfect_bundle(radius=2.,focus_z=50.,rings=10,arms=48,tilt=0.):
    # A prescribed exact spherical eikonal at a transparent physical plane.
    # This is a boundary construction oracle, not a purported manufacturable lens.
    xy=[[0.,0.]]
    for ir in range(1,rings+1):
        rr=radius*ir/rings
        for theta in np.arange(arms)*2*np.pi/arms:
            xy.append([rr*np.cos(theta),rr*np.sin(theta)])
    xy=np.asarray(xy)
    traces=[]
    q=np.array([focus_z*np.tan(tilt),0.,focus_z])
    chief_r=np.linalg.norm(q)
    for x,y in xy:
        p=np.array([x,y,0.]);delta=q-p;R=np.linalg.norm(delta)
        traces.append(TraceResult(np.array([p]),np.array([delta/R]),np.array([chief_r-R]),np.array([1.])))
    system=OpticalSystem((Surface(0.,0.,radius),),{'air':ConstantIndex()},image_z=focus_z)
    bundle=RayBundle(tuple(traces),xy/radius,np.ones(len(xy))/len(xy),Field(),550.,traces[0])
    return system,bundle


class FourierTests(unittest.TestCase):
    def test_circle_analytic_mtf_and_physical_scale(self):
        x,xx,yy,d=grid(192,.12)
        pupil=(xx*xx+yy*yy<=.1**2).astype(float)
        result=fft_psf(pupil,d,.00055,padding=4)
        frequencies=np.linspace(0,350,51)
        got=np.abs(sample_axis_otf(result,frequencies))
        expected=circular_mtf(frequencies,2*.1/.00055)
        self.assertLess(np.max(np.abs(got-expected)),.006)
        self.assertAlmostEqual(result['psf'].sum()*result['dx_mm']*result['dy_mm'],1.,places=12)
        self.assertAlmostEqual(result['raw_integral'],.00055**2*result['pupil_squared_integral'],places=16)

    def test_rectangular_pupil_separable_otf_axes(self):
        x,xx,yy,d=grid(192,.12)
        p=((abs(xx)<.08)&(abs(yy)<.04)).astype(float)
        result=fft_psf(p,d,.00055,padding=4)
        f=np.array([20.,50.,90.])
        ix=np.flatnonzero(np.any(p,axis=0));iy=np.flatnonzero(np.any(p,axis=1))
        cutoff_x=len(ix)*d/.00055;cutoff_y=len(iy)*d/.00055
        self.assertTrue(np.allclose(sample_axis_otf(result,f,'x').real,np.maximum(0,1-f/cutoff_x),atol=2e-12))
        self.assertTrue(np.allclose(sample_axis_otf(result,f,'y').real,np.maximum(0,1-f/cutoff_y),atol=2e-12))

    def test_tilt_translation_and_complex_otf_sign(self):
        x,xx,yy,d=grid(96,.12);mask=(xx*xx+yy*yy<=.1**2)
        base=fft_psf(mask,d,.00055,padding=4)
        shift=5*base['dx_mm']
        tilted=fft_psf(mask*np.exp(2j*np.pi*xx*shift/.00055),d,.00055,padding=4)
        self.assertTrue(np.allclose(tilted['psf'],np.roll(base['psf'],5,axis=1),rtol=1e-10,atol=1e-9))
        fx=base['frequency_x'];expected=base['otf']*np.exp(-2j*np.pi*fx[None,:]*shift)
        self.assertLess(np.max(abs(tilted['otf']-expected)),1e-12)
        self.assertLess(np.max(abs(tilted['mtf']-base['mtf'])),1e-12)

    def test_piston_amplitude_and_defocus(self):
        x,xx,yy,d=grid(96,.12);mask=(xx*xx+yy*yy<=.1**2)
        a=fft_psf(mask,d,.00055,padding=4)
        b=fft_psf(2*mask*np.exp(.7j),d,.00055,padding=4)
        self.assertAlmostEqual(b['raw_integral']/a['raw_integral'],4.,places=12)
        self.assertLess(np.max(abs(a['psf']-b['psf'])),1e-9)
        defocused=fft_psf(mask*np.exp(2j*np.pi*2*(xx*xx+yy*yy)/.1**2),d,.00055,padding=4)
        self.assertLess(defocused['psf'].max(),a['psf'].max()/5)
        self.assertLess(np.abs(sample_axis_otf(defocused,[80.]))[0],np.abs(sample_axis_otf(a,[80.]))[0])

    def test_clipped_and_annular_diffraction_keep_actual_pupil(self):
        x,xx,yy,d=grid(128,.12);rho=np.sqrt(xx*xx+yy*yy)
        mask=(rho<=.1)&(rho>=.05)&(xx>-.04)
        r=fft_psf(mask,d,.00055,padding=4)
        f=np.array([30.,70.,120.])
        self.assertGreater(np.max(abs(np.abs(sample_axis_otf(r,f))-circular_mtf(f,2*.1/.00055))),.05)
        self.assertGreater(np.max(abs(np.abs(sample_axis_otf(r,f,'x'))-np.abs(sample_axis_otf(r,f,'y')))),.02)
        self.assertTrue(np.allclose(r['otf'],fft_psf(abs(mask),d,.00055,padding=4)['otf']))

    def test_negative_otf_is_not_erased(self):
        n=64;psf=np.zeros((n,n));psf[n//2,n//2-4]=.5;psf[n//2,n//2+4]=.5
        r=otf_from_psf(psf,.01)
        value=sample_axis_otf(r,[12.5])[0]
        self.assertLess(value.real,-.999999)
        self.assertAlmostEqual(abs(value),1.,places=12)

    def test_padding_is_not_pupil_sampling_or_larger_window(self):
        x,xx,yy,d=grid(64,.12);p=(xx*xx+yy*yy<.1**2)
        a=fft_psf(p,d,.00055,padding=2);b=fft_psf(p,d,.00055,padding=4)
        self.assertEqual(a['input_grid_shape'],b['input_grid_shape'])
        self.assertEqual(a['physical_period_mm'],b['physical_period_mm'])
        self.assertAlmostEqual(a['dx_mm']/b['dx_mm'],2.)
        self.assertAlmostEqual(a['raw_integral'],b['raw_integral'],places=16)
        self.assertTrue(b['linear_autocorrelation_supported'])
        c=fft_psf(p,d,.00055,padding=1)
        self.assertFalse(c['linear_autocorrelation_supported'])

    def test_direct_rs_reference_same_boundary_near_focus(self):
        x,xx,yy,d=grid(100,.06);p=(xx*xx+yy*yy<.05**2).astype(complex)
        points=np.array([[0.,0.],[.001,0.],[0.,.002]])
        exact=focal_direct_reference(p,x,x,50.,.00055,points)
        approx=np.array([np.sum(p*np.exp(-2j*np.pi/.00055*(qx*xx+qy*yy)))*d*d/(1j*.00055) for qx,qy in points])
        self.assertTrue(np.allclose(abs(exact)**2,abs(approx)**2,rtol=1e-4))

    def test_direct_rs_grid_and_nonuniform_weights_agree(self):
        u=np.linspace(-1,1,40,endpoint=False);v=np.linspace(-.5,.5,32,endpoint=False)
        uu,vv=np.meshgrid(u,v);p=np.exp(-(uu**2+vv**2)*3)
        target=np.array([[0,0,10.],[.01,.02,10.]])
        a=direct_rs(p,u,v,target,.00055)
        b=direct_rs(p,uu,vv,target,.00055,area_weights_mm2=(u[1]-u[0])*(v[1]-v[0]))
        self.assertTrue(np.allclose(a,b))

    def test_bad_inputs_are_explicit(self):
        with self.assertRaises(InputError): fft_psf(np.zeros((8,8)),.01,.00055)
        with self.assertRaises(InputError): fft_psf(np.ones((7,8)),.01,.00055)
        r=fft_psf(np.ones((8,8)),.01,.00055)
        with self.assertRaises(ModelApplicabilityError): sample_axis_otf(r,[1e9])
        self.assertEqual(sample_axis_otf(r,[1e9],support_cutoff=100.)[0],0j)


class AfocalTests(unittest.TestCase):
    def test_solid_angle_and_angle_density_jacobian(self):
        # A wavelength/grid ratio deliberately broad enough to distinguish
        # direction cosines, angles, and solid-angle measure.
        x,xx,yy,d=grid(64,.001)
        p=np.exp(-((xx/.0003)**2+(yy/.0005)**2))
        r=afocal_fft(p,d,.00055,padding=4)
        valid=r['propagating_mask'];g=r['gamma']
        integration=np.sum(r['h_solid_angle'][valid]/g[valid])*r['dl']*r['dm']
        self.assertAlmostEqual(integration,1.,places=12)
        self.assertAlmostEqual(r['cell_power'].sum(),1.,places=12)
        self.assertGreater(abs(np.sum(r['h_solid_angle'])*r['dl']*r['dm']-1),.01)
        self.assertAlmostEqual(angular_otf(r,[0.])[0].real,1.,places=12)
        f=np.array([1.,2.]);actual=angular_otf(r,f)
        ll=np.broadcast_to(r['l'][None,:],r['cell_power'].shape)
        false=np.array([np.sum(r['cell_power']*np.exp(-2j*np.pi*v*ll)) for v in f])
        self.assertGreater(np.max(abs(actual-false)),.001)
        self.assertTrue(np.allclose(r['h_angle_density'],r['h_solid_angle']*r['angular_chart_jacobian']))

    def test_chart_jacobian_against_direction_derivatives(self):
        a=.4;b=-.3;h=1e-6
        def direction(a,b):
            d=np.array([np.tan(a),np.tan(b),1.]);return d/np.linalg.norm(d)
        da=(direction(a+h,b)-direction(a-h,b))/(2*h)
        db=(direction(a,b+h)-direction(a,b-h))/(2*h)
        self.assertAlmostEqual(angular_chart_jacobian(a,b),np.linalg.norm(np.cross(da,db)),places=9)

    def test_raw_electric_power_scaling_not_only_normalized_psf(self):
        u=np.arange(-32,32)*.02;xx,yy=np.meshgrid(u,u);p=np.exp(-(xx*xx+yy*yy)/.1**2)
        a=afocal_fft(p,.02,.00055,padding=4);b=afocal_fft(2*p,.02,.00055,padding=4)
        self.assertAlmostEqual(b['electric_field_wave_power_w']/a['electric_field_wave_power_w'],4.,places=12)
        self.assertTrue(np.allclose(a['h_solid_angle'],b['h_solid_angle']))
        boundary_power=np.sum(p*p)*.02**2*1e-6/(2*376.730313668)
        self.assertAlmostEqual(a['electric_field_wave_power_w']/boundary_power,1.,places=5)


class RayBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.settings=EngineSettings(sampling=SamplingSettings(wave_grid_size=48,fft_padding=4))

    def test_exact_spherical_eikonal_zero_opd_on_canonical_plane(self):
        sys,bundle=perfect_bundle()
        result=analyze_wave(sys,bundle,self.settings,direct_check_points=[[0.,0.]])
        self.assertLess(np.max(abs(result['pupil_opd_mm'][result['pupil_mask']])),3e-14)
        self.assertTrue(np.allclose(result['mtf_s'],result['diffraction_mtf_s'],atol=1e-12))
        self.assertEqual(result['removed_phase_terms'],['piston'])
        self.assertEqual(result['vendor_equivalence'],'unverified')
        self.assertTrue(result['sampling_diagnostics']['phase_sampling_necessary_check_passed'])
        self.assertGreater(result['boundary_origin_mm'][2],sys.surfaces[-1].z)
        coarse=analyze_wave(sys,bundle,self.settings,padding=1,frequencies=[0.,10.])
        self.assertEqual(coarse['status'],'underresolved')
        self.assertFalse(coarse['sampling_diagnostics']['linear_autocorrelation_sampling_check_passed'])

    def test_fixed_boundary_amplitude_focus_sweep(self):
        sys,bundle=perfect_bundle(radius=5.,focus_z=30.)
        baseline=analyze_wave(sys,bundle,self.settings)
        shifted=analyze_wave(sys,bundle,self.settings,image_z=30.2,boundary_state=baseline['boundary_state'])
        reset=analyze_wave(sys,bundle,self.settings,image_z=30.2,boundary_state=baseline['boundary_state'],illumination_policy='uniform_fourier_pupil_each_focus')
        self.assertIs(shifted['boundary_state'],baseline['boundary_state'])
        self.assertGreater(np.max(abs(shifted['pupil_opd_mm'][shifted['pupil_mask']])),1e-4)
        self.assertGreater(np.ptp(abs(shifted['pupil'][shifted['pupil_mask']])),1e-5)
        self.assertLess(np.ptp(abs(reset['pupil'][reset['pupil_mask']])),1e-12)
        self.assertLess(shifted['psf'].max(),baseline['psf'].max())
        with self.assertRaises(InputError):
            analyze_wave(sys,bundle,self.settings,mode='afocal',boundary_state=baseline['boundary_state'])

    def test_image_before_full_surface_boundary_is_rejected(self):
        sys,bundle=perfect_bundle()
        with self.assertRaises(ModelApplicabilityError):
            analyze_wave(sys,bundle,self.settings,image_z=.0001)

    def test_numerical_failure_not_clipping(self):
        sys,bundle=perfect_bundle()
        traces=list(bundle.traces);traces[2]=TraceResult(status='aim_failed')
        broken=replace(bundle,traces=tuple(traces))
        with self.assertRaises(ModelApplicabilityError): analyze_wave(sys,broken,self.settings)

    def test_clipped_triangle_preserves_hole_in_original_topology(self):
        sys,bundle=perfect_bundle(rings=8,arms=32)
        traces=list(bundle.traces)
        # Remove interior region away from chief; successful-only Delaunay
        # would silently bridge this physical aperture hole.
        uv=bundle.pupil_xy
        holes=(uv[:,0]>.25)&(uv[:,0]<.65)&(abs(uv[:,1])<.2)
        for i in np.flatnonzero(holes): traces[i]=TraceResult(status='aperture_blocked')
        clipped=replace(bundle,traces=tuple(traces))
        a=analyze_wave(sys,bundle,self.settings);b=analyze_wave(sys,clipped,self.settings)
        self.assertLess(b['pupil_mask'].sum(),a['pupil_mask'].sum()*.98)
        self.assertTrue(np.all(b['pupil_mask']<=a['pupil_mask']))

    def test_folded_boundary_map_rejected(self):
        sys,bundle=perfect_bundle(rings=5,arms=24)
        traces=list(bundle.traces)
        i=1+2*24+3
        t=traces[i];p=t.positions.copy();p[0,:2]*=-1
        traces[i]=TraceResult(p,t.directions,t.opl,t.powers)
        with self.assertRaises(ModelApplicabilityError):
            make_boundary_state(sys,replace(bundle,traces=tuple(traces)),self.settings)

    def test_afocal_uniform_plane_and_no_ima(self):
        sys,bundle=perfect_bundle(rings=6,arms=32)
        traces=[]
        for t in bundle.traces:
            traces.append(TraceResult(t.positions,np.array([[0.,0.,1.]]),np.array([0.]),np.array([1.])))
        b=replace(bundle,traces=tuple(traces),chief=traces[0])
        s=replace(sys,image_z=None)
        r=analyze_wave(s,b,self.settings,mode='afocal',frequencies=[0.,100.,1000.])
        self.assertAlmostEqual(r['cell_power'].sum(),1.,places=12)
        self.assertTrue(np.allclose(r['h_solid_angle'],r['same_pupil_diffraction_h_solid_angle']))
        self.assertEqual(r['mode'],'afocal')
        self.assertAlmostEqual(r['mtf_s'][0],1.,places=12)

    def test_full_surface_bound_does_not_use_sparse_ray_max(self):
        sys,bundle=perfect_bundle(radius=2.)
        curved=replace(sys.surfaces[-1],surface_type='even_asphere',coefficients=(.2,))
        altered=replace(sys,surfaces=(curved,))
        state=make_boundary_state(altered,bundle,self.settings)
        self.assertGreater(state.origin[2],.2*curved.clear_radius**2)

    def test_actual_trace_oblique_plane_wave_preserves_incident_phase(self):
        from sequential_optics.geometry import trace_bundle
        from sequential_optics.materials import Material
        # A transparent finite-aperture plate must not invent wavefront error
        # when unequal launch coordinates are illuminated by an oblique wave.
        air=Material('only G known',wavelengths_nm=(550.,),indices=(1.,))
        system=OpticalSystem((Surface(0.,0.,2.),Surface(1.,0.,2.)),{'air':air})
        bundle=trace_bundle(system,Field(.08,.03),settings=self.settings,rings=5,arms=24)
        result=analyze_wave(system,bundle,self.settings,mode='afocal')
        self.assertLess(np.max(np.abs(result['pupil_opd_mm'][result['pupil_mask']])),2e-12)
        self.assertTrue(np.allclose(result['h_solid_angle'],result['same_pupil_diffraction_h_solid_angle'],rtol=1e-8,atol=1e-8))

    def test_actual_refracting_lens_single_G_index_only(self):
        from sequential_optics.geometry import trace_bundle
        from sequential_optics.paraxial import paraxial_properties
        from sequential_optics.materials import Material
        air=Material('air G',wavelengths_nm=(550.,),indices=(1.,))
        glass=Material('synthetic G',wavelengths_nm=(550.,),indices=(1.5,))
        system=OpticalSystem((Surface(0.,1/50.,1.,medium_after='glass'),Surface(3.,-1/50.,1.)),{'air':air,'glass':glass})
        p=paraxial_properties(system,settings=self.settings)
        system=replace(system,image_z=system.surfaces[-1].z+p['bfl_mm'])
        bundle=trace_bundle(system,Field(),settings=self.settings,rings=5,arms=24)
        result=analyze_wave(system,bundle,self.settings)
        self.assertEqual(bundle.wavelength_nm,550.)
        self.assertTrue(np.all(np.isfinite(result['mtf_s'])))
        self.assertAlmostEqual(result['mtf_s'][0],1.,places=12)
        self.assertGreater(result['sampling_diagnostics']['successful_ray_count'],100)


if __name__=='__main__': unittest.main()
