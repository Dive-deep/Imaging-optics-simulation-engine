"""Independent C1/C2 physical controls; deliberately no production edits."""
from pathlib import Path
import json
import sys
from dataclasses import replace
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sequential_optics import Surface, OpticalSystem, Field
from sequential_optics.geometry import aim_ray, trace_ray, validate_geometry
from sequential_optics.imaging import PixelGrid
from sequential_optics.materials import Material
from sequential_optics.radiometry import trace_transmission
from sequential_optics.source_wave import wave_to_source_kernel


def independent_vertical_integral(J, width, height, xb, yb):
    """Area from y-interval integration, not polygon clipping/area code.

    The transformed rectangle satisfies two linear strip inequalities.
    At x, intersect their y intervals with the detector y interval. Split
    quadrature at intersections of bound lines: the integrand is piecewise
    linear and each segment is integrated independently.
    """
    lower = [(0., yb[0])]
    upper = [(0., yb[1])]
    left, right = xb
    for (a, b), extent in zip(J, (width, height)):
        if abs(b) < 1e-14:
            if abs(a) < 1e-14:
                continue
            q = sorted((-extent / (2*a), extent / (2*a)))
            left, right = max(left, q[0]), min(right, q[1])
        else:
            m = -a/b
            q = sorted((-extent/(2*b), extent/(2*b)))
            lower.append((m, q[0]))
            upper.append((m, q[1]))
    if left >= right:
        return 0.
    lines = lower + upper
    breaks = [left, right]
    for i, (a, b) in enumerate(lines):
        for c, d in lines[i+1:]:
            if a != c:
                x = (d-b)/(a-c)
                if left < x < right:
                    breaks.append(x)
    breaks = sorted(set(breaks))
    f = lambda x: max(0., min(a*x+b for a,b in upper)-max(a*x+b for a,b in lower))
    return sum(quad(f, a, b, epsabs=1e-22, epsrel=1e-12)[0]
               for a,b in zip(breaks,breaks[1:]))


def detector_controls():
    rng = np.random.default_rng(91263)
    comparisons = []
    for index in range(12):
        theta, phi = rng.uniform(.05, 1.1), rng.uniform(-np.pi, np.pi)
        c = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        s = np.array([1.,0.,0.])-c[0]*c
        s /= np.linalg.norm(s)
        basis = np.array([s,np.cross(c,s),c])
        np.testing.assert_allclose(basis@basis.T,np.eye(3),atol=6e-16)
        J = basis[:2,:2]
        assert np.linalg.det(J)>0
        dx, dy = rng.uniform(.002,.015,2)
        shift = np.array([5.,-7.]) if index%2 else np.array([1000.,-1200.])
        # The cropped rectangle is deliberately asymmetric, avoiding a
        # half-plane symmetry oracle and preventing a hidden normalize-to-one.
        xe = shift[0]+np.linspace(-.31*dx,.74*dx,13)
        ye = shift[1]+np.linspace(-.65*dy,.19*dy,11)
        source = {'wavelength_nm':550.,'wave':dict(
            mode='focal',illumination='physical_source_triangle_boundary_amplitude',
            output_basis=basis,image_reference_point_mm=np.r_[shift,20.],
            dx_mm=dx,dy_mm=dy,x_mm=np.array([0.]),y_mm=np.array([0.]),
            psf=np.array([[1/(dx*dy)]]))}
        grid = PixelGrid(xe,ye)
        kernel,_ = wave_to_source_kernel(source,grid,allow_local_detector_projection=True)
        expected = np.empty(grid.shape)
        lx,ly=xe-shift[0],ye-shift[1]
        area=dx*dy/np.linalg.det(J)
        for iy in range(len(ly)-1):
            for ix in range(len(lx)-1):
                expected[iy,ix]=independent_vertical_integral(J,dx,dy,lx[ix:ix+2],ly[iy:iy+2])/area
        error=float(np.max(abs(expected-kernel.probability_mass)))
        np.testing.assert_allclose(kernel.probability_mass,expected,atol=3e-16,rtol=2e-12)
        assert .1 < kernel.captured_fraction < .9
        comparisons.append(dict(theta_rad=float(theta),phi_rad=float(phi),
            reference_captured=float(expected.sum()),actual_captured=kernel.captured_fraction,
            max_pixel_error=error,convex=True,projection_determinant=float(np.linalg.det(J))))
    return comparisons


def coincident_power_controls():
    mats={'air':Material('air',n_constant=1.,absorption_constant_per_mm=0.),
          'g':Material('glass',n_constant=1.5,absorption_constant_per_mm=.027)}
    original=OpticalSystem((Surface(0.,.02,10.,medium_after='g'),
                            Surface(3.,-.02,10.)),mats,stop_radius=5.,image_z=52.5)
    # Duplicate the front and back boundaries, so both absorbing and air-side
    # zero-distance events occur. No new refractive jump is introduced.
    duplicate=replace(original,surfaces=(original.surfaces[0],original.surfaces[0],
                                        original.surfaces[1],original.surfaces[1]))
    assert validate_geometry(duplicate)['status']=='valid'
    cases=[]
    for field in (Field(),Field(.03,.07)):
        for xy in ((0.,0.),(.7,.1),(-.3,.8)):
            ray=aim_ray(original,field,xy)
            a,b=trace_ray(original,ray),trace_ray(duplicate,ray)
            assert a.success and b.success
            np.testing.assert_allclose(a.positions,b.positions[[0,2]],atol=0.,rtol=0.)
            np.testing.assert_allclose(a.opl,b.opl[[0,2]],atol=0.,rtol=0.)
            for mode in ('unpolarized_jones_power','scalar_interface_average'):
                p,q=(trace_transmission(s,ray,t,mode=mode,missing_absorption='error')
                     for s,t in ((original,a),(duplicate,b)))
                keys=('transmission','absorbed_fraction','reflected_fraction')
                delta=max(abs(p[k]-q[k]) for k in keys)
                assert delta<2e-15,(mode,p,q)
                assert abs(q['ledger_residual'])<2e-15
                cases.append(dict(field=[field.x,field.y],pupil=list(xy),mode=mode,
                                  max_power_delta=delta,transmission=q['transmission']))
    return cases


if __name__=='__main__':
    result={'source':str(ROOT/'src'),'detector_controls':detector_controls(),
            'coincident_power_controls':coincident_power_controls(),'status':'passed'}
    target=Path(__file__).with_name('postfix_physics_c_evidence.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'detector_cases':len(result['detector_controls']),
                      'power_cases':len(result['coincident_power_controls']),
                      'max_pixel_error':max(x['max_pixel_error'] for x in result['detector_controls']),
                      'max_power_delta':max(x['max_power_delta'] for x in result['coincident_power_controls'])}))
