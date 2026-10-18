"""Real sequential rays: actual A2--A12 sag, first forward hit and vector Snell.

Coordinates/lengths are mm. ``intersect`` includes t=0 for coincident dummy
surfaces. A surface is finite inside its mechanical (otherwise clear) cylinder;
its active aperture is applied AFTER selecting its first geometric hit. Powers
are unchanged lens-metric weights here; radiometry handles physical losses.

The bounded interval solver isolates roots in increasing t order. A tangent or
unresolved earlier interval is never skipped to select a later simple root.
This is numerical interval bounding with outward rounding, not a claim of a
formally verified arbitrary-coefficient root solver.
"""
from dataclasses import replace
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import brentq

from .models import Ray, TraceResult, RayBundle
from .settings import EngineSettings
from .errors import InputError, GeometryError, NumericalError, TotalInternalReflection


def _settings(settings):
    return settings if settings is not None else EngineSettings()


def _radial(surface, q, derivatives=False):
    q = np.asarray(q, dtype=float)
    gamma = (1.0 + surface.conic) * surface.curvature**2
    disc = 1.0 - gamma*q
    if np.any(disc < 0):
        raise GeometryError('sag_domain: conic square-root argument is negative')
    root = np.sqrt(disc)
    z = surface.curvature*q/(1.0+root)
    # np.polynomial uses stable nested evaluation with increasing coefficients.
    a = np.asarray((0.0, *surface.coefficients))
    z = z + np.polynomial.polynomial.polyval(q, a)
    if not derivatives:
        return z
    if np.any(root == 0):
        raise GeometryError('grazing_degenerate: conic has infinite slope')
    b = surface.curvature/root
    c = surface.curvature*gamma/root**3
    if len(a) > 1:
        b = b + 2*np.polynomial.polynomial.polyval(q, np.polynomial.polynomial.polyder(a))
    if len(a) > 2:
        c = c + 4*np.polynomial.polynomial.polyval(q, np.polynomial.polynomial.polyder(a, 2))
    return z, b, c


def sag(surface, x, y):
    """Local sag; scalar or broadcast NumPy inputs, stable at c=0 and r=0."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    return _radial(surface, x*x+y*y)


def gradient(surface, x, y):
    """Sag gradient, final dimension (x,y)."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    _, b, _ = _radial(surface, x*x+y*y, True)
    return np.stack((b*x, b*y), axis=-1)


def hessian(surface, x, y):
    """Sag Hessian, final dimensions (2,2); includes A2 in vertex curvature."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    _, b, c = _radial(surface, x*x+y*y, True)
    return np.stack((np.stack((b+c*x*x, c*x*y), axis=-1),
                     np.stack((c*x*y, b+c*y*y), axis=-1)), axis=-2)


def _outward(lo, hi):
    return (float(np.nextafter(lo, -np.inf)), float(np.nextafter(hi, np.inf)))


def _poly_interval(coefficients, qlo, qhi, derivative=0):
    lo = hi = 0.0
    for m, a in enumerate(coefficients, 1):
        if m < derivative or a == 0:
            continue
        a *= math.factorial(m)/math.factorial(m-derivative)
        v0, v1 = a*qlo**(m-derivative), a*qhi**(m-derivative)
        lo += min(v0, v1); hi += max(v0, v1)
    return _outward(lo, hi)


def _sag_interval(surface, qlo, qhi):
    c = surface.curvature
    gamma = (1+surface.conic)*c*c
    d0, d1 = max(0., 1-gamma*qlo), max(0., 1-gamma*qhi)
    v0 = c*qlo/(1+math.sqrt(d0)); v1 = c*qhi/(1+math.sqrt(d1))
    p0, p1 = _poly_interval(surface.coefficients, qlo, qhi)
    return _outward(min(v0, v1)+p0, max(v0, v1)+p1)


def sag_bounds(surface, radius=None):
    """Conservative local z interval on a physical circular surface domain."""
    radius = (surface.mechanical_radius or surface.clear_radius) if radius is None else float(radius)
    if radius <= 0: raise InputError('Sag-bound radius must be positive')
    qmax = radius*radius
    gamma = (1+surface.conic)*surface.curvature**2
    if gamma > 0: qmax = min(qmax, 1/gamma)
    return _sag_interval(surface, 0., qmax)


@lru_cache(maxsize=256)
def _validate_surfaces(surfaces,tol,max_iterations,root_samples):
    for j,s in enumerate(surfaces):
        gamma=(1+s.conic)*s.curvature**2
        if gamma*s.clear_radius**2>=1:
            raise GeometryError(f'sag_domain: surface {j} clear aperture includes non-smooth/undefined conic domain')
    counts=[]
    for j,(front,back) in enumerate(zip(surfaces,surfaces[1:])):
        radius=min(front.clear_radius,back.clear_radius)
        delta=back.z-front.z
        # Subtract equal conics symbolically; otherwise independent monotonic
        # bounds are conservative. Polynomial differences retain cancellation.
        same_conic=(front.curvature==back.curvature and front.conic==back.conic)
        coeff=tuple((back.coefficients[k] if k<len(back.coefficients) else 0.)-
                    (front.coefficients[k] if k<len(front.coefficients) else 0.) for k in range(6))
        fbase=replace(front,coefficients=(),surface_type='standard')
        bbase=replace(back,coefficients=(),surface_type='standard')
        stack=[(0.,radius*radius,0)];visits=0
        while stack:
            a,b,depth=stack.pop();visits+=1
            if visits>root_samples*max_iterations:
                raise NumericalError(f'geometry_unresolved: adjacent pair {j},{j+1} interval budget exhausted')
            pl,ph=_poly_interval(coeff,a,b)
            if same_conic: cl=ch=0.
            else:
                fl,fh=_sag_interval(fbase,a,b);bl,bh=_sag_interval(bbase,a,b)
                cl,ch=bl-fh,bh-fl
            lower,upper=delta+cl+pl,delta+ch+ph
            if lower>=-tol: continue
            midpoint=(a+b)/2
            r=math.sqrt(midpoint)
            gap=delta+float(sag(back,r,0)-sag(front,r,0))
            if upper < -tol or gap < -tol:
                raise GeometryError(f'surface_overlap: adjacent pair {j},{j+1} crosses inside shared clear aperture')
            if depth>=max_iterations or b-a<=tol*tol:
                raise NumericalError(f'geometry_unresolved: adjacent pair {j},{j+1} separation not resolved')
            stack.append((midpoint,b,depth+1));stack.append((a,midpoint,depth+1))
        counts.append(visits)
    return tuple(counts)


def validate_geometry(system,settings=None):
    """Reject conic-domain faults or crossed sequential surfaces over aperture.

    Radial interval bounds cover the entire common CLEAR radius, not just axis
    and rim. Exactly coincident dummy/refracting surfaces remain permitted.
    Results are cached by immutable surface values and numerical tolerances.
    Mechanical rim shapes outside the declared clear surface are not invented.
    """
    settings=_settings(settings)
    counts=_validate_surfaces(tuple(system.surfaces),settings.ray.intersection_tol_mm,
                              settings.ray.max_iterations,settings.ray.root_samples)
    return {'status':'valid','method':'radial interval separation on common clear apertures',
            'pair_interval_counts':counts}


def _q_interval(o, d, a, b):
    pa=o[:2]+a*d[:2]; pb=o[:2]+b*d[:2]
    qa, qb = float(pa@pa), float(pb@pb)
    aa=float(d[:2]@d[:2])
    lo=min(qa,qb)
    if aa:
        t=-float(o[:2]@d[:2])/aa
        if a <= t <= b:
            p=o[:2]+t*d[:2];lo=min(lo,float(p@p))
    return max(0.,np.nextafter(lo,-np.inf)), np.nextafter(max(qa,qb),np.inf)


def _bound_interval(surface, ray, settings):
    o,d=ray.origin,ray.direction;tol=settings.ray.intersection_tol_mm
    radius=surface.mechanical_radius or surface.clear_radius
    gamma=(1+surface.conic)*surface.curvature**2
    if gamma>0: radius=min(radius,math.sqrt(1/gamma))
    a=float(d[:2]@d[:2]); b=2*float(o[:2]@d[:2]); c=float(o[:2]@o[:2])-radius**2
    lo,hi=0.,np.inf
    if a==0:
        if c>tol*max(1.,radius): raise GeometryError('escaped: ray misses physical surface cylinder')
    else:
        disc=b*b-4*a*c
        if disc<0: raise GeometryError('escaped: ray misses physical surface cylinder')
        r=math.sqrt(max(0.,disc))
        lo=max(lo,(-b-r)/(2*a));hi=min(hi,(-b+r)/(2*a))
    zl,zh=sag_bounds(surface,radius);zl+=surface.z;zh+=surface.z
    if abs(d[2])>settings.ray.direction_tol:
        t0,t1=(zl-o[2])/d[2],(zh-o[2])/d[2]
        lo=max(lo,min(t0,t1)-tol);hi=min(hi,max(t0,t1)+tol)
    elif o[2]<zl-tol or o[2]>zh+tol:
        raise GeometryError('escaped: ray misses surface sag range')
    if hi < lo-tol or hi < -tol:
        raise GeometryError('escaped: no forward intersection bounding interval')
    if not np.isfinite(hi):
        raise GeometryError('grazing_degenerate: unbounded coplanar surface intersection')
    return max(0.,lo),max(0.,hi)


def _function_bounds(surface,ray,a,b):
    o,d=ray.origin,ray.direction
    qlo,qhi=_q_interval(o,d,a,b)
    gamma=(1+surface.conic)*surface.curvature**2
    if gamma>0: qhi=min(qhi,1/gamma)
    sl,sh=_sag_interval(surface,qlo,qhi)
    za=o[2]+a*d[2]-surface.z;zb=o[2]+b*d[2]-surface.z
    fl,fh=_outward(min(za,zb)-sh,max(za,zb)-sl)
    r0,r1=math.sqrt(max(0.,1-gamma*qlo)),math.sqrt(max(0.,1-gamma*qhi))
    if min(r0,r1)==0:
        return fl,fh,-np.inf,np.inf
    c0,c1=surface.curvature/(2*r0),surface.curvature/(2*r1)
    pl,ph=_poly_interval(surface.coefficients,qlo,qhi,1)
    gl,gh=min(c0,c1)+pl,max(c0,c1)+ph
    da=2*float((o[:2]+a*d[:2])@d[:2]);db=2*float((o[:2]+b*d[:2])@d[:2])
    products=[gl*da,gl*db,gh*da,gh*db]
    dl,dh=_outward(d[2]-max(products),d[2]-min(products))
    return fl,fh,dl,dh


def intersect(surface, ray, settings=None):
    """Return (physical point, first forward distance).

    Raises GeometryError for escape/unsupported grazing, NumericalError for an
    unresolved first-root interval. No aperture mask is applied by this helper.
    """
    settings=_settings(settings);tol=settings.ray.intersection_tol_mm
    o,d=ray.origin,ray.direction
    lo,hi=_bound_interval(surface,ray,settings)
    def f(t):
        p=o+t*d
        return float(p[2]-surface.z-sag(surface,p[0],p[1]))
    def accepted(t):
        p=o+max(0.,t)*d
        residual=f(max(0.,t))
        if abs(residual)>tol*8:
            raise NumericalError('intersection_unresolved: original sag residual failed')
        g=gradient(surface,p[0],p[1]);normal=np.array([-g[0],-g[1],1.])
        if abs(float(normal@d))/np.linalg.norm(normal)<=settings.ray.direction_tol:
            raise GeometryError('grazing_degenerate: tangential first intersection')
        return p,float(max(0.,t))
    # Analytic plane is exact and also avoids a coplanar continuum of roots.
    if surface.curvature==0 and not any(surface.coefficients):
        if abs(d[2])<=settings.ray.direction_tol:
            raise GeometryError('grazing_degenerate: ray parallel to surface plane')
        t=(surface.z-o[2])/d[2]
        if lo-tol<=t<=hi+tol and t>=-tol: return accepted(t)
        raise GeometryError('escaped: plane intersection lies outside physical cylinder')
    # Generic first-hit root isolation. An explicit stack preserves earliest t.
    stack=[(lo,hi,0)]; visits=0
    budget=settings.ray.root_samples*settings.ray.max_iterations
    while stack:
        a,b,depth=stack.pop();visits+=1
        if visits>budget:
            raise NumericalError('intersection_unresolved: interval work budget exhausted')
        fl,fh,dl,dh=_function_bounds(surface,ray,a,b)
        if fl>tol or fh < -tol: continue
        fa,fb=f(a),f(b)
        monotonic=dl>0 or dh<0
        if monotonic:
            # A small residual alone is NOT a root certificate near a tangent:
            # f=(t-t0)^2 can be tiny well before the double root. Endpoints are
            # accepted only at exact floating zero; all other simple roots need
            # an actual sign-changing bracket.
            if fa==0: return accepted(a)
            if fa*fb<0:
                try:
                    t=brentq(f,a,b,xtol=tol,rtol=4*np.finfo(float).eps,maxiter=settings.ray.max_iterations)
                except (RuntimeError,ValueError) as exc:
                    raise NumericalError('intersection_unresolved: bracket solve failed') from exc
                return accepted(t)
            if fb==0: return accepted(b)
            continue
        if b-a<=tol or depth>=settings.ray.max_iterations:
            # A sign-changing multiple root can also be unresolved here; never
            # classify this interval as empty or skip it to a later root.
            raise GeometryError('grazing_degenerate: non-simple first-root interval')
        m=(a+b)/2
        if m==a or m==b:
            raise NumericalError('intersection_unresolved: floating-point interval exhausted')
        stack.append((m,b,depth+1));stack.append((a,m,depth+1))
    raise GeometryError('escaped: no intersection on physical sag branch')


def refract(direction, normal, n1, n2):
    """Vector Snell; normal orientation is automatically opposed to incidence."""
    d=np.asarray(direction,float);normal=np.asarray(normal,float)
    if d.shape!=(3,) or normal.shape!=(3,) or not np.all(np.isfinite(d)) or not np.all(np.isfinite(normal)):
        raise InputError('Snell vectors must be finite 3-vectors')
    dscale=float(np.max(np.abs(d)));nscale=float(np.max(np.abs(normal)))
    if min(n1,n2)<=0 or not np.isfinite(n1+n2) or dscale==0 or nscale==0:
        raise InputError('Snell indices/vectors must be positive and nonzero')
    # Public callers may supply any finite vector scale, including values whose
    # naive sum-of-squares norm overflows. Normalize scaled vectors first.
    d=d/dscale;d=d/np.linalg.norm(d)
    normal=normal/nscale;normal=normal/np.linalg.norm(normal)
    if np.dot(d,normal)>0: normal=-normal
    cos_i=-float(np.dot(d,normal));eta=n1/n2
    k=1-eta*eta*(1-cos_i*cos_i)
    roundoff=32*np.finfo(float).eps*max(1.,eta*eta)
    if k < -roundoff: raise TotalInternalReflection('tir: no real transmitted Snell branch')
    if abs(k)<=roundoff:
        raise NumericalError('critical_angle_unresolved: transmitted ray is grazing')
    out=eta*d+(eta*cos_i-math.sqrt(k))*normal
    oscale=float(np.max(np.abs(out)))
    if not np.all(np.isfinite(out)) or not np.isfinite(oscale) or oscale==0:
        raise NumericalError('snell_unresolved: nonfinite or degenerate transmitted vector')
    out=out/oscale;out=out/np.linalg.norm(out)
    if np.dot(out,normal)>=0:
        raise NumericalError('snell_unresolved: transmitted hemisphere check failed')
    return out


def trace_ray(system,ray,wavelength_nm=None,settings=None,stop_at=None,ignore_apertures=False):
    """History after each refracted surface; never appends an image-plane hit.

    A failed surface is excluded from history. Missing material data remains an
    input error rather than being mislabeled a blocked physical ray.
    """
    settings=_settings(settings);wavelength=settings.wavelength(wavelength_nm)
    validate_geometry(system,settings)
    last=len(system.surfaces)-1 if stop_at is None else int(stop_at)
    if not 0<=last<len(system.surfaces): raise InputError('stop_at index is invalid')
    p=ray.origin.copy();d=ray.direction.copy();opl=ray.opl
    positions=[];directions=[];opls=[];powers=[]
    def result(status='ok',message=''):
        return TraceResult(np.asarray(positions,float).reshape(-1,3),np.asarray(directions,float).reshape(-1,3),np.asarray(opls),np.asarray(powers),status,message)
    n1=system.index(system.ambient_material,wavelength)
    for j,surface in enumerate(system.surfaces[:last+1]):
        n2=system.index(surface.medium_after,wavelength)
        try:
            hit,distance=intersect(surface,Ray(p,d),settings)
            radius=system.stop_semi_diameter if j==system.stop_index else surface.clear_radius
            if not ignore_apertures and np.hypot(hit[0],hit[1])>radius+settings.ray.intersection_tol_mm:
                return result('aperture_blocked',f'Surface {j}: active aperture blocks first physical hit')
            g=gradient(surface,hit[0],hit[1]);normal=np.array([-g[0],-g[1],1.])
            # An identical medium does not physically refract; its dummy surface
            # still has aperture/OPL significance and must be intersected.
            dout=d.copy() if n1==n2 else refract(d,normal,n1,n2)
            if dout[2]<=settings.ray.direction_tol:
                return result('unsupported_backward_ray',f'Surface {j}: sequential forward model cannot continue')
        except TotalInternalReflection as exc:
            return result('tir',f'Surface {j}: {exc}')
        except (GeometryError,NumericalError) as exc:
            status=str(exc).split(':',1)[0]
            return result(status,f'Surface {j}: {exc}')
        opl+=n1*distance
        positions.append(hit);directions.append(dout);opls.append(opl);powers.append(ray.power)
        p=hit;d=dout;n1=n2
    return result()


def launch_plane_z(system,settings=None):
    """One fixed plane before the entire first finite surface, all fields."""
    settings=_settings(settings)
    lo,_=sag_bounds(system.surfaces[0])
    return system.surfaces[0].z+lo-settings.sampling.boundary_clearance_mm


def _stop_matrix(system,wavelength,launch_z):
    """Reduced-coordinate matrix from launch plane to STOP BEFORE refraction."""
    n=system.index(system.ambient_material,wavelength);z=launch_z;matrix=np.eye(2)
    for j,s in enumerate(system.surfaces[:system.stop_index+1]):
        matrix=np.array([[1.,(s.z-z)/n],[0.,1.]])@matrix
        if j==system.stop_index: break
        n2=system.index(s.medium_after,wavelength)
        cv=s.curvature+2*(s.coefficients[0] if s.coefficients else 0.)
        matrix=np.array([[1.,0.],[-(n2-n)*cv,1.]])@matrix
        n=n2;z=s.z
    return matrix


def aim_ray(system,field,pupil_xy,wavelength_nm=None,settings=None):
    """Solve launch variables so actual STOP sag hit equals normalized target.

    Lens-metric angle rays share a plane-wave eikonal n*d dot origin. Finite
    object rays share zero phase at the object point. Only one locally regular
    aiming branch is supported; a singular Newton map is an explicit failure.
    The stopping condition also verifies the STOP's physical radial predicate.
    This matters at the exact rim when aim_tolerance exceeds intersection_tol:
    a small aiming residual must not turn a valid boundary sample into clipping.
    """
    settings=_settings(settings);wavelength=settings.wavelength(wavelength_nm)
    if field.kind=='angle' and system.object_z is not None:
        raise InputError('Finite object_z requires a height field; finite-object angle conversion is not defined')
    pupil=np.asarray(pupil_xy,float)
    if pupil.shape!=(2,) or not np.all(np.isfinite(pupil)) or np.linalg.norm(pupil)>1+settings.ray.direction_tol:
        raise InputError('Normalized STOP pupil coordinate must lie in unit disk')
    target=pupil*system.stop_semi_diameter
    n0=system.index(system.ambient_material,wavelength)
    stop=system.surfaces[system.stop_index]
    target3=np.array([*target,stop.z+float(sag(stop,*target))])
    if field.kind=='height':
        if system.object_z is None: raise InputError('A height field requires finite object_z')
        origin=np.array([field.x,field.y,system.object_z])
        if target3[2]<=origin[2]: raise InputError('Object must lie before STOP for forward aiming')
        direct=(target3[:2]-origin[:2])/(target3[2]-origin[2])
        matrix=_stop_matrix(system,wavelength,system.object_z)
        if abs(matrix[0,1])<=settings.ray.direction_tol:
            seed=direct
        else: seed=(target-matrix[0,0]*origin[:2])/(matrix[0,1]*n0)
        def build(v): return Ray(origin,np.array([v[0],v[1],1.]),opl=0.)
    else:
        z=launch_plane_z(system,settings)
        direction=np.array([math.tan(field.x),math.tan(field.y),1.]);direction/=np.linalg.norm(direction)
        slope=direction[:2]/direction[2]
        direct=target-(target3[2]-z)*slope
        matrix=_stop_matrix(system,wavelength,z)
        seed=direct if abs(matrix[0,0])<=settings.ray.direction_tol else (target-matrix[0,1]*n0*slope)/matrix[0,0]
        def build(v):
            origin=np.array([v[0],v[1],z])
            return Ray(origin,direction,opl=n0*float(direction@origin))
    def residual(v):
        launch=build(v)
        # STOP aiming constrains incidence position. Refraction AT the STOP is
        # not part of this residual: otherwise a valid incident ray undergoing
        # TIR there would be mislabeled as an aiming solver failure.
        if system.stop_index:
            tr=trace_ray(system,launch,wavelength,settings,stop_at=system.stop_index-1,ignore_apertures=True)
            if not tr.success: raise NumericalError(f'aim_unresolved: geometry-only trace ended {tr.status}')
            incoming=Ray(tr.last_position,tr.last_direction)
        else: incoming=launch
        try: hit,_=intersect(stop,incoming,settings)
        except (GeometryError,NumericalError) as exc:
            raise NumericalError(f'aim_unresolved: STOP incidence {exc}') from exc
        return hit[:2]-target
    seed=np.asarray(seed,float)
    try: r=residual(seed)
    except NumericalError:
        seed=np.asarray(direct,float);r=residual(seed)
    for _ in range(settings.ray.aim_max_iterations):
        actual_stop_xy=target+r
        inside_stop=np.hypot(actual_stop_xy[0],actual_stop_xy[1])<=system.stop_semi_diameter+settings.ray.intersection_tol_mm
        if np.linalg.norm(r,np.inf)<=settings.ray.aim_tolerance_mm and inside_stop:
            return build(seed)
        h=settings.ray.derivative_step
        jac=np.zeros((2,2))
        for k in range(2):
            step=np.zeros(2);step[k]=h
            plus=minus=None
            try: plus=residual(seed+step)
            except NumericalError: pass
            try: minus=residual(seed-step)
            except NumericalError: pass
            if plus is not None and minus is not None: jac[:,k]=(plus-minus)/(2*h)
            elif plus is not None: jac[:,k]=(plus-r)/h
            elif minus is not None: jac[:,k]=(r-minus)/h
            else: raise NumericalError('aim_unresolved: cannot form local aiming derivative')
        if not np.all(np.isfinite(jac)) or np.linalg.cond(jac)>1/settings.ray.direction_tol:
            raise NumericalError('multi_branch_aim: singular or ill-conditioned STOP map')
        step=np.linalg.solve(jac,-r)
        accepted=False
        for backtrack in range(settings.ray.aim_max_iterations):
            candidate=seed+step*(0.5**backtrack)
            try: rr=residual(candidate)
            except NumericalError: continue
            if np.linalg.norm(rr)<np.linalg.norm(r):
                seed,r=candidate,rr;accepted=True;break
        if not accepted: raise NumericalError('aim_unresolved: damped Newton made no progress')
    raise NumericalError('aim_unresolved: aiming iteration budget exhausted')


def trace_bundle(system,field,wavelength_nm=None,settings=None,rings=None,arms=None,include_boundary=None):
    """Gauss-Legendre quadrature in rho², periodic azimuth; weights sum to one.

    This is uniform STOP-area lens-metric sampling, NOT source power, nor an
    assumption of uniform exit-pupil illumination. By default the physical STOP
    rim and chief are additionally stored with ZERO quadrature weight: these
    topology samples do not change geometric moments, but prevent diffraction
    reconstruction from incorrectly using the interior Gauss hull as aperture.
    Count is rings*arms + arms + 1 when boundary samples are enabled.
    """
    settings=_settings(settings);wavelength=settings.wavelength(wavelength_nm)
    rings=settings.sampling.pupil_rings if rings is None else rings
    arms=settings.sampling.pupil_arms if arms is None else arms
    include_boundary=getattr(settings.sampling,'include_pupil_boundary',True) if include_boundary is None else include_boundary
    if not isinstance(include_boundary,(bool,np.bool_)):
        raise InputError('include_boundary must be boolean')
    if int(rings)!=rings or int(arms)!=arms or rings<1 or arms<3:
        raise InputError('Pupil rings >=1 and arms >=3 must be integers')
    nodes,radial_weights=np.polynomial.legendre.leggauss(int(rings))
    points=[];weights=[];traces=[]
    def one(pupil):
        try: ray=aim_ray(system,field,pupil,wavelength,settings)
        except NumericalError as exc: return TraceResult(status=str(exc).split(':',1)[0],message=str(exc))
        return trace_ray(system,ray,wavelength,settings)
    for i,(node,w) in enumerate(zip(nodes,radial_weights)):
        radius=math.sqrt((node+1)/2)
        for a in range(int(arms)):
            theta=2*math.pi*(a+0.5*(i%2))/arms
            pupil=np.array([radius*math.cos(theta),radius*math.sin(theta)])
            points.append(pupil);weights.append(w/(2*arms));traces.append(one(pupil))
    chief=one(np.zeros(2))
    if include_boundary:
        points.append(np.zeros(2));weights.append(0.);traces.append(chief)
        for a in range(int(arms)):
            theta=2*math.pi*a/arms
            pupil=np.array([math.cos(theta),math.sin(theta)])
            points.append(pupil);weights.append(0.);traces.append(one(pupil))
    return RayBundle(tuple(traces),np.asarray(points),np.asarray(weights),field,wavelength,chief)
