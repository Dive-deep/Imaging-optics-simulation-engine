"""First-order quantities only; real rays are computed in geometry.py.

The reduced state is (height, n_absolute*angle), so every matrix determinant is
one. Vertex power includes A2. Report EFL uses the explicit settings reference
material (default: system ambient). Selecting an air-reference material is a
report normalization, not a validated proprietary adapter.
"""
from dataclasses import replace
import numpy as np
from .settings import EngineSettings
from .errors import InputError, NumericalError
from .geometry import _stop_matrix


def paraxial_properties(system,wavelength_nm=None,settings=None):
    settings=settings or EngineSettings();w=settings.wavelength(wavelength_nm)
    n0=system.index(system.ambient_material,w);n=n0;matrix=np.eye(2)
    z=system.surfaces[0].z
    for surface in system.surfaces:
        matrix=np.array([[1.,(surface.z-z)/n],[0.,1.]])@matrix
        n2=system.index(surface.medium_after,w)
        cv=surface.curvature+2*(surface.coefficients[0] if surface.coefficients else 0.)
        matrix=np.array([[1.,0.],[-(n2-n)*cv,1.]])@matrix
        n=n2;z=surface.z
    A,B,C,D=matrix.ravel()
    reference=settings.report_reference_material or system.ambient_material
    if reference not in system.materials:
        raise InputError(f'Report-reference material {reference!r} is absent from the prescription')
    na=system.index(reference,w)
    # This declared conditioning threshold has dimensions 1/mm. It is a
    # numerical diagnostic, not an optical claim that tiny powers equal zero.
    threshold=settings.ray.direction_tol/max(1.,abs(system.surfaces[-1].z-system.surfaces[0].z))
    afocal=abs(C)<=threshold
    efl=None if afocal else -n/C
    report=None if afocal else -na/C
    bfl=None if afocal else -n*A/C
    image_z=None if bfl is None else z+bfl
    magnification=None
    finite_matrix=None
    if system.object_z is not None:
        finite_matrix=matrix@np.array([[1.,(system.surfaces[0].z-system.object_z)/n0],[0.,1.]])
        ao,bo,co,do=finite_matrix.ravel()
        if abs(do)<=settings.ray.direction_tol:
            image_z=None
        else:
            distance=-n*bo/do;image_z=z+distance;magnification=ao+distance*co/n
    return dict(matrix=matrix,determinant=float(np.linalg.det(matrix)),
                efl_image_mm=efl,efl_report_mm=report,bfl_mm=bfl,
                bfl_from_last_surface_mm=bfl,paraxial_image_z=image_z,
                magnification=magnification,finite_object_matrix=finite_matrix,
                n_object=n0,n_image=n,report_reference_material=reference,
                report_reference_index=na,wavelength_nm=w,
                status='afocal_or_ill_conditioned' if afocal else 'ok',
                power_abs_per_mm=-C,conditioning_threshold_per_mm=threshold,
                convention='absolute-index reduced ray; BFL from last physical surface')


def resolve_aperture(system,epd_mm=None,f_number=None,settings=None):
    """Resolve nominal infinite-conjugate EPD/F-number into physical STOP.

    F-number means abs(report EFL)/EPD, not finite-conjugate working F/#.
    STOP clear radius is a fixed physical bound and is never silently enlarged.
    """
    settings=settings or EngineSettings()
    if epd_mm is not None and f_number is not None:
        raise InputError('Choose EPD or nominal F-number, not both')
    if epd_mm is None and f_number is None: return system
    w=settings.wavelength()
    if f_number is not None:
        if not np.isfinite(f_number) or f_number<=0: raise InputError('F-number must be positive')
        efl=paraxial_properties(system,w,settings)['efl_report_mm']
        if efl is None: raise NumericalError('A finite nominal F-number requires finite EFL')
        epd_mm=abs(efl)/f_number
    if not np.isfinite(epd_mm) or epd_mm<=0: raise InputError('EPD must be positive')
    # Parallel on-axis ray at first vertex y_EP reaches STOP with A*y_EP.
    prefix=_stop_matrix(system,w,system.surfaces[0].z)
    scale=abs(prefix[0,0])
    if scale<=settings.ray.direction_tol:
        raise NumericalError('EPD conversion is singular: STOP at paraxial focus')
    radius=scale*epd_mm/2
    if radius>system.surfaces[system.stop_index].clear_radius+settings.ray.intersection_tol_mm:
        raise InputError('Resolved physical STOP exceeds its fixed clear aperture')
    return replace(system,stop_radius=min(radius,system.surfaces[system.stop_index].clear_radius))
