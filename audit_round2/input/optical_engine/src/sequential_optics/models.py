"""Shared value types. Lens prescriptions use mm, radians and absolute n."""
from dataclasses import dataclass, field
from typing import Any
import math
from numbers import Integral
import numpy as np
from .errors import InputError
from .values import FrozenDict

@dataclass(frozen=True)
class Surface:
    z: float
    curvature: float
    clear_radius: float
    conic: float = 0.
    coefficients: tuple = ()
    mechanical_radius: float | None = None
    medium_after: str = 'air'
    surface_type: str = 'standard'
    name: str = ''
    inactive_height_mm: float | None = None

    def __post_init__(self):
        if self.surface_type not in ('standard','even_asphere'):
            raise InputError('Only standard and even_asphere surfaces are supported')
        if not all(math.isfinite(v) for v in (self.z,self.curvature,self.clear_radius,self.conic,*self.coefficients)):
            raise InputError('Surface numbers must be finite; represent a plane with curvature=0')
        if self.clear_radius <= 0: raise InputError('clear_radius must be positive')
        if len(self.coefficients)>6: raise InputError('Only A2 through A12 supported')
        if self.surface_type=='standard' and any(self.coefficients):
            raise InputError('Standard surface cannot have asphere coefficients')
        if self.mechanical_radius is not None and (not math.isfinite(self.mechanical_radius) or self.mechanical_radius<self.clear_radius):
            raise InputError('Mechanical radius must be finite and >= clear radius')
        if self.inactive_height_mm is not None and (not math.isfinite(self.inactive_height_mm) or self.inactive_height_mm<0):
            raise InputError('Inactive-region height metadata must be finite and nonnegative')
        object.__setattr__(self,'coefficients',tuple(self.coefficients))

@dataclass(frozen=True)
class OpticalSystem:
    surfaces: tuple
    materials: dict
    object_z: float | None = None
    image_z: float | None = None
    stop_index: int = 0
    stop_radius: float | None = None
    ambient_material: str = 'air'
    reference_temperature_c: float = 20.
    temperature_c: float | None = None
    name: str = ''

    def __post_init__(self):
        object.__setattr__(self,'surfaces',tuple(self.surfaces))
        object.__setattr__(self,'materials',FrozenDict(self.materials))
        if not self.surfaces: raise InputError('At least one surface is required')
        if isinstance(self.stop_index,bool) or not isinstance(self.stop_index,Integral):raise InputError('STOP index must be an integer')
        if not 0 <= self.stop_index < len(self.surfaces): raise InputError('Invalid STOP index')
        if any(b.z<a.z for a,b in zip(self.surfaces,self.surfaces[1:])):
            raise InputError('Forward sequential vertices must be nondecreasing in z')
        for key in [self.ambient_material]+[s.medium_after for s in self.surfaces]:
            if key not in self.materials: raise InputError(f'Missing material {key}')
        if self.stop_radius is not None and (not math.isfinite(self.stop_radius) or self.stop_radius<=0 or self.stop_radius>self.surfaces[self.stop_index].clear_radius):
            raise InputError('STOP radius must be positive and inside its clear aperture')
        for x in (self.object_z,self.image_z,self.temperature_c,self.reference_temperature_c):
            if x is not None and not math.isfinite(x): raise InputError('Positions and temperatures must be finite')

    @property
    def stop_semi_diameter(self):
        return self.stop_radius if self.stop_radius is not None else self.surfaces[self.stop_index].clear_radius

    def index(self,material_id,wavelength_nm):
        return self.materials[material_id].index(wavelength_nm,self.temperature_c)

@dataclass(frozen=True)
class Field:
    x: float = 0.
    y: float = 0.
    kind: str = 'angle'

    def __post_init__(self):
        if self.kind not in ('angle','height'): raise InputError('Field kind must be angle or height')
        if not math.isfinite(self.x) or not math.isfinite(self.y): raise InputError('Fields must be finite')
        if self.kind=='angle' and (abs(self.x)>=np.pi/2 or abs(self.y)>=np.pi/2):
            raise InputError('Forward field angles must be inside (-pi/2,pi/2)')

@dataclass(frozen=True)
class Ray:
    origin: np.ndarray
    direction: np.ndarray
    opl: float = 0.
    power: float = 1.

    def __post_init__(self):
        o=np.asarray(self.origin,dtype=float);d=np.asarray(self.direction,dtype=float)
        if o.shape!=(3,) or d.shape!=(3,) or not np.all(np.isfinite(o)) or not np.all(np.isfinite(d)) or np.max(np.abs(d))==0:
            raise InputError('Ray origin/direction must be finite 3-vectors with nonzero direction')
        scaled=d/np.max(np.abs(d));unit=scaled/np.linalg.norm(scaled)
        o=o.copy();o.setflags(write=False);unit.setflags(write=False)
        object.__setattr__(self,'origin',o);object.__setattr__(self,'direction',unit)
        if not np.isfinite(self.opl) or not np.isfinite(self.power) or self.power<0:
            raise InputError('Ray OPL must be finite and power nonnegative')

@dataclass
class TraceResult:
    positions: np.ndarray = field(default_factory=lambda: np.empty((0,3)))
    directions: np.ndarray = field(default_factory=lambda: np.empty((0,3)))
    opl: np.ndarray = field(default_factory=lambda: np.empty(0))
    powers: np.ndarray = field(default_factory=lambda: np.empty(0))
    status: str = 'ok'
    message: str = ''

    @property
    def success(self): return self.status=='ok'
    @property
    def last_position(self): return self.positions[-1]
    @property
    def last_direction(self): return self.directions[-1]
    @property
    def last_opl(self): return float(self.opl[-1])
    @property
    def last_power(self): return float(self.powers[-1])

@dataclass
class RayBundle:
    traces: tuple
    pupil_xy: np.ndarray
    weights: np.ndarray
    field: Field
    wavelength_nm: float
    chief: TraceResult

    @property
    def success_mask(self): return np.array([t.success for t in self.traces],dtype=bool)
