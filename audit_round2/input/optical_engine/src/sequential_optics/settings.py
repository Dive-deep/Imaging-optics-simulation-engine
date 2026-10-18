"""User-editable run defaults; optical prescriptions carry no global settings."""
from dataclasses import dataclass, field, asdict
from pathlib import Path
import json
import math
from numbers import Integral
from .errors import InputError
from .values import FrozenDict

@dataclass(frozen=True)
class RaySettings:
    intersection_tol_mm: float = 1e-10
    direction_tol: float = 1e-12
    root_samples: int = 128
    max_iterations: int = 80
    aim_tolerance_mm: float = 1e-9
    aim_max_iterations: int = 30
    derivative_step: float = 1e-5

@dataclass(frozen=True)
class SamplingSettings:
    pupil_rings: int = 8
    pupil_arms: int = 24
    include_pupil_boundary: bool = True
    wave_grid_size: int = 64
    fft_padding: int = 4
    boundary_clearance_mm: float = 1e-3
    max_frequency_lp_mm: float = 120.0
    frequency_samples: int = 121
    max_angular_frequency_cycles_rad: float = 10000.0
    angular_frequency_samples: int = 121

@dataclass(frozen=True)
class EngineSettings:
    primary_wavelength: str = 'G'
    report_reference_material: str | None = None
    inactive_rim_policy: str = 'analytic_sag_opaque'
    wavelengths_nm: dict = field(default_factory=lambda: {'R':650.0,'G':550.0,'B':450.0})
    ray: RaySettings = field(default_factory=RaySettings)
    sampling: SamplingSettings = field(default_factory=SamplingSettings)

    def __post_init__(self):
        object.__setattr__(self,'wavelengths_nm',FrozenDict(self.wavelengths_nm))
        if self.inactive_rim_policy!='analytic_sag_opaque':
            raise InputError('Only analytic_sag_opaque inactive rim is implemented; physical stepped rims need a separate profile')
        if self.report_reference_material is not None and (not isinstance(self.report_reference_material,str) or not self.report_reference_material):
            raise InputError('report_reference_material must be a nonempty ID or None=ambient')
        if self.primary_wavelength not in self.wavelengths_nm:
            raise InputError('Primary wavelength label is absent from wavelengths_nm')
        for v in self.wavelengths_nm.values():
            if not math.isfinite(v) or v <= 0:
                raise InputError('Wavelengths must be finite positive vacuum nm')
        for k,v in asdict(self.ray).items():
            if not math.isfinite(v) or v <= 0:
                raise InputError(f'ray.{k} must be positive')
        for k,v in asdict(self.sampling).items():
            if k=='include_pupil_boundary':
                if not isinstance(v,bool):raise InputError('include_pupil_boundary must be boolean')
                continue
            if not math.isfinite(v) or v <= 0:
                raise InputError(f'sampling.{k} must be positive')
        for name, container, keys in [('ray',self.ray,('root_samples','max_iterations','aim_max_iterations')),
                                     ('sampling',self.sampling,('pupil_rings','pupil_arms','wave_grid_size','fft_padding','frequency_samples','angular_frequency_samples'))]:
            for key in keys:
                value=getattr(container,key)
                if isinstance(value,bool) or not isinstance(value,Integral):
                    raise InputError(f'{name}.{key} must be an integer')
        if self.sampling.pupil_arms<3 or self.sampling.wave_grid_size<8 or self.sampling.wave_grid_size%2:
            raise InputError('pupil_arms must be >=3; wave_grid_size must be even and >=8')
        if min(self.sampling.frequency_samples,self.sampling.angular_frequency_samples)<2:
            raise InputError('Frequency sample counts must be >=2')
        if self.ray.direction_tol>=1:raise InputError('direction_tol must be <1')

    def wavelength(self, value=None):
        if value is None: value = self.primary_wavelength
        if isinstance(value,str):
            if value not in self.wavelengths_nm: raise InputError(f'Unknown wavelength label: {value}')
            value = self.wavelengths_nm[value]
        value = float(value)
        if not math.isfinite(value) or value <= 0: raise InputError('Wavelength must be positive vacuum nm')
        return value

    @classmethod
    def from_dict(cls,data):
        d=dict(data)
        if 'ray' in d: d['ray']=RaySettings(**d['ray'])
        if 'sampling' in d: d['sampling']=SamplingSettings(**d['sampling'])
        return cls(**d)

    @classmethod
    def from_json(cls,path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding='utf-8')))

    def to_json(self,path):
        Path(path).write_text(json.dumps(asdict(self),indent=2)+'\n',encoding='utf-8')
