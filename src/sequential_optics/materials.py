"""Supplied material data with explicit wavelength/index/temperature conventions.

The ray-engine boundary always receives absolute phase index at vacuum nm.
This module does not infer dispersion or thermo-optic data from Abbe number,
CTE, or a generic material name. See calculation_spec_ko.md, part 2.
"""
from dataclasses import asdict, dataclass
import math
from typing import Any

import numpy as np

from .errors import InputError, MissingDataError


def _finite(value, label):
    if not math.isfinite(float(value)):
        raise InputError(f'{label} must be finite')
    return float(value)


def _axis(axis, values, label, positive_axis=False):
    x = np.asarray(axis, dtype=float)
    y = np.asarray(values, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) == 0:
        raise InputError(f'{label}: matching nonempty one-dimensional axes required')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise InputError(f'{label}: data must be finite')
    if np.any(np.diff(x) <= 0) or (positive_axis and np.any(x <= 0)):
        raise InputError(f'{label}: axis must be strictly increasing' + (' and positive' if positive_axis else ''))
    return x, y


def _interpolate(x, axis, values, label):
    """Piecewise linear, no extrapolation, including single measured samples."""
    a, b = np.asarray(axis, float), np.asarray(values, float)
    # Equality tolerance handles harmless floating point basis conversions only.
    tol = 8 * np.finfo(float).eps * max(1., abs(x), abs(a[0]), abs(a[-1]))
    if x < a[0] - tol or x > a[-1] + tol:
        raise MissingDataError(f'{label}: requested {x:g} outside supplied [{a[0]:g}, {a[-1]:g}]')
    return float(np.interp(np.clip(x, a[0], a[-1]), a, b))


@dataclass(frozen=True)
class Material:
    """Declarative supplied material model; fields round-trip through JSON.

    Choose exactly one index model: ``n_constant`` (explicitly achromatic),
    ``wavelengths_nm/indices`` (including a single known wavelength), or
    ``sellmeier_b/sellmeier_c_um2``. ``temperature_indices`` optionally replaces
    the reference index table by a rectangular [temperature,wavelength] table.
    ``dn_dT_per_k`` is a supplied ABSOLUTE-index derivative; wavelength-table
    derivatives may instead use ``dn_dT_values_per_k`` on ``wavelengths_nm``.

    ``index_basis`` = absolute or reference_air. ``wavelength_basis`` = vacuum
    or reference_air, independently. The reference-air conversion requires a
    supplied ``reference_air_index``; it is a fixed reference-condition value,
    not an atmospheric model. A relative n(T) table is intentionally rejected:
    convert every temperature row to absolute first using its own air data.

    CTE models: ``linear_reference`` uses L(T)=1+alpha*(T-T_ref);
    ``instantaneous`` integrates alpha=(1/L)dL/dT exponentially;
    ``reference_derivative`` integrates alpha_ref=dL/dT / L(T_ref).
    The length ratio always uses the geometry baseline, which can differ from
    the material reference temperature. All CTE values are 1/K, not ppm/K.

    Omitted absorption is unknown, never implicitly zero. Loss-free ray and
    lens-MTF analyses need not call absorption at all.
    """
    name: str
    n_constant: float | None = None
    wavelengths_nm: tuple = ()
    indices: tuple = ()
    sellmeier_b: tuple = ()
    sellmeier_c_um2: tuple = ()
    reference_temperature_c: float = 20.
    dn_dT_per_k: float | None = None
    dn_dT_values_per_k: tuple = ()
    temperatures_c: tuple = ()
    temperature_indices: tuple = ()
    index_temperature_invariant: bool = False
    index_basis: str = 'absolute'
    wavelength_basis: str = 'vacuum'
    reference_air_index: float | None = None
    wavelength_range_nm: tuple | None = None
    temperature_range_c: tuple | None = None
    cte_per_k: float | None = None
    cte_model: str = 'linear_reference'
    cte_reference_temperature_c: float | None = None
    cte_temperatures_c: tuple = ()
    cte_values_per_k: tuple = ()
    absorption_constant_per_mm: float | None = None
    absorption_wavelengths_nm: tuple = ()
    absorption_values_per_mm: tuple = ()
    absorption_temperature_invariant: bool = False
    abbe_number: float | None = None
    source: str = ''
    model_id: str = 'supplied_data'

    def __post_init__(self):
        for key in ('wavelengths_nm', 'indices', 'sellmeier_b', 'sellmeier_c_um2',
                    'dn_dT_values_per_k', 'temperatures_c', 'cte_temperatures_c',
                    'cte_values_per_k', 'absorption_wavelengths_nm', 'absorption_values_per_mm'):
            object.__setattr__(self, key, tuple(getattr(self, key)))
        object.__setattr__(self, 'temperature_indices', tuple(tuple(r) for r in self.temperature_indices))
        models = int(self.n_constant is not None) + bool(self.indices or self.temperature_indices) + bool(self.sellmeier_b)
        if models != 1:
            raise InputError('Choose exactly one supplied index model')
        _finite(self.reference_temperature_c, 'reference temperature')
        if self.n_constant is not None and _finite(self.n_constant, 'index') <= 0:
            raise InputError('Phase index must be positive')
        if self.index_basis not in ('absolute', 'reference_air') or self.wavelength_basis not in ('vacuum', 'reference_air'):
            raise InputError('Unsupported index_basis or wavelength_basis')
        if self.index_basis == 'reference_air' or self.wavelength_basis == 'reference_air':
            if self.reference_air_index is None:
                raise MissingDataError('Reference-air convention requires supplied reference_air_index')
            if _finite(self.reference_air_index, 'reference_air_index') <= 0:
                raise InputError('reference_air_index must be positive')
        if self.indices:
            _, y = _axis(self.wavelengths_nm, self.indices, 'index table', True)
            if np.any(y <= 0): raise InputError('Phase indices must be positive')
        if self.sellmeier_b:
            if len(self.sellmeier_b) != len(self.sellmeier_c_um2) or not all(math.isfinite(v) for v in (*self.sellmeier_b, *self.sellmeier_c_um2)):
                raise InputError('Sellmeier B,C pairs must be finite and equal length')
        elif self.sellmeier_c_um2:
            raise InputError('Sellmeier C without B')
        if self.temperature_indices:
            if self.index_basis != 'absolute' or self.wavelength_basis != 'vacuum':
                raise InputError('Convert temperature tables to absolute n and vacuum wavelengths before import')
            _axis(self.temperatures_c, range(len(self.temperatures_c)), 'temperature axis')
            _axis(self.wavelengths_nm, range(len(self.wavelengths_nm)), 'wavelength axis', True)
            a = np.asarray(self.temperature_indices, float)
            if a.shape != (len(self.temperatures_c), len(self.wavelengths_nm)) or not np.all(np.isfinite(a)) or np.any(a <= 0):
                raise InputError('temperature_indices must be positive rectangular [T,wavelength] data')
        elif self.temperatures_c:
            raise InputError('Temperature axis without index values')
        thermals = bool(self.temperature_indices) + (self.dn_dT_per_k is not None) + bool(self.dn_dT_values_per_k) + bool(self.index_temperature_invariant)
        if thermals > 1: raise InputError('Choose one index temperature model')
        if self.dn_dT_per_k is not None: _finite(self.dn_dT_per_k, 'dn/dT')
        if self.dn_dT_values_per_k:
            _axis(self.wavelengths_nm, self.dn_dT_values_per_k, 'dn/dT table', True)
        if self.cte_model not in ('linear_reference', 'instantaneous', 'reference_derivative'):
            raise InputError('Unknown CTE definition')
        if self.cte_per_k is not None: _finite(self.cte_per_k, 'CTE')
        if self.cte_reference_temperature_c is not None: _finite(self.cte_reference_temperature_c, 'CTE reference')
        if self.cte_values_per_k:
            _axis(self.cte_temperatures_c, self.cte_values_per_k, 'CTE table')
            if self.cte_per_k is not None or self.cte_model == 'linear_reference':
                raise InputError('Tabular CTE requires instantaneous or reference_derivative model, without scalar CTE')
        elif self.cte_temperatures_c: raise InputError('CTE axis without values')
        if self.absorption_constant_per_mm is not None:
            if _finite(self.absorption_constant_per_mm, 'absorption') < 0: raise InputError('Absorption cannot be negative')
            if self.absorption_values_per_mm: raise InputError('Choose scalar or table absorption')
        if self.absorption_values_per_mm:
            _, y = _axis(self.absorption_wavelengths_nm, self.absorption_values_per_mm, 'absorption table', True)
            if np.any(y < 0): raise InputError('Absorption cannot be negative')
        elif self.absorption_wavelengths_nm: raise InputError('Absorption wavelengths without values')
        for attr in ('wavelength_range_nm', 'temperature_range_c'):
            r = getattr(self, attr)
            if r is not None:
                if len(r) != 2 or not all(math.isfinite(v) for v in r) or r[1] < r[0]:
                    raise InputError(f'{attr} must contain ordered finite bounds')
                object.__setattr__(self, attr, tuple(r))

    def _wavelength(self, wavelength_nm):
        w = _finite(wavelength_nm, 'vacuum wavelength_nm')
        if w <= 0: raise InputError('Wavelength must be positive')
        if self.wavelength_range_nm is not None and not self.wavelength_range_nm[0] <= w <= self.wavelength_range_nm[1]:
            raise MissingDataError(f'{self.name}: vacuum wavelength outside valid model range')
        return w

    def _temperature(self, temperature_c):
        t = self.reference_temperature_c if temperature_c is None else _finite(temperature_c, 'temperature')
        if self.temperature_range_c is not None and not self.temperature_range_c[0] <= t <= self.temperature_range_c[1]:
            raise MissingDataError(f'{self.name}: temperature outside valid model range')
        return t

    def index(self, wavelength_nm, temperature_c=None):
        """Absolute phase index; missing unrelated CTE/transmission is irrelevant."""
        w, t = self._wavelength(wavelength_nm), self._temperature(temperature_c)
        argument = w / self.reference_air_index if self.wavelength_basis == 'reference_air' else w
        if self.temperature_indices:
            values = [_interpolate(w, self.wavelengths_nm, row, f'{self.name} wavelength') for row in self.temperature_indices]
            return _interpolate(t, self.temperatures_c, values, f'{self.name} temperature')
        if self.n_constant is not None:
            n = float(self.n_constant)
        elif self.indices:
            n = _interpolate(argument, self.wavelengths_nm, self.indices, f'{self.name} wavelength')
        else:
            l2 = (argument * 1e-3) ** 2
            den = l2 - np.asarray(self.sellmeier_c_um2)
            if np.any(den == 0): raise InputError('Sellmeier pole at requested wavelength')
            n2 = 1. + float(np.sum(np.asarray(self.sellmeier_b) * l2 / den))
            if not math.isfinite(n2) or n2 <= 0: raise InputError('Sellmeier index is not positive real')
            n = math.sqrt(n2)
        if self.index_basis == 'reference_air': n *= self.reference_air_index
        dt = t - self.reference_temperature_c
        if dt != 0:
            if self.dn_dT_per_k is not None: n += self.dn_dT_per_k * dt
            elif self.dn_dT_values_per_k:
                # Derivatives are absolute dn/dT; their wavelength axis has the
                # same declared wavelength convention as the dispersion table.
                n += _interpolate(argument, self.wavelengths_nm, self.dn_dT_values_per_k, 'dn/dT wavelength') * dt
            elif not self.index_temperature_invariant:
                raise MissingDataError(f'{self.name}: non-reference temperature requires n(T) or explicit index_temperature_invariant')
        if not math.isfinite(n) or n <= 0: raise InputError('Computed index must be positive finite')
        return float(n)

    def absorption(self, wavelength_nm, temperature_c=None):
        """Internal POWER attenuation coefficient [1/mm], distinct from Fresnel."""
        w, t = self._wavelength(wavelength_nm), self._temperature(temperature_c)
        if self.absorption_constant_per_mm is None and not self.absorption_values_per_mm:
            raise MissingDataError(f'{self.name}: absorption/transmission data not supplied')
        if t != self.reference_temperature_c and not self.absorption_temperature_invariant:
            raise MissingDataError(f'{self.name}: thermal absorption unavailable; explicitly select absorption_temperature_invariant if intended')
        if self.absorption_constant_per_mm is not None: return float(self.absorption_constant_per_mm)
        # Absorption tables are always explicitly vacuum wavelength tables.
        return _interpolate(w, self.absorption_wavelengths_nm, self.absorption_values_per_mm, 'absorption wavelength')

    def _cte_integral(self, start, end):
        if self.cte_per_k is not None: return self.cte_per_k * (end - start)
        if not self.cte_values_per_k: raise MissingDataError(f'{self.name}: thermal geometry requires supplied CTE')
        lo, hi = sorted((start, end))
        a = self.cte_temperatures_c
        low = _interpolate(lo, a, self.cte_values_per_k, 'CTE temperature')
        high = _interpolate(hi, a, self.cte_values_per_k, 'CTE temperature')
        inner = [x for x in a if lo < x < hi]
        xs = np.array([lo, *inner, hi])
        ys = np.array([low, *[_interpolate(x, a, self.cte_values_per_k, 'CTE temperature') for x in inner], high])
        # Trapezoids exactly integrate the chosen piecewise-linear CTE model.
        integral = float(np.sum((ys[1:] + ys[:-1]) * np.diff(xs) / 2))
        return integral if end >= start else -integral

    def expansion_scale(self, temperature_c, baseline_temperature_c=None):
        """Length ratio L(T)/L(T0), independent of thermo-optic requirements."""
        t = _finite(temperature_c, 'temperature')
        t0 = self.reference_temperature_c if baseline_temperature_c is None else _finite(baseline_temperature_c, 'baseline temperature')
        if t == t0: return 1.
        # temperature_range_c is the supplied material model's common validity
        # envelope. Geometry-only operation must not silently extrapolate it.
        self._temperature(t)
        self._temperature(t0)
        ref = self.reference_temperature_c if self.cte_reference_temperature_c is None else self.cte_reference_temperature_c
        if self.cte_model == 'instantaneous':
            try: scale = math.exp(self._cte_integral(t0, t))
            except OverflowError as exc: raise InputError('Expansion scale overflow') from exc
        else:
            l0 = 1. + self._cte_integral(ref, t0)
            lt = 1. + self._cte_integral(ref, t)
            if l0 <= 0 or lt <= 0: raise InputError('CTE length model gives nonpositive length')
            scale = lt / l0
        if not math.isfinite(scale) or scale <= 0: raise InputError('Expansion scale must be positive finite')
        return scale

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        return cls(**data)

    @classmethod
    def from_internal_transmission(cls, *, name, wavelength_nm, transmission, thickness_mm, **kwargs):
        """Create supplied absorption table from INTERNAL transmission.

        Inputs transmission are fractions (0,1], not external transmittance.
        Interpolation is linear in mu=-ln(tau)/thickness, not in tau.
        Index data must still be provided in kwargs.
        """
        d = _finite(thickness_mm, 'transmission sample thickness')
        w, tau = _axis(wavelength_nm, transmission, 'internal transmission', True)
        if d <= 0 or np.any(tau <= 0) or np.any(tau > 1):
            raise InputError('Internal transmission requires d>0 and 0<tau<=1')
        return cls(name=name, absorption_wavelengths_nm=tuple(w), absorption_values_per_mm=tuple(-np.log(tau) / d), **kwargs)
