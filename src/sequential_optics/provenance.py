"""Trace identity and integrity guards, without requiring unrelated material data.

Only ``trace_bundle`` seals engine-produced bundles, at their creation time.
Do not recreate a seal at a consuming API: that would bless old rays with a new
prescription. Changing only IMA is supported because tracing ends at the last
physical surface. FFT grid, padding and reporting frequencies are also analysis
choices, not trace inputs. Effective absolute indices are evaluated only at the
bundle's actual wavelength and temperature; no unused RGB/CTE/absorption is read.

Hand-constructed bundles have no engine provenance. Low-level algorithms may
accept them by explicit opt-in for analytic boundary/measurement studies, but
that is never an assertion that those rays solve the supplied prescription.
"""
from dataclasses import asdict, replace
import hashlib
import json
import numpy as np

from .errors import InputError
from .models import TraceProvenance
from .settings import EngineSettings


def _json_hash(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False, default=float).encode()).hexdigest()


def trace_configuration_fingerprint(system, field, wavelength_nm, settings=None):
    """Physical trace inputs, explicitly excluding IMA-only refocus.

    Material identity here means its medium ID and the effective absolute phase
    index used by this monochromatic trace. Adding unrelated dispersion, CTE or
    loss data does not invalidate geometrical rays; physical source amplitudes
    are recomputed by source_wave rather than reused from this trace identity.
    """
    settings = settings or EngineSettings()
    wavelength = settings.wavelength(wavelength_nm)
    media = sorted({system.ambient_material, *(s.medium_after for s in system.surfaces)})
    payload = dict(version='sequential_trace_v1', surfaces=[asdict(s) for s in system.surfaces],
                   object_z=system.object_z, stop_index=int(system.stop_index),
                   stop_radius=system.stop_semi_diameter, ambient_material=system.ambient_material,
                   temperature_c=system.temperature_c,
                   reference_temperature_c=system.reference_temperature_c,
                   material_indices=[(m, system.index(m, wavelength)) for m in media],
                   field=asdict(field), wavelength_nm=wavelength, ray_settings=asdict(settings.ray),
                   boundary_clearance_mm=settings.sampling.boundary_clearance_mm,
                   inactive_rim_policy=settings.inactive_rim_policy)
    return _json_hash(payload)


def bundle_content_fingerprint(bundle):
    """Hash all trace arrays, status, chief, field, wavelength and pupil measure."""
    digest = hashlib.sha256()
    digest.update(_json_hash(dict(field=asdict(bundle.field), wavelength_nm=float(bundle.wavelength_nm))).encode())

    def array(values):
        a = np.ascontiguousarray(values, dtype=np.dtype('<f8'))
        digest.update(repr(a.shape).encode()); digest.update(a.tobytes())

    array(bundle.pupil_xy); array(bundle.weights)
    digest.update(str(len(bundle.traces)).encode())
    for trace in (*bundle.traces, bundle.chief):
        digest.update(str(trace.status).encode())
        for values in (trace.positions, trace.directions, trace.opl, trace.powers):
            array(values)
    return digest.hexdigest()


def _seal_traced_bundle(system, bundle, settings):
    """Private producer hook; never call this to repair a stale input bundle."""
    return replace(bundle, provenance=TraceProvenance(
        trace_configuration_fingerprint(system, bundle.field, bundle.wavelength_nm, settings),
        bundle_content_fingerprint(bundle)))


def validate_bundle_integrity(bundle, *, allow_unverified=False):
    """Reject mutations of an engine bundle; return an explicit provenance status."""
    provenance = bundle.provenance
    if provenance is None:
        if allow_unverified:
            return 'unverified_external_bundle'
        raise InputError('Ray bundle has no trace provenance; retrace the prescription, or explicitly '
                         'opt in to allow_unverified_bundle for a supplied analytic boundary')
    if not isinstance(provenance, TraceProvenance) or provenance.version != 'sequential_trace_v1':
        raise InputError('Unsupported ray bundle trace provenance')
    if provenance.content_fingerprint != bundle_content_fingerprint(bundle):
        raise InputError('Ray bundle trace provenance mismatch: recorded rays, field, wavelength or pupil data changed; retrace')
    return 'verified_engine_trace'


def require_bundle_match(system, bundle, settings=None, *, allow_unverified=False):
    """Check trace-time identity before any state is built from current inputs."""
    status = validate_bundle_integrity(bundle, allow_unverified=allow_unverified)
    if status == 'verified_engine_trace':
        actual = trace_configuration_fingerprint(system, bundle.field, bundle.wavelength_nm, settings)
        if actual != bundle.provenance.configuration_fingerprint:
            raise InputError('Ray bundle trace provenance mismatch: prescription, material index, temperature, '
                             'STOP, physical source point/conjugate or ray settings changed; retrace '
                             '(only IMA motion may reuse the same traced rays)')
    return status
