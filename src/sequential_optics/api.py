"""Dependency-aware orchestration; low-level modules remain independently reusable."""
from dataclasses import dataclass, field, asdict
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from .models import Field
from .settings import EngineSettings
from .errors import OpticsError, InputError, MissingDataError, ModelApplicabilityError
from .values import FrozenDict

SUPPORTED_METRICS=frozenset({'first_order','spot','mtf','ray_fans','field_curvature','mapping','common_focus'})

@dataclass(frozen=True)
class RunRequest:
    fields: tuple = field(default_factory=lambda:(Field(),))
    wavelengths: tuple | None = None
    metrics: tuple = ('first_order','spot','mtf')
    output_mode: str = 'focal'
    temperatures_c: tuple | None = None
    thermal_options: dict = field(default_factory=dict)
    include_ray_history: bool = True
    fail_fast: bool = False

    def __post_init__(self):
        object.__setattr__(self,'fields',tuple(self.fields))
        object.__setattr__(self,'metrics',tuple(self.metrics))
        if self.wavelengths is not None:object.__setattr__(self,'wavelengths',tuple(self.wavelengths))
        if self.temperatures_c is not None:object.__setattr__(self,'temperatures_c',tuple(self.temperatures_c))
        thermal_options=dict(self.thermal_options)
        allowed_thermal={'lens_groups','anchor_fraction','layout','stop_policy','image_policy','object_policy','validate_geometry','validation_samples','tolerance_mm'}
        if set(thermal_options)-allowed_thermal:
            raise InputError(f'Unknown thermal options: {sorted(set(thermal_options)-allowed_thermal)}')
        if thermal_options.get('lens_groups') is not None:
            from .thermal import LensGroup
            thermal_options['lens_groups']=tuple(LensGroup.from_dict(g) if isinstance(g,dict) else g for g in thermal_options['lens_groups'])
            for group in thermal_options['lens_groups']:
                if not isinstance(group,LensGroup):
                    raise InputError('thermal lens_groups must contain LensGroup objects or dictionaries')
                group.validate()
        object.__setattr__(self,'thermal_options',FrozenDict(thermal_options))
        if not self.fields: raise InputError('At least one field is required')
        if any(not isinstance(f,Field) for f in self.fields):raise InputError('fields must contain Field objects')
        if self.output_mode not in ('focal','afocal'):raise InputError('output_mode must be focal or afocal')
        unknown=set(self.metrics)-SUPPORTED_METRICS
        if unknown:raise InputError(f'Unknown batch metrics {sorted(unknown)}; radiometry/imaging use their explicit APIs')
        if self.wavelengths is not None and not self.wavelengths:raise InputError('Empty wavelength request')
        if self.temperatures_c is not None and not self.temperatures_c:raise InputError('Empty temperature request')

def _attempt(function,fail_fast):
    try:
        value=function()
        return {'status':'ok','analysis_status':value.get('status','computed') if isinstance(value,dict) else 'computed','value':value}
    except OpticsError as exc:
        if fail_fast:raise
        status='missing_data' if isinstance(exc,MissingDataError) else ('not_applicable' if isinstance(exc,ModelApplicabilityError) else 'failed')
        return {'status':status,'error_type':type(exc).__name__,'message':str(exc)}

def _bundle_record(bundle):
    return {'pupil_xy':bundle.pupil_xy,'weights':bundle.weights,'chief':bundle.chief,
            'traces':bundle.traces,'success_count':int(bundle.success_mask.sum()),'total_count':len(bundle.traces)}

def evaluate(system,request=None,settings=None):
    """Evaluate only requested dependencies. Missing spectral/thermal data stay local.

    No source, transmission or CTE input is requested for nominal G spot/FFT.
    `status=ok` means computation returned, not validated vendor equivalence.
    Numerical/applicability flags inside each value remain part of the result.
    """
    from .geometry import trace_bundle
    from .paraxial import paraxial_properties
    from . import metrics as gm
    from .wave import analyze_wave
    from .io import system_to_dict, json_value
    request=request or RunRequest();settings=settings or EngineSettings()
    wavelengths=[settings.wavelength(w) for w in (request.wavelengths if request.wavelengths is not None else (None,))]
    started=time.perf_counter()
    prescription=system_to_dict(system)
    identity=json.dumps(json_value({'prescription':prescription,'settings':asdict(settings),'request':asdict(request)}),sort_keys=True,allow_nan=False)
    out={'schema':'sequential-optics-result/1','engine_version':'0.2.1',
         'input_sha256':hashlib.sha256(identity.encode()).hexdigest(),
         'settings':asdict(settings),'request':asdict(request),'input':prescription,
         'conventions':{'length':'mm','wavelength':'vacuum_nm','index':'absolute',
                        'vendor_equivalence':'untested','model':'real rays + independent scalar diffraction',
                        'inactive_rim_policy':settings.inactive_rim_policy,
                        'inactive_height_effect':'preserved metadata, not a stepped-rim ray intersection model'},
         'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'states':[]}
    temps=request.temperatures_c if request.temperatures_c is not None else (None,)
    for temperature in temps:
        state={'requested_temperature_c':temperature,
               'temperature_c':system.temperature_c if temperature is None else temperature,
               'index_temperature_convention':('supplied_nominal_per_material_reference' if temperature is None and system.temperature_c is None else 'explicit_uniform_temperature'),
               'wavelengths':[]};out['states'].append(state)
        if temperature is None: current=system
        else:
            from .thermal import thermal_system
            transformed=_attempt(lambda:thermal_system(system,temperature,**request.thermal_options),request.fail_fast)
            if transformed['status']!='ok': state.update(transformed);continue
            current=transformed['value']
            # A thermal transformation may return a state + policy ledger.
            if isinstance(current,dict):current=current['system']
            state['thermal_policy']=request.thermal_options
        state['status']='ok'
        for wavelength in wavelengths:
            wr={'wavelength_nm':wavelength,'fields':[]};state['wavelengths'].append(wr)
            if 'first_order' in request.metrics:
                wr['first_order']=_attempt(lambda:paraxial_properties(current,wavelength,settings),request.fail_fast)
            bundles=[]
            need_rays=bool(set(request.metrics)-{'first_order'})
            for field_value in request.fields:
                fr={'field':asdict(field_value)};wr['fields'].append(fr)
                if not need_rays:continue
                br=_attempt(lambda:trace_bundle(current,field_value,wavelength,settings),request.fail_fast)
                if br['status']!='ok':fr['rays']=br;continue
                bundle=br['value'];bundles.append(bundle)
                fr['rays']={'status':'ok','successful':int(bundle.success_mask.sum()),'total':len(bundle.traces)}
                if request.include_ray_history:fr['ray_history']=_bundle_record(bundle)
                if 'spot' in request.metrics:
                    fn=(lambda:gm.spot_diagram(current,bundle,settings=settings)) if request.output_mode=='focal' else (lambda:gm.angular_spot(bundle))
                    fr['spot']=_attempt(fn,request.fail_fast)
                if 'mtf' in request.metrics:
                    angular_frequencies=(np.linspace(0,settings.sampling.max_angular_frequency_cycles_rad,settings.sampling.angular_frequency_samples)
                                         if request.output_mode=='afocal' else None)
                    fr['mtf']=_attempt(lambda:analyze_wave(current,bundle,settings,mode=request.output_mode,frequencies=angular_frequencies),request.fail_fast)
                if 'mapping' in request.metrics:
                    fr['mapping']=_attempt(lambda:gm.chief_mapping(current,[field_value],wavelength_nm=wavelength,angular=request.output_mode=='afocal',settings=settings),request.fail_fast)
                if 'ray_fans' in request.metrics:
                    fr['ray_fans']=_attempt(lambda:gm.ray_fans(current,field_value,wavelength_nm=wavelength,settings=settings,angular=request.output_mode=='afocal'),request.fail_fast)
                if 'field_curvature' in request.metrics:
                    def curvature():
                        reference_z=None
                        if request.output_mode=='afocal':
                            from .geometry import sag_bounds
                            last=current.surfaces[-1]
                            reference_z=last.z+sag_bounds(last)[1]+settings.sampling.boundary_clearance_mm
                        return gm.differential_field_curvature(current,field_value,wavelength_nm=wavelength,settings=settings,reference_z=reference_z)
                    fr['field_curvature']=_attempt(curvature,request.fail_fast)
            if 'common_focus' in request.metrics:
                if request.output_mode=='afocal':wr['common_focus']={'status':'not_applicable','message':'Moving IMA does not alter angular spread; select a physical compensator separately.'}
                elif len(bundles)!=len(request.fields):wr['common_focus']={'status':'failed','message':'At least one requested field has no traced bundle'}
                else:wr['common_focus']=_attempt(lambda:gm.common_geometric_refocus(current,bundles,settings=settings),request.fail_fast)
    out['runtime']['seconds']=time.perf_counter()-started
    return out
