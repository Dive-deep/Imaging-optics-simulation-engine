"""Portable JSON prescriptions and lossless finite-array result serialization."""
from dataclasses import asdict, is_dataclass
from pathlib import Path
import json
import math
import numpy as np
from .models import OpticalSystem, Surface
from .settings import EngineSettings
from .errors import InputError

SCHEMA_VERSION = 'sequential-optics/1'

def system_to_dict(system):
    return {
        'schema':SCHEMA_VERSION,
        'system':{k:v for k,v in asdict(system).items() if k not in ('surfaces','materials')},
        'surfaces':[asdict(s) for s in system.surfaces],
        'materials':{k:v.to_dict() for k,v in system.materials.items()},
    }

def system_from_dict(data):
    from .materials import Material
    if data.get('schema') != SCHEMA_VERSION:
        raise InputError(f'Expected prescription schema {SCHEMA_VERSION}')
    args=dict(data.get('system',{}))
    if 'surfaces' in args or 'materials' in args:
        raise InputError('Surfaces/materials belong at the document top level')
    return OpticalSystem(surfaces=tuple(Surface(**s) for s in data['surfaces']),
                         materials={k:Material.from_dict(v) for k,v in data['materials'].items()},**args)

def load_system(path):
    return system_from_dict(json.loads(Path(path).read_text(encoding='utf-8')))

def save_system(system,path):
    Path(path).write_text(json.dumps(system_to_dict(system),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def from_sequential_rows(rows,materials,*,first_vertex_z_mm=0.,object_z=None,image_z=None,
                         stop_index=0,stop_radius=None,ambient_material='air',name='',reference_temperature_c=20.):
    """Convert familiar sequential radius/thickness/diameter rows to absolute z.

    Every row is a physical surface (OBJ/IMA are separate). ``thickness_mm``
    is distance to the next vertex; last-row thickness locates IMA if provided.
    On a front lens row that naturally equals lens center thickness; on a rear
    row it equals the next air gap. radius_mm=0 or None denotes a plane. Supply
    radius_mm OR curvature_per_mm, never both. This is not a .ZMX parser.
    """
    if not rows:raise InputError('At least one sequential row is required')
    z=float(first_vertex_z_mm);surfaces=[];derived_image=None
    for i,row in enumerate(rows):
        r=dict(row)
        allowed={'radius_mm','curvature_per_mm','thickness_mm','active_diameter_mm','mechanical_diameter_mm',
                 'conic','coefficients','medium_after','surface_type','name','inactive_height_mm'}
        if set(r)-allowed:raise InputError(f'Unknown sequential row keys: {sorted(set(r)-allowed)}')
        if 'radius_mm' in r and 'curvature_per_mm' in r:raise InputError('Specify radius OR curvature, not both')
        if 'curvature_per_mm' in r:c=float(r['curvature_per_mm'])
        else:
            radius=r.get('radius_mm')
            c=0. if radius is None or radius==0 else 1./float(radius)
        typ={'Standard':'standard','Even Asphere':'even_asphere'}.get(r.get('surface_type'),r.get('surface_type','standard'))
        coefficients=r.get('coefficients',())
        if isinstance(coefficients,dict):
            if set(coefficients)-{f'A{k}' for k in range(2,13,2)}:raise InputError('Only A2,A4,...,A12 supported')
            coefficients=tuple(coefficients.get(f'A{k}',0.) for k in range(2,13,2))
        if 'active_diameter_mm' not in r:raise InputError('Each row needs active_diameter_mm')
        surfaces.append(Surface(z=z,curvature=c,clear_radius=float(r['active_diameter_mm'])/2,
              mechanical_radius=None if r.get('mechanical_diameter_mm') is None else float(r['mechanical_diameter_mm'])/2,
              conic=r.get('conic',0.),coefficients=coefficients,medium_after=r.get('medium_after',ambient_material),
              surface_type=typ,name=r.get('name',f'S{i+1}'),inactive_height_mm=r.get('inactive_height_mm')))
        if i<len(rows)-1 and 'thickness_mm' not in r:raise InputError('Every nonfinal row needs thickness_mm')
        if 'thickness_mm' in r:
            t=float(r['thickness_mm'])
            if not math.isfinite(t) or t<0:raise InputError('Forward sequential thickness must be finite and nonnegative')
            z+=t
            if i==len(rows)-1:derived_image=z
    if image_z is not None and derived_image is not None:
        raise InputError('Specify final thickness OR explicit image_z, not two independent IMA definitions')
    return OpticalSystem(tuple(surfaces),materials,object_z,image_z if image_z is not None else derived_image,
            stop_index=stop_index,stop_radius=stop_radius,ambient_material=ambient_material,
            reference_temperature_c=reference_temperature_c,name=name)

def json_value(value):
    """Complex arrays retain real/imag; invalid numeric values are explicit, not NaN JSON."""
    if isinstance(value,np.ndarray):
        if np.iscomplexobj(value):
            return {'encoding':'complex-array','real':json_value(value.real),'imag':json_value(value.imag)}
        return json_value(value.tolist())
    if isinstance(value,np.generic): return json_value(value.item())
    if is_dataclass(value): return json_value(asdict(value))
    if isinstance(value,dict): return {str(k):json_value(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [json_value(v) for v in value]
    if isinstance(value,complex): return {'encoding':'complex','real':json_value(value.real),'imag':json_value(value.imag)}
    if isinstance(value,float) and not math.isfinite(value):
        return {'encoding':'nonfinite','value':str(value)}
    if isinstance(value,(str,int,float,bool)) or value is None: return value
    raise TypeError(f'Unsupported result type {type(value).__name__}')

def save_results(result,path):
    Path(path).write_text(json.dumps(json_value(result),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
