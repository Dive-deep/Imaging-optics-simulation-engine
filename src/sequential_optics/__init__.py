"""Independent sequential optical engine; not a Zemax numerical clone."""
from .models import Surface, OpticalSystem, Field, Ray, RayBundle, TraceResult
from .settings import EngineSettings, RaySettings, SamplingSettings
from .errors import *
from .materials import Material
from .api import RunRequest, evaluate

__version__='0.2.1'
