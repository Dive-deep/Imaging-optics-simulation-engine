"""Explicit errors: missing data, invalid physics, unsupported model, numerical failure."""

class OpticsError(Exception):
    pass

class InputError(OpticsError, ValueError):
    pass

class MissingDataError(InputError):
    pass

class GeometryError(OpticsError):
    pass

class NumericalError(OpticsError):
    pass

class ModelApplicabilityError(OpticsError):
    pass

class TotalInternalReflection(OpticsError):
    pass
