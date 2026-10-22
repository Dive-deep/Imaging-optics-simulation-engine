"""Bulk centered-surface tracing, designed for candidate optimization.

Surface columns: vertex z, curvature, conic, A4,A6,A8,A10,A12,clear radius.
Units are mm; coefficients multiply r**4 through r**12 without normalization.
This module never drops failed rays. Consumers must penalize every failed ray.
STOP/launch geometry and image-plane projection belong to the caller.
"""
from pathlib import Path
import ctypes
import numpy as np

HERE=Path(__file__).resolve().parent
STATUS={0:'ok',1:'escaped_clear_surface',2:'total_internal_reflection',3:'intersection_or_numerical_failure',4:'unsupported_backward_ray'}
_lib=None

def _library():
    global _lib
    if _lib is None:
        path=HERE/'_raytrace_kernel.so'
        if not path.exists():raise RuntimeError('Build native/_raytrace_kernel.so with python native/build.py before tracing')
        _lib=ctypes.CDLL(str(path))
        _lib.trace_rays.argtypes=[ctypes.c_int,ctypes.c_int]+[ctypes.c_void_p]*12
        _lib.trace_rays.restype=None
    return _lib

def trace(surfaces,origins,directions,indices):
    surfaces=np.ascontiguousarray(surfaces,dtype=np.float64)
    origins=np.ascontiguousarray(origins,dtype=np.float64)
    directions=np.ascontiguousarray(directions,dtype=np.float64)
    if surfaces.ndim!=2 or surfaces.shape[1]!=9:raise ValueError('surfaces shape must be (S,9)')
    if origins.ndim!=2 or origins.shape[1]!=3 or directions.shape!=origins.shape:raise ValueError('origins and directions must have shape (N,3)')
    N=len(origins);S=len(surfaces)
    indices=np.asarray(indices,dtype=np.float64)
    if indices.shape==(S+1,):indices=np.broadcast_to(indices,(N,S+1))
    indices=np.ascontiguousarray(indices)
    if indices.shape!=(N,S+1):raise ValueError('indices must have shape (S+1,) or (N,S+1)')
    if not np.all(np.isfinite(surfaces)) or np.any(surfaces[:,8]<=0):raise ValueError('surface parameters must be finite, with positive clear radii')
    if not np.all(np.isfinite(origins)) or not np.all(np.isfinite(directions)):raise ValueError('ray vectors must be finite')
    if not np.all(np.isfinite(indices)) or np.any(indices<=0):raise ValueError('indices must be finite and positive')
    pos=np.empty((N,3));direction=np.empty((N,3));opl=np.empty(N)
    status=np.empty(N,dtype=np.int32);failed=np.empty(N,dtype=np.int32)
    history=np.full((N,S,3),np.nan);dh=np.full((N,S,3),np.nan);oh=np.full((N,S),np.nan)
    args=[surfaces,origins,directions,indices,pos,direction,opl,status,failed,history,dh,oh]
    _library().trace_rays(N,S,*[ctypes.c_void_p(x.ctypes.data) for x in args])
    return {'pos':pos,'dir':direction,'opl':opl,'status':status,'failed_surface':failed,
            'history':history,'direction_history':dh,'opl_history':oh}

def sag(surfaces,r):
    """Convenience reference sag evaluation, broadcast r against surface rows."""
    s=np.asarray(surfaces,float);r=np.asarray(r,float);q=r*r
    c=s[...,1];k=s[...,2]
    with np.errstate(invalid='ignore'):
        z=c*q/(1+np.sqrt(1-(1+k)*c*c*q))
    for j in range(3,8):z=z+s[...,j]*q**(j-1)
    return z
