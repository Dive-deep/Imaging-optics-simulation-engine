"""Synthetic, unoptimized specimens. These numbers are NOT material catalog data."""
from dataclasses import replace
from sequential_optics import OpticalSystem, Surface
from sequential_optics.materials import Material

def synthetic_triplet(*,thermal=False,spectral=False,aspheric=False):
    """Three separated lenses; analytic ABCD expected EFL=84.1490162331 mm."""
    mats={'air':Material('ideal_air',n_constant=1.,index_temperature_invariant=True,
                         cte_per_k=0.,absorption_constant_per_mm=0.,absorption_temperature_invariant=True)}
    for name,n,beta,alpha in [('crown',1.50,1e-5,8e-6),('flint',1.62,2e-5,7e-6)]:
        extra={'dn_dT_per_k':beta,'cte_per_k':alpha} if thermal else {}
        wavelengths=(450.,550.,650.) if spectral else (550.,)
        indices=(n+.006,n,n-.004) if spectral else (n,)
        mats[name]=Material('synthetic_'+name,wavelengths_nm=wavelengths,indices=indices,
                            absorption_constant_per_mm=0.,absorption_temperature_invariant=True,**extra)
    raw=[(0.,1/50.,6.,'crown'),(3.,-1/50.,6.,'air'),
         (5.,-1/40.,5.,'flint'),(7.,1/40.,5.,'air'),
         (9.,1/45.,6.,'crown'),(12.,-1/45.,6.,'air')]
    surfaces=tuple(Surface(z=z,curvature=c,clear_radius=a,mechanical_radius=a+.5,
                           medium_after=m,name=f'S{i+1}') for i,(z,c,a,m) in enumerate(raw))
    if aspheric:
        surfaces=(replace(surfaces[0],surface_type='even_asphere',conic=-.5,
                          coefficients=(1e-4,1e-7,-1e-10,1e-13,-1e-16,1e-19)),)+surfaces[1:]
    return OpticalSystem(surfaces,mats,object_z=None,image_z=91.92430937353885,
                         stop_index=2,stop_radius=4.,name='synthetic_unoptimized_triplet')

def synthetic_collimator(**kwargs):
    return replace(synthetic_triplet(**kwargs),object_z=-78.52124279495672,image_z=None,
                   name='synthetic_unoptimized_finite_object_collimator')

def ideal_plane(clear_radius=10., image_z=20., object_z=-20.):
    air=Material('ideal_air',n_constant=1.,index_temperature_invariant=True,
                 cte_per_k=0.,absorption_constant_per_mm=0.,absorption_temperature_invariant=True)
    return OpticalSystem((Surface(0.,0.,clear_radius),),{'air':air},
                         image_z=image_z,object_z=object_z,name='ideal_plane_aperture')
