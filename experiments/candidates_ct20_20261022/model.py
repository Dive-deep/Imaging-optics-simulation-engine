"""Explicit high-NA seed model; mm, vacuum-nominal nm, radian angles, float64.
Normal-line dispersion uses prior user nd/Vd and a declared three-term Cauchy fit.
No core engine files or previous numerical evidence are modified.
"""
from pathlib import Path
import sys, json
import numpy as np
E=Path(__file__).resolve().parent
if str(E) not in sys.path:sys.path.insert(0,str(E))
MATERIALS={'PMMA':(1.492,57.47),'PC':(1.586,32.0),'FDS18W':(1.94594,17.98),'TAFD45':(1.95375,32.32)}
WAVELENGTHS=np.array([486.,588.,656.])
ORDERS=np.array([4,6,8,10,12])
NORM_RADIUS=25.

def dispersion_coefficients(name):
    nd,vd=MATERIALS[name];d,F,C,g=np.array([587.5618,486.1327,656.2725,435.8343])/1000
    f=lambda l:np.array([1.,l**-2,l**-4])
    delta=(nd-1)/vd;pgf=.6438-.001682*vd
    return np.linalg.solve(np.array([f(d),f(F)-f(C),f(g)-f(F)]),np.array([nd,delta,pgf*delta]))
COEFFS={name:dispersion_coefficients(name) for name in MATERIALS}
def index_normal_line(name,wavelength_nm):
    if name.lower()=='air':return np.ones_like(np.asarray(wavelength_nm,float))
    l=np.asarray(wavelength_nm,float)/1000;a,b,c=COEFFS[name]
    return a+b/l**2+c/l**4

def native_surfaces(design, *, aperture_scale=1.):
    a=np.zeros((6,9));a[:,0]=np.r_[0.,np.cumsum(design['axial_distances_mm'][:5])]
    a[:,1]=design['curvatures'];a[:,2]=design['conics']
    a[:,3:8]=np.asarray(design['asphere_sag_coefficients_mm'])/NORM_RADIUS**ORDERS
    a[:,8]=np.asarray(design['clear_radii_mm'])*aperture_scale
    for j in range(6):
        if design['materials'][j//2] not in ('PMMA','PC') and (a[j,2]!=0 or np.any(a[j,3:8]!=0)):
            raise ValueError('Glass surfaces must be spherical standard with k=0 and all aspheric coefficients=0')
    return a

def sag(surfaces,r):
    s=np.asarray(surfaces);r=np.asarray(r);r2=r*r
    c=s[...,1,None];k=s[...,2,None]
    disc=1-(1+k)*c*c*r2
    base=c*r2/(1+np.sqrt(np.maximum(disc,1e-30)))
    return base+np.sum(s[...,3:8,None]*r[None,None,:]**ORDERS[None,:,None],axis=-2),disc

def polar_samples(rings=5,arms=12,rim=True):
    # Equal-area midpoint annuli with integration weights; chief and rim zero weights.
    radii=np.sqrt((np.arange(rings)+.5)/rings)
    ang=np.arange(arms)*2*np.pi/arms
    xy=np.stack([np.cos(ang)[None,:]*radii[:,None],np.sin(ang)[None,:]*radii[:,None]],axis=-1).reshape(-1,2)
    xy=np.vstack([np.zeros((1,2)),xy]);weights=np.r_[0.,np.full(rings*arms,1/(rings*arms))]
    if rim:
      xy=np.vstack([xy,np.c_[np.cos(ang),np.sin(ang)]]);weights=np.r_[weights,np.zeros(arms)]
    return xy,weights

def trace_design(design,fields=(0.,.5,.8,1.),wavelengths=WAVELENGTHS,*,rings=5,arms=12,pupil_fraction=1.,aperture_scale=1.,xy=None):
    from native.raytrace import trace
    if xy is None:xy,weights=polar_samples(rings,arms)
    else:xy=np.asarray(xy,float);weights=np.ones(len(xy))/len(xy)
    fields=np.asarray(fields,float);wavelengths=np.asarray(wavelengths,float)
    wf=np.array([(wi,fi) for wi in range(len(wavelengths)) for fi in range(len(fields))]);n=len(xy)
    theta=np.radians(12)*fields[wf[:,1]]
    origins=np.zeros((len(wf),n,3));origins[:,:,:2]=xy*design['stop_radius_mm']*pupil_fraction;origins[:,:,2]=-design['stop_distance_mm']
    directions=np.zeros_like(origins);directions[:,:,1]=np.sin(theta)[:,None];directions[:,:,2]=np.cos(theta)[:,None]
    indices=np.ones((len(wf),n,7))
    for j,name in enumerate(design['materials']):indices[:,:,2*j+1]=index_normal_line(name,wavelengths[wf[:,0]])[:,None]
    out=trace(native_surfaces(design,aperture_scale=aperture_scale),origins.reshape(-1,3),directions.reshape(-1,3),indices.reshape(-1,7))
    for key,val in list(out.items()):
      if isinstance(val,np.ndarray) and len(val)==len(wf)*n:out[key]=val.reshape((len(wavelengths),len(fields),n)+val.shape[1:])
    zimage=float(np.sum(design['axial_distances_mm']))
    dt=(zimage-out['pos'][...,2])/np.where(abs(out['dir'][...,2])>1e-10,out['dir'][...,2],1e-10)
    out['image_xy']=out['pos'][...,:2]+out['dir'][...,:2]*dt[...,None]
    out['image_forward']=dt>=-1e-9
    out['pupil_xy']=xy*design['stop_radius_mm']*pupil_fraction;out['weights']=weights;out['fields']=fields;out['wavelengths_nm']=wavelengths
    return out

def first_order(design,wavelength_nm=588.):
    s=native_surfaces(design);n=np.ones(7)
    for j,name in enumerate(design['materials']):n[2*j+1]=index_normal_line(name,wavelength_nm)
    m=np.eye(2)
    for j in range(6):
      m=np.array([[1.,0.],[-(n[j+1]-n[j])*s[j,1],1.]])@m
      if j<5:m=np.array([[1.,design['axial_distances_mm'][j]/n[j+1]],[0.,1.]])@m
    power=-m[1,0];efl=1/power if abs(power)>1e-15 else np.inf;bfl=-m[0,0]/m[1,0] if abs(m[1,0])>1e-15 else np.inf
    return {'efl_mm':float(efl),'paraxial_bfl_mm':float(bfl),'matrix':m.tolist()}

def geometry_checks(design,points=101):
    s=native_surfaces(design);sags=[];domains=[]
    for row in s:
      rr=np.linspace(0,row[8],points);ss,dd=sag(row[None,:],rr);sags.append(ss[0]);domains.append(float(dd.min()))
    edge=[]
    for j in range(5):
      rr=np.linspace(0,min(s[j,8],s[j+1,8]),points)
      ss,_=sag(s[[j,j+1]],rr);edge.append(float(np.min(s[j+1,0]+ss[1]-s[j,0]-ss[0])))
    return {'minimum_conic_domains':domains,'minimum_adjacent_surface_separation_mm':edge,
      'stop_front_min_clearance_mm':float(design['stop_distance_mm']+np.min(sags[0])),
      'image_back_min_clearance_mm':float(design['axial_distances_mm'][5]-np.max(sags[5])),
      'center_thickness_mm':[design['axial_distances_mm'][j] for j in (0,2,4)],
      'axis_air_gaps_mm':[design['axial_distances_mm'][j] for j in (1,3)]}

def export_prescription(d):
    s=native_surfaces(d);out=[]
    for j,row in enumerate(s):
      mat=d['materials'][j//2];plastic=mat in ('PMMA','PC')
      out.append({'surface':j+1,'lens':j//2+1,'type':'even_asphere' if plastic else 'standard','vertex_z_mm':float(row[0]),
       'radius_mm':None if row[1]==0 else float(1/row[1]),'curvature_per_mm':float(row[1]),
       'thickness_to_next_surface_mm':float(d['axial_distances_mm'][j]),'material_after':mat if j%2==0 else 'air',
       'clear_semi_diameter_mm':float(row[8]),'conic':float(row[2]) if plastic else 0.,
       'aspheric_coefficients':{f'A{order}_mm_pow_{1-order}':float(row[3+i]) for i,order in enumerate(ORDERS)} if plastic else {}})
    return out
