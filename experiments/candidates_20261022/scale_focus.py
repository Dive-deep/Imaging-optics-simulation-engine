"""Common-image-plane scan followed by explicit similarity rescaling of a prescription.
Only front-diameter50 designs are eligible: the first clear radius is reset to25mm,
and all rays are retraced to verify that this final aperture did not clip light.
"""
import copy
import numpy as np
from model import trace_design,ORDERS

def normalize_height(d):
 if d['aperture_mode']!='front_diameter50':raise ValueError('Similarity normalization would change fixedEPD50')
 q=trace_design(d,fields=(1.,),wavelengths=(588.,),xy=np.zeros((1,2)))
 if q['status'][0,0,0]!=0:raise ValueError('Chief ray failed')
 h=q['image_xy'][0,0,0,1];scale=6.4/h
 if not np.isfinite(scale) or scale<=0 or not .8<scale<1.2:raise ValueError('Excessive image height normalization')
 new=copy.deepcopy(d)
 new['curvatures']=(np.array(d['curvatures'])/scale).tolist()
 new['asphere_sag_coefficients_mm']=(np.array(d['asphere_sag_coefficients_mm'])/scale**(ORDERS-1)).tolist()
 new['axial_distances_mm']=(np.array(d['axial_distances_mm'])*scale).tolist()
 new['clear_radii_mm']=(np.array(d['clear_radii_mm'])*scale).tolist();new['clear_radii_mm'][0]=25.
 new['stop_distance_mm']*=scale;new['stop_radius_mm']*=scale
 return new,float(scale)
