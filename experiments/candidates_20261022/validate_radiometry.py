"""Independent analytic and original-engine checks for candidate radiometry."""
from pathlib import Path
import copy, json, math, sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'src'))
from model import trace_design, native_surfaces, index_normal_line, first_order
from radiometry import evaluate_design, fresnel_transmission, _moments
from native.validate import reference
from native.raytrace import trace


def main():
    records=[]
    def check(name,passed,**data):
        records.append({'name':name,'passed':bool(passed),**data})
    # Independent closed-form six plane-interface power transmission.
    plane={'seed_id':'analytic_plane_control','materials':['PMMA','PC','TAFD45'],
           'curvatures':[0.]*6,'conics':[0.]*6,
           'asphere_sag_coefficients_mm':[[0.]*5 for _ in range(6)],
           'clear_radii_mm':[50.]*6,'axial_distances_mm':[3.,1.,3.,1.,3.,3.],
           'stop_radius_mm':20.,'stop_distance_mm':1.,'aperture_mode':'control'}
    t=trace_design(plane,rings=4,arms=12)
    bare,evidence=fresnel_transmission(plane,t)
    theory=[]
    for wavelength in t['wavelengths_nm']:
        expected=1.
        for material in plane['materials']:
            n=float(index_normal_line(material,wavelength))
            expected*=(1-((n-1)/(n+1))**2)**2
        theory.append(expected)
    error=float(np.max(np.abs(bare[:,0]-np.asarray(theory)[:,None])))
    check('plane_six_interface_unpolarized_transmission',error<1e-14,max_absolute_error=error)
    check('plane_snell_invariant',evidence['maximum_surface_snell_sine_residual']<1e-7,
          maximum_snell_residual=evidence['maximum_surface_snell_sine_residual'])
    # Weighted statistic has unequal spectral weights and a rejected bright ray.
    xy=np.array([[0.,0.],[2.,0.],[100.,100.]])
    m=_moments(xy,np.array([.2,.3,.5]),np.array([True,True,False]))
    check('weighted_conditional_moment_preserves_lost_mass',
          abs(m['transmitted_quadrature_weight']-.5)<1e-15 and
          np.max(np.abs(m['centroid_xy_mm']-[1.2,0.]))<1e-15 and
          abs(m['rms_radius_mm']-math.sqrt(.96))<1e-15)
    # Material normal-line anchors are actual nd and Vd, not 550nm indices.
    from model import MATERIALS
    for material,(nd,vd) in MATERIALS.items():
        calc_nd=float(index_normal_line(material,587.5618))
        nF=float(index_normal_line(material,486.1327));nC=float(index_normal_line(material,656.2725))
        calc_v=(calc_nd-1)/(nF-nC)
        check('material_'+material+'_nd_Vd',abs(calc_nd-nd)<5e-15 and abs(calc_v-vd)<2e-11,
              nd_error=calc_nd-nd,Vd_error=calc_v-vd)
    # Compare actual optimized high-NA design to independent original Python
    # intersections, gradient and Snell code, including chief and full-rim rays.
    path=HERE/'results/pilot/seed_0083.json'
    design=json.loads(path.read_text())['design']
    surfaces=native_surfaces(design)
    points=np.array([[0.,0.],[0.,.5],[0.,1.],[0.,-1.],[1.,0.],[-1.,0.]])
    max_position=max_direction=max_opl=max_image=0.;mismatches=[];cases=0
    for wavelength in (486.,588.,656.):
        ns=np.ones(7)
        for j,material in enumerate(design['materials']):ns[2*j+1]=index_normal_line(material,wavelength)
        for fraction in (0.,.5,.8,1.):
            angle=math.radians(12)*fraction;direction=np.array([0.,math.sin(angle),math.cos(angle)])
            for point in points:
                cases+=1;origin=np.r_[point*design['stop_radius_mm'],-design['stop_distance_mm']]
                got=trace(surfaces,origin[None],direction[None],ns)
                ref=reference(surfaces,origin,direction,ns)
                if (int(got['status'][0]),int(got['failed_surface'][0]))!=ref[:2]:
                    mismatches.append({'wavelength':wavelength,'field':fraction,'pupil':point.tolist()});continue
                count=len(ref[2])
                if count:
                    max_position=max(max_position,float(np.max(np.abs(got['history'][0,:count]-ref[2]))))
                    max_direction=max(max_direction,float(np.max(np.abs(got['direction_history'][0,:count]-ref[3]))))
                    max_opl=max(max_opl,float(np.max(np.abs(got['opl_history'][0,:count]-ref[4]))))
                if ref[0]==0:
                    z=sum(design['axial_distances_mm'])
                    a=got['pos'][0,:2]+(z-got['pos'][0,2])*got['dir'][0,:2]/got['dir'][0,2]
                    b=ref[2][-1,:2]+(z-ref[2][-1,2])*ref[3][-1,:2]/ref[3][-1,2]
                    max_image=max(max_image,float(np.max(np.abs(a-b))))
    check('high_NA_actual_design_original_python_reference',not mismatches and max_position<1e-8 and max_direction<1e-9 and max_opl<1e-8 and max_image<1e-8,
          ray_cases=cases,status_mismatches=mismatches,max_position_difference_mm=max_position,
          max_direction_difference=max_direction,max_opl_difference_mm=max_opl,max_image_difference_mm=max_image)
    # Independent first-order check from one infinitesimal ray through original code.
    ns=np.ones(7)
    for j,material in enumerate(design['materials']):ns[2*j+1]=index_normal_line(material,588.)
    h=1e-5
    ref=reference(surfaces,np.array([0.,h,-design['stop_distance_mm']]),np.array([0.,0.,1.]),ns)
    slope=ref[3][-1,1]/ref[3][-1,2]
    efl=-h/slope;bfl=ref[2][-1,2]-surfaces[-1,0]-ref[2][-1,1]/slope
    fo=first_order(design)
    check('ABCD_against_infinitesimal_original_real_ray',abs(fo['efl_mm']-efl)<2e-6 and abs(fo['paraxial_bfl_mm']-bfl)<2e-6,
          efl_difference_mm=fo['efl_mm']-efl,paraxial_bfl_difference_mm=fo['paraxial_bfl_mm']-bfl)
    # High-NA complete-rim symmetry and geometric RI's radiance/Jacobian law.
    audit=evaluate_design(design)
    green=audit['by_wavelength_nm']['588.0']['fields']
    na=green['0.0']['image_space_NA']
    check('complete_rim_NA_is_not_one_ray_maximum',na['complete_rim_transmitted'] and na['valid_rim_ray_count']==48 and abs(na['valid_rim_maximum']-na['valid_rim_minimum'])<1e-9,
          rim_min=na['valid_rim_minimum'],rim_max=na['valid_rim_maximum'])
    check('onaxis_RI_normalized_both_models',green['0.0']['RI_geometric']==1. and green['0.0']['RI_uncoated_fresnel']==1.)
    check('all_native_failures_retained',sum(audit['all_requested_native_status_counts'].values())==12*817)
    # Independently evaluate r derivative at a third, smaller step to audit RI.
    h=2.5e-6;theta=math.radians(12)
    q=trace_design(design,[(theta-h)/theta,(theta+h)/theta],(588.,),xy=[[0.,0.]])
    dr=(q['image_xy'][0,1,0,1]-q['image_xy'][0,0,0,1])/(2*h)
    f0=audit['chief_mapping_evidence']['axis_chief_slope_mm_per_rad'][1]
    ratio=green['1.0']['geometric_transmitted_pupil_fraction']/green['0.0']['geometric_transmitted_pupil_fraction']
    ri=ratio*f0*f0*math.cos(theta)*math.sin(theta)/(green['1.0']['chief_height_mm']*dr)
    check('RI_independent_smaller_chief_derivative_step',abs(ri-green['1.0']['RI_geometric'])<2e-7,
          independent_RI=ri,reported_RI=green['1.0']['RI_geometric'])
    # Invalid rim is never relabelled successful and does not certify nominal EPD.
    clipped=copy.deepcopy(design);clipped['clear_radii_mm'][4]=1.
    clipped_audit=evaluate_design(clipped,rings=4,arms=12)
    check('vignetting_has_no_NA_or_EPD_false_certificate',clipped_audit['aperture']['EPD_mm'] is None and
          not clipped_audit['constraint_results']['onaxis_588nm_NA_at_least0p7']['passed'] and
          clipped_audit['all_invalid_ray_count_including_chief_rim']>0)
    plane_audit=evaluate_design(plane,rings=4,arms=12)
    check('afocal_plane_stack_does_not_pass_infinite_BFL',
          not plane_audit['constraint_results']['positive_finite_paraxial_EFL']['passed'] and
          not plane_audit['constraint_results']['paraxial_BFL_at_least3mm']['passed'])
    out={'schema':'independent-radiometry-controls/1','passed':all(r['passed'] for r in records),
         'checks':len(records),'records':records,'source_design':str(path.relative_to(ROOT)),
         'scope':'analyticplaneFresnel/materialanchors/weightedloss,originalPythonactualhighNA72rays,firstorder,RIderivative,NAandvignetting'}
    (HERE/'results/radiometry_controls.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2));assert out['passed']


if __name__=='__main__':main()
