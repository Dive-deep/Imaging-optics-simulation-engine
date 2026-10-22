"""Exercise the public evaluator at exact CT limits; no optical winner claim."""
from pathlib import Path
import copy, hashlib, json
import numpy as np

from radiometry import evaluate_design


def main():
    root=Path(__file__).resolve().parent
    # Plane stack deliberately fails focusing/NA. These controls validate only
    # the explicitly reported CT gates, not that this is a usable imaging lens.
    plane={'seed_id':'ct20_boundary_control','materials':['PMMA','PC','TAFD45'],
           'curvatures':[0.]*6,'conics':[0.]*6,
           'asphere_sag_coefficients_mm':[[0.]*5 for _ in range(6)],
           'clear_radii_mm':[50.]*6,'axial_distances_mm':[10.,1.,10.,1.,10.,3.],
           'stop_radius_mm':20.,'stop_distance_mm':1.,'aperture_mode':'control'}
    records=[]
    for lens,slot in enumerate((0,2,4),1):
        for name,value,upper_expected,lower_expected in (
            ('exact20',20.,True,True),
            ('nextafter20_up',float(np.nextafter(20.,np.inf)),False,True),
            ('nextafter20_down',float(np.nextafter(20.,-np.inf)),True,True),
            ('exact2',2.,True,False),
            ('nextafter2_up',float(np.nextafter(2.,np.inf)),True,True)):
            design=copy.deepcopy(plane);design['axial_distances_mm'][slot]=value
            # No optimizer bound or homotopy metadata can weaken final gates.
            design['optimization_ct_upper_mm']=35.
            result=evaluate_design(design,rings=2,arms=8)
            upper=result['constraint_results']['center_thickness_at_most20mm']['passed']
            lower=result['constraint_results']['center_thickness_strictly_over2mm']['passed']
            passed=upper is upper_expected and lower is lower_expected
            records.append({'lens':lens,'case':name,'center_thickness_mm':value,
                            'hex_value':value.hex(),'upper_gate_passed':upper,
                            'upper_expected':upper_expected,'lower_gate_passed':lower,
                            'lower_expected':lower_expected,'passed':passed,
                            'whole_design_qualification':result['status']})
    report={'schema':'ct20-exact-boundary-controls/1','passed':all(r['passed'] for r in records),
            'checks':len(records),'records':records,
            'evaluator_sha256':hashlib.sha256((root/'radiometry.py').read_bytes()).hexdigest(),
            'selection_auditor_sha256':hashlib.sha256((root/'validate_selection_physics.py').read_bytes()).hexdigest(),
            'scope':'15public-evaluator CT gate checks; exact20 accepted,nextrepresentablefloatabove20 rejected for each lens; exact2 rejected,nextfloatabove2 accepted; optimizerhomotopy35doesnotchangefinalgate. Planecontrolnotanopticalcandidate.'}
    out=root/'results/ct20_boundary_controls.json'
    with out.open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(report,indent=2));assert report['passed']


if __name__=='__main__':main()
