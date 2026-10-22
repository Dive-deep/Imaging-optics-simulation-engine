"""Small controls for the new guarded optimizer; no candidate qualification."""
from pathlib import Path
import copy,hashlib,json,time
import numpy as np
from rescue_ct20_guarded import virtual_trace_design,smooth_merit
from native.raytrace import trace
from model import trace_design,native_surfaces

HERE=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    start=time.perf_counter();base=HERE/'rescue_ct20.py';base_before=digest(base)
    d={'curvatures':[.1,.1,.1,0.,.1,.1], 'conics':[0.,0.,0.,2.,-1.,-2.],
       'clear_radii_mm':[9.,10.*(1.-1e-10),10.1,9.,9.,9.], 'unrelated':{'kept':[1,2,3]}}
    before=copy.deepcopy(d);v,meta=virtual_trace_design(d)
    actual=np.asarray(d['clear_radii_mm']);virtual=np.asarray(v['clear_radii_mm'])
    assert d==before and v is not d and v['unrelated'] is not d['unrelated']
    assert np.all(virtual>=actual)
    assert abs(virtual[0]-10.*(1.-1e-8))<1e-14
    assert virtual[1]==actual[1] and virtual[2]==actual[2]
    np.testing.assert_array_equal(virtual[3:],actual[3:]*1.25)
    assert meta['invalid_actual_domain_surface_numbers']==[3]
    assert meta['physical_radius_retained_inside_margin_surface_numbers']==[2]
    s=np.zeros((1,9));s[0,1]=.1;s[0,8]=9.
    origin=np.array([[0.,0.,-1.],[1.,0.,-1.]])
    direction=np.tile([0.,0.,1.],(2,1));indices=np.array([1.,1.5])
    statuses={}
    for name,radius in [('physical_valid',9.),('naive_1p25_invalid',11.25),
                        ('guarded',virtual[0]),('physical_inside_margin',virtual[1]),
                        ('physical_invalid_preserved',virtual[2])]:
        s[0,8]=radius;statuses[name]=trace(s,origin,direction,indices)['status'].tolist()
    assert statuses['physical_valid']==[0,0]
    assert statuses['naive_1p25_invalid']==[3,3]
    assert statuses['guarded']==[0,0]
    assert statuses['physical_inside_margin']==[0,0]
    assert statuses['physical_invalid_preserved']==[3,3]
    sample=HERE/'results/repair/seed_0075.json'
    real=json.loads(sample.read_text())['design'];real_before=copy.deepcopy(real)
    guarded,real_meta=virtual_trace_design(real)
    native_a=native_surfaces(real);native_b=native_surfaces(guarded)
    np.testing.assert_array_equal(native_a[:,:8],native_b[:,:8])
    trace_status={}
    for name,design,scale in [('physical',real,1.),('naive',real,1.25),('guarded',guarded,1.)]:
        t=trace_design(design,rings=2,arms=8,aperture_scale=scale)
        values,counts=np.unique(t['status'],return_counts=True)
        trace_status[name]={str(int(k)):int(v) for k,v in zip(values,counts)}
    assert '3' in trace_status['naive']
    assert '3' not in trace_status['guarded']
    residual=smooth_merit(real,rings=2,arms=8)
    assert residual.ndim==1 and np.all(np.isfinite(residual)) and real==real_before
    assert digest(base)==base_before
    result={'schema':'ct20-guarded-rescue-domain-controls/1','passed':True,
            'guarded_script_sha256':digest(HERE/'rescue_ct20_guarded.py'),
            'base_script_sha256_unchanged':base_before,'control_script_sha256':digest(Path(__file__)),
            'controls':['virtual clear never below physical clear','physical input deep-copy preserved',
                        'valid sphere expanded only to domain margin','valid physical radius inside margin retained',
                        'invalid actual domain preserved and rejected','planes/paraboloids/hyperboloids remain nominal1.25',
                        'native single-surface naive invalidity repaired only for valid physical domain',
                        'actual seed0075 coarse trace eliminates artificial domain failure',
                        'finite guarded residual uses unchanged physical geometry'],
            'synthetic_metadata':meta,'synthetic_native_status':statuses,
            'real_case_source':str(sample.relative_to(HERE)), 'real_case_sha256':digest(sample),
            'real_case_virtual_metadata':real_meta,'real_case_status_counts':trace_status,
            'real_case_merit_residual_length':len(residual),'seconds':time.perf_counter()-start,
            'scope':'Domain-control regression only, not optimization success or 13-field candidate qualification'}
    out=HERE/'results/rescue_guarded_domain_controls.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':True,'output':str(out),'synthetic_statuses':statuses,'real_status_counts':trace_status,'seconds':result['seconds']},indent=2))

if __name__=='__main__':main()
