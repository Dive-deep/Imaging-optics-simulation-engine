"""Continue each independent parent into the new CT<=20mm feasible region.

Clipping a starting vector to a new bound is only an optimization initialization.
Every returned design must be independently requalified; no old label is reused.
"""
from pathlib import Path
import argparse, copy, json, time, traceback
import numpy as np
from scipy.optimize import least_squares
from optimize import Parameterization, residual, quick_metrics
from model import trace_design, first_order
from aperture_fit import fit_apertures
from refine import merit, summary_proxy

HERE=Path(__file__).resolve().parent


def refine_cap(d,cap,max_nfev=100,pupil_fraction=1.):
    d=copy.deepcopy(d);d['optimization_ct_upper_mm']=float(cap)
    p=Parameterization(d,orders=5);calls=0;begin=time.perf_counter()
    def fun(x):
        nonlocal calls
        calls+=1;current=p.unpack(x)
        base=residual(current,pupil_fraction=pupil_fraction,fields=(0.,.5,.8,1.),rings=3,arms=8)
        bfl=first_order(current)['paraxial_bfl_mm']
        return np.r_[base,max(0.,3.03-bfl)*5.]
    fit=least_squares(fun,p.initial,bounds=(p.lower,p.upper),max_nfev=max_nfev,
                      x_scale='jac',ftol=2e-6,xtol=2e-6,gtol=2e-6,diff_step=1e-5)
    return p.unpack(fit.x),{'ct_upper_mm':float(cap),'pupil_fraction':pupil_fraction,
        'nfev':fit.nfev,'residual_calls':calls,'cost':float(fit.cost),
        'seconds':time.perf_counter()-begin,'message':fit.message}


def repair(record,output,max_nfev=100,polish_nfev=150):
    d=copy.deepcopy(record['design']);start=time.perf_counter()
    old_ct=np.asarray(d['axial_distances_mm'])[[0,2,4]];maximum=max(old_ct)
    count=max(1,int(np.ceil(max(0,maximum-20)/2.5)))
    caps=np.linspace(maximum,20.,count+1)[1:] if maximum>20 else np.array([20.])
    log=[]
    for cap in caps:
        d,item=refine_cap(d,cap,max_nfev);log.append(item)
        print(json.dumps({'seed':d['seed_id'],'stage':item}),flush=True)
    # Dense-enough merit uses true pupil area and common RGB centroid/OTF.
    d['optimization_ct_upper_mm']=20.
    p=Parameterization(d,orders=5);calls=0;begin=time.perf_counter()
    def fun(x):
        nonlocal calls
        calls+=1
        return merit(p.unpack(x),rings=6,arms=16,mtf_weight=.15)
    fit=least_squares(fun,p.initial,bounds=(p.lower,p.upper),max_nfev=polish_nfev,
                      x_scale='jac',ftol=1e-6,xtol=1e-6,gtol=1e-6,diff_step=5e-6)
    d=p.unpack(fit.x)
    polish={'nfev':fit.nfev,'residual_calls':calls,'cost':float(fit.cost),
            'seconds':time.perf_counter()-begin,'message':fit.message,'rings':6,'arms':16}
    q=trace_design(d,fields=(1.,),wavelengths=(588.,),xy=np.zeros((1,2)))
    delta=None
    if q['status'][0,0,0]==0:
        v=q['dir'][0,0,0]
        if abs(v[1])>1e-14:
            delta=float((6.4-q['image_xy'][0,0,0,1])*v[2]/v[1])
            if 3.<=d['axial_distances_mm'][5]+delta<=50.:d['axial_distances_mm'][5]+=delta
    d,aperture_audit=fit_apertures(d,coupled=False)
    result={'source_result':record['source_result'],'parent_input':record,
            'design':d,'ct20_continuation':log,'ct20_polish':polish,
            'common_image_plane_delta_mm':delta,'aperture_fit':aperture_audit,
            'quick_metrics':quick_metrics(d,rings=16,arms=48),
            'geometric_mtf_proxy':summary_proxy(d),
            'total_seconds':time.perf_counter()-start,
            'qualification':'pending_fresh_CT20_and_full_optical_checks'}
    target=output/(d['seed_id']+'.json');target.write_text(json.dumps(result,indent=2,allow_nan=False))
    print(json.dumps({'seed':d['seed_id'],'DONE':True,'seconds':result['total_seconds'],
                     'CT_mm':[d['axial_distances_mm'][j] for j in (0,2,4)],
                     'NA_rgb':result['quick_metrics']['NA_rgb_axis'],
                     'failed_rays':result['quick_metrics']['failed_ray_count'],
                     'AUC':result['geometric_mtf_proxy']['mean_AUC_normalized']}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-nfev',type=int,default=100);p.add_argument('--polish-nfev',type=int,default=150)
    a=p.parse_args();rows=json.loads(a.input.read_text());rows=rows.get('parents',rows) if isinstance(rows,dict) else rows
    a.output.mkdir(parents=True,exist_ok=True)
    for record in rows:
        target=a.output/(record['design']['seed_id']+'.json')
        if target.exists():continue
        try:repair(record,a.output,a.max_nfev,a.polish_nfev)
        except Exception:
            msg=traceback.format_exc();print(msg,flush=True)
            target.with_suffix('.error.txt').write_text(msg)
