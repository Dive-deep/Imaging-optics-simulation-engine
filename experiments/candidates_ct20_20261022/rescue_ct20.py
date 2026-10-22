"""Constraint recovery with smooth footprint penalties outside clear apertures.

The expanded aperture is an optimization continuation only. Published results
are traced through the actual apertures and never keep a clipped ray as valid.
"""
from pathlib import Path
import argparse,copy,json,time,traceback
import numpy as np
from scipy.optimize import least_squares
from model import trace_design,native_surfaces,sag,first_order
from optimize import Parameterization,quick_metrics
from refine import summary_proxy,FREQUENCIES
from aperture_fit import fit_apertures


def smooth_merit(d,rings=6,arms=16,na_weight=150.):
    t=trace_design(d,rings=rings,arms=arms,aperture_scale=1.25)
    img=np.nan_to_num(np.clip(t['image_xy'],-100,100),nan=100.)
    ok=(t['status']==0)&t['image_forward'];w=t['weights'];ww=w*ok;total=ww.sum(axis=(0,2))
    center=np.sum(img*ww[...,None],axis=(0,2))/np.maximum(total[:,None],1e-30)
    spot=(img-center[None,:,None,:])*.6*np.sqrt(w)[None,None,:,None]
    targety=6.4*np.sin(np.radians(12)*t['fields'])/np.sin(np.radians(12))
    chief=(img[:,:,0]-np.c_[np.zeros(4),targety][None,:,:])*.4
    height=(img[1,-1,0,1]-6.4)*100.
    phase=2*np.pi*(img-center[None,:,None,:])[...,None]*FREQUENCIES
    otf=np.sum(np.exp(-1j*phase)*ww[...,None,None],axis=(0,2))/np.maximum(total[:,None,None],1e-30)
    rim=t['dir'][:,0,-arms:];rimok=ok[:,0,-arms:]
    na=np.min(np.where(rimok,np.linalg.norm(rim[...,:2],axis=-1),0.),axis=1)
    na_penalty=np.maximum(0,.706-na)*na_weight
    footprint=np.nan_to_num(np.linalg.norm(t['history'][...,:2],axis=-1),nan=100.)
    clipping=np.maximum(0,footprint-np.asarray(d['clear_radii_mm'])+.05)*4/np.sqrt(img.shape[2])
    s=native_surfaces(d);geo=[]
    for j in range(6):
        rr=np.linspace(0,s[j,8],33);ss,dd=sag(s[[j]],rr)
        geo.extend(np.maximum(0,.015-dd[0])*10.)
        geo.extend(np.maximum(0,np.abs(ss[0])-60.)*.5)
    for j in range(5):
        rr=np.linspace(0,min(s[j,8],s[j+1,8]),65);ss,_=sag(s[[j,j+1]],rr)
        sep=s[j+1,0]+ss[1]-s[j,0]-ss[0]
        geo.extend(np.maximum(0,(.5 if j%2==0 else 1.02)-sep)*10.)
    front,_=sag(s[[0]],np.linspace(0,s[0,8],33));back,_=sag(s[[5]],np.linspace(0,s[5,8],33))
    geo.extend([max(0,.15-d['stop_distance_mm']-front.min())*10.,max(0,.15+back.max()-d['axial_distances_mm'][5])*10.])
    geo.extend([max(0,3.03-first_order(d)['paraxial_bfl_mm'])*10.])
    return np.r_[spot.ravel(),chief.ravel(),height,((1-otf.real)*.15).ravel(),
                 (np.minimum(0,otf.real)*.3).ravel(),(~ok).ravel()*5.,
                 na_penalty,clipping.ravel(),geo,np.asarray(d['asphere_sag_coefficients_mm'])[:,2:].ravel()*.0005]


def rescue(path,out,max_nfev=250):
    raw=json.loads(path.read_text());d=copy.deepcopy(raw['design']);d['optimization_ct_upper_mm']=20.
    p=Parameterization(d,orders=5);start=time.perf_counter();calls=0
    def fun(x):
        nonlocal calls
        calls+=1;return smooth_merit(p.unpack(x))
    fit=least_squares(fun,p.initial,bounds=(p.lower,p.upper),max_nfev=max_nfev,
                      x_scale='jac',ftol=1e-7,xtol=1e-8,gtol=1e-7,diff_step=3e-6)
    d=p.unpack(fit.x)
    q=trace_design(d,fields=(1.,),wavelengths=(588.,),xy=np.zeros((1,2)))
    if q['status'][0,0,0]==0:
        v=q['dir'][0,0,0]
        if abs(v[1])>1e-14:
            delta=float((6.4-q['image_xy'][0,0,0,1])*v[2]/v[1])
            if 3.<=d['axial_distances_mm'][5]+delta<=50.:d['axial_distances_mm'][5]+=delta
    d,audit=fit_apertures(d,coupled=False)
    result={'source_result':str(path),'design':d,'ct20_rescue':{
        'method':'expanded-aperture continuation with physical footprint penalty; final actual aperture retrace',
        'optimization_aperture_scale':1.25,'rings':6,'arms':16,'nfev':fit.nfev,
        'residual_calls':calls,'cost':float(fit.cost),'message':fit.message,
        'seconds':time.perf_counter()-start},'aperture_fit':audit,
        'quick_metrics':quick_metrics(d,rings=16,arms=48),'geometric_mtf_proxy':summary_proxy(d)}
    out.mkdir(exist_ok=True,parents=True);(out/path.name).write_text(json.dumps(result,indent=2,allow_nan=False))
    print(path.stem,result['ct20_rescue'],result['quick_metrics']['NA_rgb_axis'],result['quick_metrics']['failed_ray_count'],flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('inputs',nargs='+',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-nfev',type=int,default=250);a=p.parse_args()
    for path in a.inputs:
        if (a.output/path.name).exists():continue
        try:rescue(path,a.output,a.max_nfev)
        except Exception:print(traceback.format_exc(),flush=True)
