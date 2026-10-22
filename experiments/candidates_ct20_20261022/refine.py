"""Dense physical-pupil refinement with a geometric complex-OTF merit.
The geometric OTF is an optimization surrogate, never relabeled FFT MTF.
"""
from pathlib import Path
import sys,json,copy,time,argparse
import numpy as np
from scipy.optimize import least_squares
from model import trace_design,native_surfaces,sag,first_order,geometry_checks
from optimize import Parameterization,quick_metrics
from aperture_fit import fit_apertures
E=Path(__file__).resolve().parent
FREQUENCIES=np.array([1.3,2.9,4.7,6.8,9.1,12.3,16.2,20.])

def merit(design,rings=6,arms=16,mtf_weight=.15):
 t=trace_design(design,rings=rings,arms=arms)
 img=np.nan_to_num(np.clip(t['image_xy'],-100,100),nan=100.);ok=(t['status']==0)&t['image_forward'];w=t['weights']
 # ONE centroid per field across all RGB; chromatic displacement remains.
 ww=ok*w;total=ww.sum(axis=(0,2));center=np.sum(img*ww[...,None],axis=(0,2))/np.maximum(total[:,None],1e-30)
 residuals=(img-center[None,:,None,:])*.6*np.sqrt(w)[None,None,:,None]
 targety=6.4*np.sin(np.radians(12)*t['fields'])/np.sin(np.radians(12))
 target=np.c_[np.zeros(4),targety]
 chiefs=(img[:,:,0]-target[None,:,:])*.4
 exact_height=(img[1,-1,0,1]-6.4)*80.
 # Uniform incident RGB weights; no wavelength-specific centroid shift.
 phases=2*np.pi*(img-center[None,:,None,:])[...,None]*FREQUENCIES
 complex_otf=np.sum(np.exp(-1j*phases)*ww[...,None,None],axis=(0,2))/np.maximum(total[:,None,None],1e-30)
 otf_res=(1-complex_otf.real)*mtf_weight
 # Avoid exploiting high-frequency magnitude rebounds past an OTF zero.
 negative=np.minimum(0,complex_otf.real)*mtf_weight*2
 rim=t['dir'][:,0,-arms:];rimok=ok[:,0,-arms:]
 na=np.min(np.where(rimok,np.sqrt(np.sum(rim[...,:2]**2,axis=-1)),0.),axis=1)
 napen=np.maximum(0,.710-na)*12.
 s=native_surfaces(design);geo=[]
 for j in range(6):
  rr=np.linspace(0,s[j,8],33);ss,dd=sag(s[[j]],rr);geo.extend(np.maximum(0,.015-dd[0])*10.)
  geo.extend(np.maximum(0,np.abs(ss[0])-60.)*.5)
 for j in range(5):
  rr=np.linspace(0,min(s[j,8],s[j+1,8]),65);ss,_=sag(s[[j,j+1]],rr);sep=s[j+1,0]+ss[1]-s[j,0]-ss[0]
  geo.extend(np.maximum(0,(.5 if j%2==0 else 1.02)-sep)*10.)
 front,_=sag(s[[0]],np.linspace(0,s[0,8],33));back,_=sag(s[[5]],np.linspace(0,s[5,8],33))
 geo.extend([max(0,.15-design['stop_distance_mm']-front.min())*10.,max(0,.15+back.max()-design['axial_distances_mm'][5])*10.])
 fo=first_order(design);geo.extend([max(0,3.02-fo['paraxial_bfl_mm'])*5.,max(0,design['stop_distance_mm']-max(1,fo['efl_mm']))*.1])
 return np.r_[residuals.ravel(),chiefs.ravel(),exact_height,otf_res.ravel(),negative.ravel(),(~ok).ravel()*3.,napen,geo,np.asarray(design['asphere_sag_coefficients_mm'])[:,2:].ravel()*.0005]

def summary_proxy(d,rings=16,arms=48):
 t=trace_design(d,rings=rings,arms=arms);xy=t['image_xy'];valid=(t['status']==0)&t['image_forward'];w=t['weights']*valid;tot=w.sum(axis=(0,2));c=np.sum(xy*w[...,None],axis=(0,2))/np.maximum(tot[:,None],1e-30)
 freq=np.linspace(0,20,81);phase=2*np.pi*(xy-c[None,:,None,:])[...,None]*freq
 otf=(np.exp(-1j*phase)*w[...,None,None]).sum(axis=(0,2))/np.maximum(tot[:,None,None],1e-30)
 return {'method':'geometric_complex_OTF_equal_incident_RGB_common_centroid','frequencies_lp_mm':freq.tolist(),'fields':t['fields'].tolist(),
 'MTF_t':abs(otf[:,1]).tolist(),'MTF_s':abs(otf[:,0]).tolist(),'real_OTF_t':otf[:,1].real.tolist(),'real_OTF_s':otf[:,0].real.tolist(),
 'mean_AUC_normalized':float(np.trapezoid(abs(otf),freq,axis=-1).mean()/20),'minimum_AUC_normalized':float(np.trapezoid(abs(otf),freq,axis=-1).min()/20),
 'MTF20_minimum':float(abs(otf[:,:,-1]).min()),'negative_OTF_points':int(np.sum(otf.real<0)),
 'point_count':int(otf.size),'all_rays_valid':bool(valid.all())}

def refine_result(input_path,output_dir,max_nfev=250,mtf_weight=.15,rings=6,arms=16):
 raw=json.loads(Path(input_path).read_text());original=raw.get('design',raw);d,audit=fit_apertures(original,coupled=False)
 p=Parameterization(d,orders=5);calls=[0];start=time.perf_counter()
 def fun(x):calls[0]+=1;return merit(p.unpack(x),rings=rings,arms=arms,mtf_weight=mtf_weight)
 fit=least_squares(fun,p.initial,bounds=(p.lower,p.upper),max_nfev=max_nfev,x_scale='jac',ftol=1e-6,xtol=1e-6,gtol=1e-6,diff_step=5e-6)
 d=p.unpack(fit.x)
 # Final common image plane: solve the green chief full-field height exactly.
 q=trace_design(d,fields=(1.,),wavelengths=(588.,),xy=np.array([[0.,0.]]));v=q['dir'][0,0,0];delta=(6.4-q['image_xy'][0,0,0,1])*v[2]/v[1]
 if q['status'][0,0,0]==0 and d['axial_distances_mm'][5]+delta>=3.:d['axial_distances_mm'][5]+=float(delta)
 d,final_ap=fit_apertures(d,coupled=False)
 report={'source_result':str(input_path),'initial_design':original,'initial_aperture_fit':audit,'design':d,'final_aperture_fit':final_ap,
 'refinement':{'merit':'common RGB geometric complex OTF + area weighted RMS + rim/NA + geometry; not FFT optimization','frequencies_lp_mm':FREQUENCIES.tolist(),'rings':rings,'arms':arms,'mtf_weight':mtf_weight,'max_nfev':max_nfev,'nfev':fit.nfev,'residual_calls':calls[0],'cost':float(fit.cost),'seconds':time.perf_counter()-start,'message':fit.message,'image_plane_delta_mm':float(delta)},'quick_metrics':quick_metrics(d,rings=16,arms=48),'geometric_mtf_proxy':summary_proxy(d)}
 out=Path(output_dir);out.mkdir(parents=True,exist_ok=True);name=d['seed_id'];f=out/(name+'.json');f.write_text(json.dumps(report,indent=2,allow_nan=False))
 print(json.dumps({'seed':name,'seconds':report['refinement']['seconds'],'cost':report['refinement']['cost'],'proxyAUC':report['geometric_mtf_proxy']['mean_AUC_normalized'],'minMTF20':report['geometric_mtf_proxy']['MTF20_minimum'],'rmsG':report['quick_metrics']['rms_mm'][1],'NA':report['quick_metrics']['NA_rgb_axis'],'all_valid':report['geometric_mtf_proxy']['all_rays_valid']}),flush=True)
 return report

def main():
 ap=argparse.ArgumentParser();ap.add_argument('inputs',type=Path,nargs='+');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--max-nfev',type=int,default=250);ap.add_argument('--mtf-weight',type=float,default=.15);ap.add_argument('--rings',type=int,default=6);ap.add_argument('--arms',type=int,default=16);a=ap.parse_args()
 for f in a.inputs:
  if (a.output/f.name).exists():continue
  try:refine_result(f,a.output,a.max_nfev,a.mtf_weight,a.rings,a.arms)
  except Exception as e:
   import traceback;traceback.print_exc();a.output.mkdir(parents=True,exist_ok=True);(a.output/(f.stem+'.error.txt')).write_text(str(e))
if __name__=='__main__':main()
