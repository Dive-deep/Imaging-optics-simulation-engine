"""Retraced small STOP-radius scan at fixed first-lens diameter and fixed CT.

No lens is rescaled and no pupil ray count is reduced. A smaller STOP is allowed
only in the user-authorized first-lens-opening50 alternative, with actual NA
and throughput freshly tested. EPD50 prescriptions are never changed here.
"""
from pathlib import Path
import copy,json,sys,time
import numpy as np
from model import trace_design
from qualify import qualify
from prepare_candidates import eligible,score

HERE=Path(__file__).resolve().parent


def scan(paths):
    output=HERE/'results/stop_scan';output.mkdir(exist_ok=True)
    qualified=HERE/'results/qualified_stop_scan';qualified.mkdir(exist_ok=True)
    temporary=HERE/'results/stop_scan_trials';temporary.mkdir(exist_ok=True)
    logs=[]
    for source in paths:
        raw=json.loads(source.read_text());d0=raw['design'];sid=d0['seed_id']
        if d0['aperture_mode']!='front_diameter50':continue
        candidates=[];trials=[];begin=time.perf_counter()
        offsets=np.r_[0.,np.arange(-.02,-.601,-.02),np.arange(.02,.601,.02)]
        for delta in offsets:
            d=copy.deepcopy(d0);d['stop_radius_mm']+=float(delta)
            # Fast full-circle rim and all4fields check; final acceptance is13fields.
            t=trace_design(d,rings=16,arms=48)
            ok=(t['status']==0)&t['image_forward']
            na=np.min(np.linalg.norm(t['dir'][:,0,-48:,:2],axis=-1),axis=1)
            entry={'delta_radius_mm':float(delta),'all4field_rays_valid':bool(ok.all()),'rim_NA_min':float(na.min())}
            trials.append(entry)
            if not ok.all() or na.min()<.7:continue
            record={'source_result':str(source),'design':d,'ct20_stop_scan':{
                'initial_radius_mm':d0['stop_radius_mm'],'delta_radius_mm':float(delta),
                'first_lens_clear_diameter_mm':50.,'fixed_CT_and_surface_prescription':True,
                'scope':'small STOP scan within first-opening50 alternative; no per-field/per-color refocus'}}
            temp=temporary/(sid+'.json');temp.write_text(json.dumps(record,indent=2))
            q=qualify(temp,temporary/'qualified') if (temporary/'qualified').exists() else None
            if q is None:
                (temporary/'qualified').mkdir();q=qualify(temp,temporary/'qualified')
            entry['all_constraints_passed']=eligible(q)
            if eligible(q):candidates.append((score(q),record,q))
        if candidates:
            _,record,q=max(candidates,key=lambda row:row[0]);dest=output/(sid+'.json')
            # Keep the best recovered checkpoint of each distinct parent.
            if (qualified/dest.name).exists() and score(json.loads((qualified/dest.name).read_text()))>=score(q):
                chosen='existing_better'
            else:
                dest.write_text(json.dumps(record,indent=2));q['source_result']=str(dest)
                (qualified/dest.name).write_text(json.dumps(q,indent=2,allow_nan=False));chosen=record['ct20_stop_scan']
        else:chosen=None
        logs.append({'source':str(source),'seed':sid,'trials':trials,'chosen':chosen,'seconds':time.perf_counter()-begin})
        print(sid,len(candidates),'valid STOP samples',chosen,flush=True)
        (output/'scan_log.json').write_text(json.dumps(logs,indent=2))


if __name__=='__main__':scan([Path(x) for x in sys.argv[1:]])
