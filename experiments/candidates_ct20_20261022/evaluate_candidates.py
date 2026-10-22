"""Resume per-candidate scalar wave work against frozen design SHA256 values."""
from pathlib import Path
import argparse, hashlib, json, time
from prepare_candidates import canonical_hash
from wave_evaluation import evaluate_design


def main():
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--worker',type=int,default=0);p.add_argument('--workers',type=int,default=1)
    a=p.parse_args();records=json.loads(a.input.read_text())['candidates'];a.output.mkdir(parents=True,exist_ok=True)
    for r in records[a.worker::a.workers]:
        cid=r['candidate_id'];target=a.output/(cid+'.json');sha=canonical_hash(r['design'])
        if target.exists() and json.loads(target.read_text()).get('design_sha256')==sha:
            print(cid,'already complete',flush=True);continue
        start=time.perf_counter()
        try:
            wave=evaluate_design(r['design'],rings=16,arms=48,quadrature=(128,256),fft_diagnostics=True)
            wave['candidate_id']=cid
            result={'candidate_id':cid,'design_sha256':sha,'wave':wave,
                    'wall_seconds':time.perf_counter()-start,
                    'source_sha256':hashlib.sha256((Path(__file__).parent/'wave_evaluation.py').read_bytes()).hexdigest()}
            temporary=target.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2,allow_nan=False));temporary.replace(target)
            print(cid,round(result['wall_seconds'],2),[x['status'] for x in wave['monochromatic']],flush=True)
        except Exception:
            import traceback
            msg=traceback.format_exc();(a.output/(cid+'.error.txt')).write_text(msg);print(cid,msg,flush=True)


if __name__=='__main__':main()
