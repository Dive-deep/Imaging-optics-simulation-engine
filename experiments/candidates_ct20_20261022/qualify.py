"""Dense sampled qualification and reversible clear-aperture fitting for starts."""
import json,argparse
from pathlib import Path
from aperture_fit import fit_apertures
from radiometry import evaluate_design
from refine import summary_proxy

def qualify(path,output):
 raw=json.loads(path.read_text());d=raw['design'];d,audit=fit_apertures(d,coupled=False)
 metrics=evaluate_design(d,image_height_tolerance_mm=.001)
 result={'source_result':str(path),'design':d,'aperture_fit':audit,'metrics':metrics,'geometric_mtf_proxy':summary_proxy(d)}
 target=output/path.name;target.write_text(json.dumps(result,indent=2,allow_nan=False))
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('inputs',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 for source in a.inputs:
  if (a.output/source.name).exists():continue
  try:
   r=qualify(source,a.output);print(source.stem,r['metrics']['status'],r['metrics']['failed_constraints'],flush=True)
  except Exception as e:
   import traceback;traceback.print_exc();(a.output/(source.stem+'.error.txt')).write_text(str(e))
