"""CLI: python -m sequential_optics prescription.json --settings settings.json."""
import argparse
import json
from pathlib import Path
from .api import RunRequest, evaluate
from .models import Field
from .settings import EngineSettings
from .io import load_system, save_results

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prescription')
    p.add_argument('--settings')
    p.add_argument('--request',help='RunRequest JSON, fields are x/y/kind objects')
    p.add_argument('--output',default='result.json')
    a=p.parse_args()
    settings=EngineSettings.from_json(a.settings) if a.settings else EngineSettings()
    req={}
    if a.request:
        req=json.loads(Path(a.request).read_text(encoding='utf-8'))
        if 'fields' in req:req['fields']=tuple(Field(**f) for f in req['fields'])
    result=evaluate(load_system(a.prescription),RunRequest(**req),settings)
    save_results(result,a.output)
    print(a.output)

if __name__=='__main__':main()
