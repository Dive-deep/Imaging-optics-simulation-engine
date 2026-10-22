"""Explicit native candidate-ray kernel build; no fast-math."""
from pathlib import Path
import hashlib,json,subprocess

HERE=Path(__file__).resolve().parent
def build():
    source=HERE/'raytrace.c';target=HERE/'_raytrace_kernel.so'
    command=['gcc','-shared','-fPIC','-O3','-std=c11','-fno-fast-math','-ffp-contract=off',str(source),'-lm','-o',str(target)]
    subprocess.run(command,check=True)
    report={'command':command,'compiler':subprocess.check_output(['gcc','--version'],text=True).splitlines()[0],
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'binary_sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    (HERE/'build.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
if __name__=='__main__':print(json.dumps(build(),indent=2))
