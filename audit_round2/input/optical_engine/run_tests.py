"""Run all production-path tests and retain a machine-readable test record."""
import json
from pathlib import Path
import sys
import time
import unittest

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'examples'))

class RecordedResult(unittest.TextTestResult):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs); self.records=[];self.starts={}
    def startTest(self,test):
        self.starts[test.id()]=time.perf_counter();super().startTest(test)
    def _record(self,test,status,detail=None):
        self.records.append({'test':test.id(),'status':status,
             'seconds':time.perf_counter()-self.starts.get(test.id(),time.perf_counter()),'detail':detail})
    def addSuccess(self,test):
        self._record(test,'passed');super().addSuccess(test)
    def addFailure(self,test,err):
        self._record(test,'failed',self._exc_info_to_string(err,test));super().addFailure(test,err)
    def addError(self,test,err):
        self._record(test,'error',self._exc_info_to_string(err,test));super().addError(test,err)
    def addSkip(self,test,reason):
        self._record(test,'skipped',reason);super().addSkip(test,reason)

def write_test_log(report):
    """Human-readable rendering of the actual recorded test results."""
    lines=['Recorded unittest results (rendered from test_results.json)',
           f"Python: {report['python']}",
           f"NumPy: {report['numpy']}; SciPy: {report['scipy']}", '']
    for case in report['cases']:
        lines.append(f"{case['status'].upper()}: {case['test']} ({case['seconds']:.6f}s)")
        if case['detail']:
            lines.append(case['detail'])
    lines.extend(['', f"Tests: {report['tests_run']}; failures: {report['failures']}; "
                  f"errors: {report['errors']}; skipped: {report['skipped']}",
                  f"Elapsed: {report['seconds']:.6f}s; passed: {report['passed']}"])
    (ROOT/'results/test_run.log').write_text('\n'.join(lines)+'\n', encoding='utf-8')

def main():
    start=time.perf_counter()
    suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    result=unittest.TextTestRunner(verbosity=2,resultclass=RecordedResult).run(suite)
    import numpy,scipy
    report={'scope':'Analytic invariants, adversarial cases and independent synthetic specimens; not Zemax/measurement correlation',
      'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,
      'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'skipped':len(result.skipped),
      'passed':result.wasSuccessful(),'seconds':time.perf_counter()-start,'cases':result.records}
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'results/test_results.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    write_test_log(report)
    return 0 if result.wasSuccessful() else 1

if __name__=='__main__':raise SystemExit(main())
