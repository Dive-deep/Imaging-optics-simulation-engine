"""Package CT20 candidates and reproducible source; keep the prior delivery intact."""
from pathlib import Path
import hashlib
import json
import platform
import zipfile

import numpy
import scipy

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NAME = 'optics_initial_candidates_ct20_23_20261022'
OUTPUT = REPO.parent / 'outputs'
MANIFEST = HERE / 'results/final_artifact_manifest.json'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    data = json.loads((HERE / 'candidates.json').read_text())
    physics = json.loads((HERE / 'results/final_selection_physics_audit.json').read_text())
    native = json.loads((HERE / 'results/native_final_candidates/summary.json').read_text())
    assert len(data['candidates']) == 23
    assert physics['passed'] and physics['candidate_count'] == 23
    assert native['passed'] and native['candidate_count'] == 23
    assert data['summary']['all_distinct_eligible_parents_retained']
    for candidate in data['candidates']:
        assert all(2 < candidate['design']['axial_distances_mm'][j] <= 20 for j in (0, 2, 4))
    for name in ('candidate_report.html', 'candidate_layouts_and_mtf.pdf', 'candidate_contact_sheet.png'):
        assert (HERE / 'reports' / name).is_file()
    files = set()
    for folder in ('src', 'tests', 'examples', 'docs',
                   'experiments/candidates_20261022', 'experiments/candidates_ct20_20261022'):
        for path in (REPO / folder).rglob('*'):
            if not path.is_file() or '__pycache__' in path.parts or path.suffix in ('.pyc', '.so', '.tmp'):
                continue
            rel = path.relative_to(REPO).as_posix()
            if '/reports/plots/' in rel or rel.startswith('experiments/candidates_20261022/reports/'):
                continue
            if path == MANIFEST:
                continue
            files.add(path)
    for name in ('AGENTS.md', 'README.md', 'CODEX_START.md', 'HANDOFF.md',
                 'FIX_C_REPORT.md', 'FIX9_REPORT.md', 'VALIDATION.md',
                 'pyproject.toml', 'requirements-validation.txt', 'run_tests.py', 'LICENSE'):
        if (REPO / name).is_file():
            files.add(REPO / name)
    manifest = {
        'schema': 'ct20-candidate-delivery/1', 'candidate_count': 23,
        'all_eligible_distinct_parents_retained': True,
        'center_thickness_mm': {'exclusive_minimum': 2, 'inclusive_maximum': 20},
        'baseline_commit': 'e347d56b238f86cedd239921335245f85163b806',
        'python': platform.python_version(), 'numpy': numpy.__version__, 'scipy': scipy.__version__,
        'physics_dense_rays': physics['fresh_dense_ray_count'],
        'physics_failed_rays': physics['fresh_failed_ray_count'],
        'native_original_compared_rays': native['total_rays'],
        'native_original_passed': native['passed'],
        'wave_status_counts': data['summary']['primary_wave_status_counts'],
        'omissions': ['compiled native library and Python cache (rebuild locally)',
                      'duplicate individual report PNGs (embedded in new HTML and PDF)',
                      'prior campaign rendered reports (prior separate archive preserved)'],
        'files': {str(p.relative_to(REPO)): {'bytes': p.stat().st_size, 'sha256': sha(p.read_bytes())}
                  for p in sorted(files)},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    files.add(MANIFEST)
    readme = '''CT20 초기 후보군 23개 — 2026-10-22

최신 후보 데이터: experiments/candidates_ct20_20261022/candidates.json
생성 과정: experiments/candidates_ct20_20261022/METHOD_KO.md
전체 후보 탐색: experiments/candidates_ct20_20261022/reports/candidate_report.html
후보별 단면·RI·왜곡·spot·MTF: 같은 reports의 candidate_layouts_and_mtf.pdf
전체 단면: 같은 reports의 candidate_contact_sheet.png

모든 독립 합격 후보를 유지했습니다. 하나의 winner를 선정하지 않았습니다.
각 렌즈 축상 중심 두께는 2 mm 초과, 20 mm 이하입니다.
고NA scalar 파동 결과의 underresolved/not_applicable 상태를 확인하십시오.
기하 MTF는 별도 진단이며 회절 학습 정답으로 인증한 데이터가 아닙니다.

코드 재현: python -m pip install -e .
native 빌드: python experiments/candidates_ct20_20261022/native/build.py
그림 재생성에는 matplotlib 및 한국어 폰트가 추가로 필요합니다.
보고서는 브라우저에서 로컬 파일로 열 수 있고 외부 API나 네트워크가 필요 없습니다.
이전 후보의 코드·입력·검사 이력은 candidates_20261022에 보존되어 있습니다.
이전 HTML/PDF는 별도 기존 ZIP에 있고, 이번 ZIP에는 중복 포함하지 않았습니다.
원본 코어 src/tests는 이 작업에서 변경하지 않았습니다. DB 생성이나 GitHub push는 하지 않았습니다.
'''
    OUTPUT.mkdir(parents=True, exist_ok=True)
    archive = OUTPUT / (NAME + '.zip')
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        z.writestr(NAME + '/PACKAGE_README.txt', readme)
        for path in sorted(files):
            z.write(path, NAME + '/' + path.relative_to(REPO).as_posix())
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
    result = {'archive_path': str(archive), 'bytes': archive.stat().st_size,
              'sha256': sha(archive.read_bytes()), 'crc_check_passed': True,
              'candidate_count': 23, 'file_count': len(files) + 1}
    (OUTPUT / (NAME + '.manifest.json')).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
