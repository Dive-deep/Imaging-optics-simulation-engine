# 독립 광학 엔진 리뷰 — 2차

먼저 `CODE_REVIEW.md`를 읽는다. 광학 물리·렌즈 설계·소프트웨어의 새 검토자가 확인한 결함 9개, 정상 동작 대조, 영향과 수정 우선순위를 요약했다. 각 역할의 세부 유도·출처·재현 조건은 다음 문서에 있다.

- `physics/review.md`
- `lens_design/review.md`
- `software/review.md`

`input/optical_engine/`은 **수정하지 않은 v0.1.0 배포 사본**이다. 이 리뷰 ZIP은 수정 버전이 아니다. 포함된 기존 결과는 당시 배포의 기록이며, 이번 검토 결과와 혼동하지 않는다. 고정 배포의 SHA-256 및 개별 파일 무결성 확인은 `reproductions/verification.json`에 기록했다.

## 재현

Python 3.10 이상과 NumPy/SciPy가 필요하다. 이 README가 있는 폴더에서:

```bash
python -m pip install -e input/optical_engine
python physics/repro_angular_otf_alias.py
python physics/repro_stale_bundle.py
python lens_design/reproduce_lens_findings.py
python software/reproduce_software_findings.py
```

각 스크립트는 인접한 고정 소스를 우선 import한다. 출력 JSON은 각 역할 폴더에 보존했고, 주 검토자의 독립 재실행 출력은 `reproductions/`에 있다. 위 명령은 **현재 결함을 재현**한다. 예외가 없거나 assertion이 통과했다는 것이 엔진의 올바른 계산을 승인한다는 뜻은 아니다.

기존 전체 183개 unittest를 별도로 실행한 로그는 `software/baseline_tests.log`다. 이 테스트들은 모두 통과했지만 새 반례들을 검증하지 않았다.

## 범위

P1 4건, P2 4건, P3 1건은 실제로 실행하여 확인한 항목이다. Scalar 회절·비축 투영·실제 무효부 단차·Zemax 수치 일치처럼 이미 선언된 모델 제한은 새 구현 결함과 분리했다. 엔진 수정과 수정 후 회귀 검증은 이 산출물에 포함되지 않는다.
