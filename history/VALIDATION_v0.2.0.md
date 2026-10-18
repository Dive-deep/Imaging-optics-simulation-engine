# v0.2.0 검증 상태

2026-10-18, Python 3.12.14 / NumPy 2.3.5 / SciPy 1.17.0.

- 원래 183개 테스트 + R1–R9 회귀 52개 = **235개 통과**.
- 실패 0, 오류 0, skip 0; 전체 실행 58.136초.
- 실제 실행: 저장소 루트에서 `python3 run_tests.py`.
- 추가 통합 실행: 두 예제 스크립트와 G-only CLI 모두 exit 0.
- 상세 결과: `results/test_results.json`, `results/test_run.log`.
- 수정/재검증: `FIX9_REPORT.md`와 `reviews/fix9_*_review.md`.

새 회귀 테스트는 R1 6개, R2 12개, R3/R4 13개, R5/R6 9개, R7/R8/R9 12개다. 독립 리뷰에서 별도로 실행한 부분 suite는 이 수와 중복되므로 합산하지 않는다.

기존 원형/직사각 pupil MTF, scalar RS, Snell, 두꺼운 렌즈 EFL/BFL, 균일 팽창 similarity, ray reciprocity/eikonal, 광량 ledger, G-only/누락 데이터 격리 등의 대조를 유지했다. 독립 에이전트 3명이 수정 코드를 다시 검토했고 발견한 동일 원인의 잔여 문제를 수정했다.

핵심 문서의 Markdown→LaTeX 변환 후 수식 231개(블록45, 인라인186)가 보존됐다. 기본 물리/계산 원리 검토와 실제 구현 테스트는 구별한다.

이 검증은 원래 확인된 9개 결함의 수정 및 명시된 대조에 관한 것이다. 모든 처방의 수렴, Maxwell 정확도, Zemax/Code V 호환, 실제 측정 상관을 보증하지 않는다. `underresolved`는 학습/최적화의 정답 데이터로 채택하지 않는다. `conditional`도 수렴 및 적용성 검토를 대체하지 않는다.

이전 v0.1.0 검증 문서는 `history/VALIDATION_v0.1.0.md`, 수정 전 독립 리뷰와 소스는 `audit_round2/`에 보존했다.
