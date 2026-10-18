# v0.2.1 검증 상태

2026-10-18, Python 3.12.14 / NumPy 2.3.5 / SciPy 1.17.0.

- 기존 235개 + C1 13개 + C2 7개 + C3/C4 16개 + MTF 6개 = **277개 통과**.
- 실패·오류·skip 0. `python3 run_tests.py`, 전체 기록 시간 99.641초.
- 세 예제(`run_demo.py`, `run_source_wave.py`, `run_field_mtf.py`) 및 G-only CLI 모두 exit 0.
- `results/test_results.json`, `results/test_run.log`, `results/fix_c_integration.json`에 기록했다.
- 구현자와 다른 물리·렌즈 설계·소프트웨어 AI 에이전트가 검토했다. 첫 수정 후 발견한 rear STOP-chain aiming과 complex wavelength 입력도 수정·재검증했다.
- 독립 C1 30-ray chain, 전력12조건, C2 기울어진 convex cell12개 및 별도 소프트웨어 대조가 통과했다. 전체 suite 수와 중복 합산하지 않는다.
- Markdown/LaTeX의 231개 수식은 종류·순서·원문이 앞뒤 공백을 제외하고 동일하며 이전 수식도 그대로 보존했다.

자세한 범위는 `FIX_C_REPORT.md`와 `reviews/postfix_*_c.md`를 참조한다. D1 구현은 변경하지 않았다. 상용툴/측정 상관 및 모든 처방의 수렴을 인증하지 않는다. `underresolved`는 정답 데이터로 사용하지 않으며 `conditional`도 수렴/적용성 확인이 필요하다. 기존 scalar·chief-normal/local detector 모델의 한계를 유지한다.

과거 검증·리뷰 문서는 `history/`, `audit_round2/`, `FIX9_REPORT.md`에 보존했다. 현재 구현은 루트 `src/`다.
