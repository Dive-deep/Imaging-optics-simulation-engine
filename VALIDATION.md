# v0.2.1 검증 상태

2026-10-18, Python 3.12.14 / NumPy 2.3.5 / SciPy 1.17.0.

GitHub 이전 폴더에서 277개 테스트를 다시 실행했다. 원래 v0.2.1 배포의 99.641초 검증 기록은 `history/VALIDATION_v0.2.1_release.md`, `history/test_results_v0.2.1_release.json`, `history/test_run_v0.2.1_release.log`에 보존했다. 루트 results의 테스트 기록은 이번 재실행 결과다. 예제/CLI 및 독립 물리 리뷰의 기록은 원래 배포 때 얻은 것이며 이번에 다시 실행한 것으로 표시하지 않는다.

- 기존 235개 + C1 13개 + C2 7개 + C3/C4 16개 + MTF 6개 = **277개 통과**.
- 실패·오류·skip 0. `python3 run_tests.py`, 마이그레이션 폴더에서 재실행한 전체 기록 시간 90.868초.
- 세 예제(`run_demo.py`, `run_source_wave.py`, `run_field_mtf.py`) 및 G-only CLI 모두 exit 0.
- `results/test_results.json`, `results/test_run.log`, `results/fix_c_integration.json`에 기록했다.
- 구현자와 다른 물리·렌즈 설계·소프트웨어 AI 에이전트가 검토했다. 첫 수정 후 발견한 rear STOP-chain aiming과 complex wavelength 입력도 수정·재검증했다.
- 독립 C1 30-ray chain, 전력12조건, C2 기울어진 convex cell12개 및 별도 소프트웨어 대조가 통과했다. 전체 suite 수와 중복 합산하지 않는다.
- Markdown/LaTeX의 231개 수식은 종류·순서·원문이 앞뒤 공백을 제외하고 동일하며 이전 수식도 그대로 보존했다.

자세한 범위는 `FIX_C_REPORT.md`와 `reviews/postfix_*_c.md`를 참조한다. D1 구현은 변경하지 않았다. 상용툴/측정 상관 및 모든 처방의 수렴을 인증하지 않는다. `underresolved`는 정답 데이터로 사용하지 않으며 `conditional`도 수렴/적용성 확인이 필요하다. 기존 scalar·chief-normal/local detector 모델의 한계를 유지한다.

과거 검증·리뷰 문서는 `history/`, `audit_round2/`, `FIX9_REPORT.md`에 보존했다. 현재 구현은 루트 `src/`다.
