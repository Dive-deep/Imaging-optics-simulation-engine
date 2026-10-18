# 현재 작업 상태 — Sequential Optics Core v0.2.0

2026-10-18. 사용자가 Cloud 마이그레이션을 보류하고 후속 리뷰의 R1–R9 수정 및 기본 원리 문서를 요청했다. GitHub push나 Cloud 연결은 수행하지 않았다.

## 완료 상태

- R1–R9 모두 수정했다. 각 결함의 기대 동작을 검사하는 회귀 테스트 52개를 추가했다.
- 기존 183개를 포함한 전체 **235개 테스트 통과**, 실패·오류·skip 0.
- 구현자와 별도인 광학 물리학자·렌즈 설계 엔지니어·소프트웨어 개발자 역할의 에이전트가 다시 검토했다. 추가로 발견한 BoundaryState 변경과 curved STOP rim 문제도 수정·재검증했다.
- 기본 원리 문서는 `docs/core_principles_ko.md`, 같은 내용의 LaTeX 원문은 `docs/core_principles_ko.tex`다. 수식 231개(블록 45, 인라인 186)가 변환 후 보존됨을 검사했다.
- 상세 수정·증거는 `FIX9_REPORT.md`, 테스트 기록은 `results/test_results.json`, 독립 리뷰는 `reviews/fix9_*_review.md`에 있다.

## 다음 작업자가 유지할 계약

1. 실제 비구면 교점·Snell·OPL은 real ray다. EFL/BFL는 근축 정의량이며 회절은 명시적 scalar 근사다.
2. G 기본 파장은 설정으로 관리하며 불필요한 RGB·열·흡수 데이터를 요구하지 않는다.
3. 렌즈 단독 MTF와 실제 Lambertian source의 RI/영상을 구별한다. Image lp/mm와 angular cycles/rad도 구별한다.
4. 균일 온도에 따른 n(T)/자유 팽창만 현재 범위다. 제조 공차·응력·비균일 온도·GRIN·mount FEA는 포함하지 않는다.
5. `conditional`은 정확도 인증이 아니다. `underresolved` 및 수치 실패를 물리적 차단이나 0 전력으로 숨기지 않는다.
6. 소스 개발 대상은 루트 `src/`와 `tests/`다. `audit_round2/input/optical_engine/`는 수정 전 고정 사본으로 보존한다. 여기의 재현 스크립트는 버그가 존재함을 확인하므로 수정 검증은 루트 테스트를 사용한다.
7. `history/`는 v0.1.0 당시 상태다. 최신 판정은 루트 문서가 우선한다.

Zemax/Code V 및 실제 측정과 수치 상관 검증은 수행하지 않았다. 해당 실제 처방의 데이터 생성에는 ray/grid/padding/source 적분 수렴과 모델 적용성을 확인해야 한다.
