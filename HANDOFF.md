# 현재 작업 상태 — Sequential Optics Core v0.2.1

2026-10-18. C1–C4 수정, D1 그림 설명, field별 T/S MTF metadata와 예제를 완료했다. 사용자가 GitHub/Codex 이전을 승인했다. 저장소는 `Dive-deep/Imaging-optics-simulation-engine`이며 환경 설정은 `CODEX_START.md`를 따른다.

- 기존235 + 신규42 = 277개 통과. 실제 결과는 `results/test_results.json`.
- 별도 물리·렌즈 설계·소프트웨어 AI 리뷰 완료. 후속 STOP-chain/complex wavelength 문제도 수정됐다.
- `FIX_C_REPORT.md`, `VALIDATION.md`, `reviews/postfix_*_c.md`에 근거와 한계가 있다.
- D1 코드는 미수정. 성분별 `f_theta` chart와 방사 방향 F-θ의 차이를 `docs/D1_explanation_ko.md` 및 그림으로 설명한다.
- FFT MTF는 각 field의 고정 chief-normal X/Y 절단 T/S다. 메타데이터와 `docs/field_mtf_usage_ko.md`, `examples/run_field_mtf.py`를 따른다. Field-radial T/S나 임의 2D field의 Zemax 수치 동치를 주장하지 않는다.
- `docs/core_principles_ko.md/.tex`의 수식231개를 보존하면서 축 규약과 F-θ 정의 설명을 보완했다.
- `history/`는 v0.1.0/v0.2.0의 역사적 기록이다. `audit_round2/input/optical_engine`는 고정 사본이며 수정하지 않는다. 버그 존재 확인 script를 수정본의 acceptance test로 사용하지 않는다.

Real-ray 교점/Snell/OPL, 1차 EFL/BFL, scalar 회절의 역할을 구별한다. 기본 G와 미사용 선택 데이터의 독립성을 유지한다. Lens-only MTF와 실제 Lambertian source RI/영상은 별도 모델이다. 균일 n(T)/팽창 이외 응력·GRIN·mount FEA·제조 공차는 범위 밖이다. 수치 실패를 물리적 차단으로 감추지 않고 `underresolved`/`conditional`을 정확도 인증으로 해석하지 않는다.

다음 작업은 루트 `src/`, `tests/`에서 수행한다. 남는 개선사항은 D1 정의 변경, 영상 합성의 quality 전파 등 기존 독립 리뷰와 수정 보고서에 명시되어 있다. 이 저장소에서 후속 작업을 이어가며 이전 작업을 반복하지 않는다.
