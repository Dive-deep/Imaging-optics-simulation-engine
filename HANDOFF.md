# 결상 계산 Engine&DB — 현재 작업 상태

엔진 기준은 Sequential Optics Core v0.2.1이다. 새 프로젝트 대화에는 [시작 프롬프트](docs/new_chat_prompt_ko.txt)를 붙여넣고 [대화 인계 문서](docs/conversation_handoff_20261022_ko.md)를 읽는다. 저장소가 대화 원문이나 미업로드 속도 최적화 파일까지 자동 전달하지는 않는다.

> 이 저장소의 AGENTS.md, CODEX_START.md, HANDOFF.md와 docs/conversation_handoff_20261022_ko.md를 읽고 「결상 계산 Engine&DB」를 이어가자. 최신 CT20 후보 23개 전부를 유지하고, 완료·미검증·미완료 상태를 구분해 복원한 뒤 내 다음 지시에 따라 진행해줘.

## 2026-10-22: 각 렌즈 CT20 제한 후보군

- 최신 실험은 `experiments/candidates_ct20_20261022/`다. 각 렌즈의 **축상 중심 두께 2 < CT ≤ 20 mm** 조건으로 이전 기하 합격 부모 38개를 다시 탐색했다. 합격한 서로 다른 부모 **23개를 모두 유지**하며 전체 순위나 단일 winner를 정하지 않았다. 재료 순서 16종, 유리 위치별 앞/중간/뒤 10/4/9개다.
- 실제 CT 범위는 **2.0317579927–19.9999970367 mm**다. 첫 렌즈 전면 clear diameter 50 mm 대안 22개와 EPD 50 mm 대안 1개를 구분한다. 축상 complete-rim NA≥0.7은 RGB 각각 검사했고, 주변 field의 모든 방향 cone NA≥0.7을 강제하지 않았다. 6.4 mm 상고는 588 nm의 12° chief 기준이며 공통 상면을 사용한다. CT 상한은 총 길이나 주변 두께 상한이 아니다.
- [후보 JSON](experiments/candidates_ct20_20261022/candidates.json), [방법·범위](experiments/candidates_ct20_20261022/METHOD_KO.md), [동결 처방과 hash](experiments/candidates_ct20_20261022/results/frozen_selection.json), [최종 집계](experiments/candidates_ct20_20261022/results/final_summary.json)가 후속 작업의 기준이다. 전체 23개 그림은 [HTML](experiments/candidates_ct20_20261022/reports/candidate_report.html), [PDF](experiments/candidates_ct20_20261022/reports/candidate_layouts_and_mtf.pdf), [전체 단면 모음](experiments/candidates_ct20_20261022/reports/candidate_contact_sheet.png)에 연결한다. 최종 전달 묶음 경로는 [CT20 후보 23개 ZIP](https://github.com/Dive-deep/Imaging-optics-simulation-engine/releases/download/candidates-ct20-20261022/optics_initial_candidates_ct20_23_20261022.zip)이다.
- 독립 물리 검사는 **23/23개 통과, 새 광선 732,849개 실패 0**이다. Export 왕복·4097개 radial 기하 표본·저장 metric 대조를 포함하며 [검증 원자료](experiments/candidates_ct20_20261022/results/final_selection_physics_audit.json)와 [통계](experiments/candidates_ct20_20261022/results/final_selection_physics_statistics.json)에 보존했다. [원본 엔진/native 커널 대조](experiments/candidates_ct20_20261022/results/native_final_candidates/summary.json)는 **6,900개 광선 상태·수치 실패 0**, 원본 interval geometry 23/23개 통과다. 이는 검사한 처방과 표본에 대한 검증이며 고NA 파동 인증이 아니다.
- 주 파동 MTF는 scalar pupil autocorrelation 적분이며 actual FFT 진단과 구별한다. Field/파장별 T/S **276쌍 중 200 `underresolved_quadrature`, 3 `conditional`, 73 `not_applicable`**다. 적용 불가 패널에는 회절 없는 기하 MTF를 별도 진단으로 표시하고 원래 파동 실패 상태를 유지한다. 전 시야 고MTF 완성 설계나 수렴한 surrogate 학습 label로 취급하지 않는다. Sensor crop·코팅·공차·제조성 검증은 포함하지 않았다.
- 재료는 사용자 nd(d587.5618 nm)/Vd와 승인된 normal-line/Cauchy3 근사다. PMMA 1.492/57.47, PC 1.586/32, FDS18W 1.94594/17.98, TAFD45 1.95375/32.32를 유지한다. n(550 nm) 또는 제조사 실측 분산으로 바꾸지 않는다. 과거 FP16·속도 작업의 미업로드 파일은 복원되지 않았으며 이번 후보 검증은 별도 FP64 커널이다. 0.5% 최악 조건은 평균으로 완화하지 않고, 현재 확인한 것은 제한된 scalar 적분의 점별 격자 비교다.
- 원본 `src/`·`tests/`와 아래 이전 24개 실험을 수정하지 않았다. [이전 후보 24개 ZIP, 약 27 MB](https://github.com/Dive-deep/Imaging-optics-simulation-engine/releases/download/candidates-ct20-20261022/optics_initial_candidates_24_20261022.zip)도 삭제·덮어쓰기 없이 보존한다. 새 실험의 처방 hash가 다르면 이전 파동 probe 결과를 새 후보의 인증으로 재사용하지 않는다. 후보 생성 단계에서는 GitHub push나 DB 구축을 수행하지 않았다. 이후 인계 업로드 완료 여부는 실제 원격 commit으로 확인한다.
- DB·MySQL·학습된 surrogate는 아직 없다. 다음은 사용자 지시에 따른 주변 설계 스캔과 검증된 성능 label 정의다. 고NA 파동 모델 적용성과 수렴이 해결되지 않은 출력을 확정 정답으로 사용하지 않는다.

## 이전 2026-10-22: 고NA 3매 초기 후보 24개 실험

- 사용자 목표: infinity, FOV±12°, 588nm chief 상고6.4mm, 3매 air-spaced, PMMA/PC2매+구면glass1매, 첫 개구 또는EPD50mm, 축상NA≥0.7, BFL≥3mm, CT>2mm, gap≥1mm. 플라스틱은 conic/A4…A12의even asphere, 유리는standard만 허용한다.
- `experiments/candidates_20261022/`에96개 수식 기반+16개 공개특허 변환 시작점,112개 구경확대 탐색, 후속 OTF 개선·공통초점 시험을 보존했다.38개 기하 합격parent 중24개를 선정했다. 원본`src/`·`tests/`는 수정하지 않았다.
- 최종 입력은`candidates.json`, 구체적인 방법·범위·실패 해석은`METHOD_KO.md`, 후보별 광선도/spot/RI/distortion/0–20lp/mm T/S곡선은`reports/`에 있다. 처방과 검증hash는`results/frozen_selection.json` 기준이다.
- 최종24개 전부 첫 렌즈 clear diameter50mm 대안이며 실제EPD43.20–43.96mm다. 모든RGB 축상 complete-rimNA≥0.7을 확인했다. 주변 field 전체coneNA≥0.7 또는RGB 모두 상고6.4mm를 강제하지 않았다. Sensorcrop·coating·공차·제조성은 포함하지 않았다.
- 독립 물리 검사: 새764,712ray 실패0, export roundtrip·4097radial기하·저장metric대조24/24통과. 독립 원본 엔진 대조:7,200ray 상태불일치0, 원본interval geometry24/24통과. 새환경 원본277test통과. 근거는실험`results/final_selection_physics_audit.json`, `results/native_final24/`, `results/restored_core_tests.json`.
- 파동 MTFT/S는 조건부scalar autocorrelation이며 actualFFT진단과 구분한다. 일부field는boundary fold/caustic으로 적용불가, 다수는 격자미수렴이다. 모든field의 기하MTF진단은 제공하지만 이를 회절MTF성공으로 취급하지 않는다. 정확한 최종집계는`results/final_summary.json`. 아직 전 시야 고MTF 완성설계나 신뢰가 확인된 neural surrogate 학습label이 아니다.
- 새환경에서 과거 미업로드 최적화코드는 복원되지 않았다. 현재 별도FP64 nativekernel로 처방을 생성했다. 이 이전 후보 생성 단계에서는 GitHub push/DB구축을 수행하지 않았다.

## 이전 완료 기록

2026-10-18. C1–C4 수정, D1 그림 설명, field별 T/S MTF metadata와 예제를 완료했다. 사용자가 GitHub/Codex 이전을 승인했다. 저장소는 `Dive-deep/Imaging-optics-simulation-engine`이며 환경 설정은 `CODEX_START.md`를 따른다.

- 기존235 + 신규42 = 277개 통과. 실제 결과는 `results/test_results.json`.
- 별도 물리·렌즈 설계·소프트웨어 AI 리뷰 완료. 후속 STOP-chain/complex wavelength 문제도 수정됐다.
- `FIX_C_REPORT.md`, `VALIDATION.md`, `reviews/postfix_*_c.md`에 근거와 한계가 있다.
- D1 코드는 미수정. 성분별 `f_theta` chart와 방사 방향 F-θ의 차이를 `docs/D1_explanation_ko.md` 및 그림으로 설명한다.
- FFT MTF는 각 field의 고정 chief-normal X/Y 절단 T/S다. 메타데이터와 `docs/field_mtf_usage_ko.md`, `examples/run_field_mtf.py`를 따른다. Field-radial T/S나 임의 2D field의 Zemax 수치 동치를 주장하지 않는다.
- `docs/core_principles_ko.md/.tex`의 수식231개를 보존하면서 축 규약과 F-θ 정의 설명을 보완했다.
- `history/`는 v0.1.0/v0.2.0의 역사적 기록이다. `audit_round2/input/optical_engine`는 고정 사본이며 수정하지 않는다. 버그 존재 확인 script를 수정본의 acceptance test로 사용하지 않는다.

Real-ray 교점/Snell/OPL, 1차 EFL/BFL, scalar 회절의 역할을 구별한다. 기본 G와 미사용 선택 데이터의 독립성을 유지한다. Lens-only MTF와 실제 Lambertian source RI/영상은 별도 모델이다. 균일 n(T)/팽창 이외 응력·GRIN·mount FEA·제조 공차는 범위 밖이다. 수치 실패를 물리적 차단으로 감추지 않고 `underresolved`/`conditional`을 정확도 인증으로 해석하지 않는다.

코어를 수정하는 후속 작업은 루트 `src/`, `tests/`를 대상으로 하며, 후보 주변 스캔은 새 실험 경로에 분리한다. 남는 코어 개선사항은 D1 정의 변경, 영상 합성의 quality 전파 등 기존 독립 리뷰와 수정 보고서에 명시되어 있다. 현재 우선순위와 실행 범위는 새 대화의 사용자 지시에 따른다. 기존 결과·해시를 보존하며 이 저장소에서 이어간다.
