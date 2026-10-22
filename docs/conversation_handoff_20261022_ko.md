# 결상 계산 Engine&DB — 새 대화 인계

기준일: 2026-10-22 UTC. 저장소: [Dive-deep/Imaging-optics-simulation-engine](https://github.com/Dive-deep/Imaging-optics-simulation-engine).

새 프로젝트 대화에는 [붙여넣기용 시작 프롬프트](new_chat_prompt_ko.txt)를 사용한다. 짧게 시작하려면 다음 문장을 붙여넣는다.

> 프로젝트 「결상 계산 Engine&DB」를 이어가자. https://github.com/Dive-deep/Imaging-optics-simulation-engine 의 AGENTS.md, CODEX_START.md, HANDOFF.md, docs/conversation_handoff_20261022_ko.md와 CT20 실험 METHOD_KO.md를 읽어줘. 최신 기준은 각 렌즈 2 < CT ≤ 20 mm를 만족하는 독립 후보 23개 전부이며 단일 winner가 아니다. 기존 결과와 해시는 보존하고, 검증된 기하 결과와 미수렴·적용 불가 파동 결과를 구분해 현재 상태와 남은 작업을 복원한 뒤 내 다음 지시에 따라 진행해줘.

이 문서는 대화에서 이어가야 할 요구조건·결정·근거·미완료 작업을 정리한 인계 자료다. GitHub 연결만으로 이전 대화 원문이나 다른 프로젝트 대화의 기억이 자동 이관되지는 않는다. 저장소에 없는 과거 파일이나 실행은 복원된 것으로 간주하지 않는다.

## 1. 목표와 현재 작업 경계

장기 목표는 실제 광선 기반 결상 계산 엔진과 설계·성능 DB, 이후 검증된 성능 label을 사용하는 surrogate 작업으로 이어가는 것이다. 현재 단계는 고NA 3매 광학계의 다양한 초기 후보를 만들고 그 입력·출력·실패 상태·생성 이력을 보존한 것이다. **후보 주변 설계 DB, MySQL 서버/테이블, 학습된 surrogate는 아직 만들지 않았다.**

이번 후보 보고서의 범위와 별도로, 과거 사용자 요구에서 계속 유지할 사항은 다음과 같다.

- **장기 출력:** FOV, EPD, 0F/0.8F의 FFT MTF T/S 곡선 0–50 lp/mm, 0.5F/1F spot, distortion, NA, BFL, EFL, TTL, ray aberration, field curvature & distortion, grid distortion. 0.8F는 전체 field의 0.8배다. 열해석은 요청 범위에서 제외한다. 이번 후보 보고서의 MTF는 **0–20 lp/mm**이므로 일반 요구의 50 lp/mm까지 제공·검증했다고 해석하지 않는다.
- **ML 재료 입력 문헌 검토:** 렌즈 위치 A/B/C의 재료를 categorical ID로 넣는 방식, nd/Vd 수치로 넣는 방식, 파장별 n(λ)로 넣는 방식의 차이와 학술 근거를 검토해 달라는 요청이 있었다. 현재 보존된 기록만으로 최종 문헌 검토 결론이나 학습 성능을 입증하지 않으므로 **별도 검토 과제로 남긴다.** 네 재료 등록과 normal-line 계산의 완료는 이 비교 연구의 완료를 뜻하지 않는다.
- **과거 속도 목표:** STOP pupil 16×48과 FFT 128×128을 고정하고, 저장시간을 제외한 100개 계산시간을 기존 대비 절반으로 줄이는 목표와 FP16 사용 요청이 있었다. 관련 미업로드 실행·구현은 현재 환경에서 복원되지 않아 달성을 주장하지 않는다. 과거 언급된 외부 **0.136 s**는 축소 pupil과 다른 출력 범위를 사용하고 focus/ray fan/grid/2D PSF/JSON/file 처리를 포함하지 않은 수치이므로 현재 전체 작업과 직접 비교하지 않는다. 정확도와 속도 범위는 §5처럼 구분한다.

원본 엔진은 루트 `src/sequential_optics/`의 v0.2.1이며 `tests/`와 함께 이번 후보 생성 작업에서 변경하지 않았다. 새 탐색·native 커널·평가 코드는 `experiments/`에 분리했다. `audit_round2/input/optical_engine/`는 과거 고정 사본이다. 실질적인 엔진 변경 시 물리·렌즈 설계·소프트웨어 역할의 독립 리뷰를 수행하는 계약은 [AGENTS.md](../AGENTS.md)를 따른다.

## 2. 유지할 설계 요구조건

상세 입력은 [CT20 requirements.json](../experiments/candidates_ct20_20261022/requirements.json)과 [METHOD_KO.md](../experiments/candidates_ct20_20261022/METHOD_KO.md)가 기준이다.

| 항목 | 확정한 조건과 해석 |
|---|---|
| 물체·화각 | 무한 물체, ±12° = 전체 24°. 성능 field는 0/6/9.6/12° = 0/0.5/0.8/1F |
| 매수·재료 배치 | 공기 간격을 둔 3매, PMMA/PC 중 두 장 + FDS18W/TAFD45 중 한 장. 접합 없음 |
| 면 형식 | 플라스틱 standard 또는 even asphere의 conic/A4/A6/A8/A10/A12. 유리는 standard 구면, conic=0, 비구면 계수=0 |
| 각 렌즈 중심 두께 | **2 mm < CT ≤ 20 mm**. 각 매의 축상 두께이며 총 길이나 주변 두께 상한이 아님 |
| 간격·BFL | 두 축상 air gap ≥1 mm. 마지막 꼭짓점부터 실제 공통 상면까지 BFL ≥3 mm, 별도로 588 nm 근축 BFL ≥3 mm |
| 구경 | EPD 50 mm 또는 첫 렌즈 전면 clear diameter 50 mm 대안. 둘을 같은 조건으로 부르지 않음 |
| STOP | 첫 렌즈 앞, EFL보다 가까운 위치를 탐색. 실제 면 간섭과 광선 통과를 검사 |
| NA | 세 파장 각각 축상 전체 rim의 최소 image-space NA ≥0.7. 주변 field 모든 방향의 cone NA ≥0.7 조건은 아님 |
| 상면·상고 | 모든 field/RGB가 공유하는 한 평면. 588 nm의 12° chief 상고 +6.4 mm. 0.001 mm는 구현상 gate이며 사용자가 승인한 설계 허용오차가 아님 |
| 파장·가중치 | 이번 후보 실험은 486/588/656 nm, 입사 RGB power 동등 가중. 코어 예제 기본값 450/550/650 nm와 혼동하지 않음 |
| 제외 | 이번 실험은 열해석을 하지 않음. Detector crop·pixel response·실제 코팅·제조 공차·기구 제조성 인증도 없음 |

후보는 최종 우승 설계 한 개를 고르는 목적이 아니다. 합격한 서로 다른 부모의 후보를 모두 보존한다. ID와 유리 위치·재료·seed 순서의 표시는 성능 순위가 아니다. 전 시야 고MTF나 센서 안에 전체 RGB spot이 들어간다는 조건을 달성했다고 해석하지 않는다.

## 3. 재료 정의: nd / Vd와 명시적 normal-line 근사

사용자 입력 n의 최종 정의는 **nd**, 기준선 d=587.5618 nm다. n(550 nm)로 다시 해석하지 않는다. V는 **Vd**다.

| 재료 | nd | Vd |
|---|---:|---:|
| PMMA | 1.492 | 57.47 |
| PC | 1.586 | 32.00 |
| FDS18W | 1.94594 | 17.98 |
| TAFD45 | 1.95375 | 32.32 |

사용자가 채택한 모델은 `PgF = 0.6438 − 0.001682 Vd`, `nF−nC = (nd−1)/Vd`, `n(λ)=a+b/λ²+c/λ⁴`다. λ는 μm이며 nd, F−C, g−F 세 조건으로 Cauchy3 계수를 구한다. nominal 기준선은 d587.5618/F486.1327/C656.2725/g435.8343 nm다. 이번 후보는 normal line의 ΔPgF=0 가정을 사용한다.

이는 제조사의 측정 분산·Sellmeier 데이터나 Zemax Model Glass 고유식을 복원한 것이 아니다. 재료 이름을 이유로 사용자 nd/Vd를 catalog 값으로 바꾸지 않는다. Absolute index와 스펙트럼 좌표는 기존 nominal convention을 따른 것이며 원래 측정의 공기/진공 기준을 실증하거나 실제 공기→진공 변환을 수행한 것은 아니다. 열·흡수·코팅 자료도 이 두 수치에서 추정하지 않는다. 현재 체크아웃에서는 이 모델을 **별도 실험 코드와 후보 JSON**에서 확인하며, 과거 대화의 미업로드 코어 확장을 이 저장소에 존재한다고 가정하지 않는다.

## 4. 최신 CT20 결과와 이전 캠페인 보존

최신 경로는 [experiments/candidates_ct20_20261022/](../experiments/candidates_ct20_20261022/)다. 이전 탐색의 기하 합격 부모 38개를 CT20 조건으로 다시 계산했다. 합격 checkpoint 26개를 부모별로 묶어 **서로 다른 seed 23개를 전부 유지**했다. 같은 부모의 중복 checkpoint 선택에 기하 진단 점수를 사용했지만 후보군 전체의 순위나 winner는 정하지 않았다.

| 최신 후보군 항목 | 실제 결과 |
|---|---|
| 독립 후보·재료 순서 | 23개 / 16종 |
| 구경 모드 | 첫 전면 clear diameter 50 mm 대안 22개, EPD 50 mm 대안 1개 |
| 유리 위치 앞/중간/뒤 | 10 / 4 / 9개 |
| 실제 각 매 CT | 2.03175799273724–19.999997036672667 mm |
| 실제 EPD | 43.000474153825074–50.0 mm |
| RGB 축상 complete-rim NA | 0.7000361002586907–0.7715564534728826 |
| 독립 물리 검사 | 23/23개 통과, 새 광선 732,849개 실패 0 |
| 원본 엔진/native 대조 | 23개, 6,900개 광선, 상태·수치 실패 0; 원본 interval geometry 23/23개 통과 |

검증 근거는 [물리 원자료](../experiments/candidates_ct20_20261022/results/final_selection_physics_audit.json), [통계](../experiments/candidates_ct20_20261022/results/final_selection_physics_statistics.json), [원본 엔진 대조](../experiments/candidates_ct20_20261022/results/native_final_candidates/summary.json), [역할별 리뷰](../experiments/candidates_ct20_20261022/reviews/)다. 검사한 표본과 처방의 결과이며 연속 설계 공간 전체나 모든 미검사 광선의 증명은 아니다.

별도의 [경로 이전 재현 검사](../deliverables/20261022/relocated_checkout_validation.json)는 업로드 대상으로 선택한 파일을 임시 경로에 export한 뒤 native 커널을 새로 빌드하고, 동결 처방 23개의 hash 확인 및 4 field×RGB×16×48 pupil 재추적을 수행했다. 실제 광선 수는 chief/rim을 포함해 **225,492개, 실패 0**이다. 이는 기존 732,849개 독립 물리 검사를 대체하거나 합산한 결과가 아니며, 파동 수렴을 검증한 실행도 아니다.

이전 [candidates_20261022/](../experiments/candidates_20261022/)의 최종 24개·전체 탐색 이력·보고서는 보존했다. 이전 24개와 최신 CT20 23개는 다른 처방 집합이다. 이전 광선·MTF 결과를 CT20 처방의 결과로 바꿔 적지 않는다. 현재 환경에서 복원한 원본 코어의 **277개 테스트 통과** 기록은 [restored_core_tests.json](../experiments/candidates_20261022/results/restored_core_tests.json)에 있다. 이번 인계 문서 작성 때문에 이를 재실행하거나 결과를 갱신하지 않았다.

## 5. FP64 검증, 과거 속도 요청, 0.5% 조건

이전 대화의 FP16/속도 최적화 요청과 미업로드 구현 파일은 현재 새 환경에서 복원되지 않았다. 실제 확보한 출발점은 GitHub v0.2.1이다. 이번 후보는 별도 **FP64 C 커널**로 실제 교점·벡터 Snell·OPL을 계산하고 원본 엔진과 대조했다. FP16에서 후보 정확도를 검증했다거나 과거 128×128 FFT benchmark를 재현했다고 말하지 않는다. 현재 작은 native smoke의 속도 숫자를 전체 광학 평가·최적화 시간 또는 모든 유효 광선 처리량으로 일반화하지 않는다.

0.5% 조건은 평균 오차가 낮다는 뜻으로 완화하지 않는다. 속도·정밀도 변경을 검증할 때는 비교 기준, 지표, 파장·시야·주파수·실패 상태를 고정하고 **검사 대상의 최악 경우**를 확인해야 한다. 기존 기록만으로 모든 광학 성능의 물리 정확도가 0.5% 안이라고 주장할 수 없다.

현재 파동 코드가 실제 수행한 gate는 제한적이다. 같은 scalar 모델의 128→256 적분 격자 비교에서 각 field/파장의 요청 T/S 주파수 전부에 대해 `100 × |MTF256−MTF128| / |MTF256|`의 최댓값이 0.5 이하인지 검사한다. 256 결과가 0인데 차이가 0이 아니면 실패한다. 모든 쌍의 평균이나 절대 MTF 차이 0.005 기준으로 대체하지 않는다. 이 gate는 엄밀한 오차 상한·ray mesh 수렴·벡터 회절·실측/Zemax 상관 검증이 아니다.

## 6. 고NA 파동 결과는 아직 확정 label이 아님

주 경로는 physical-UV 역매핑 OPD와 pupil autocorrelation의 연속 주파수 scalar 적분이다. 고NA에서 균일 angular pupil amplitude와 chief-normal 국소 투영을 가정하며, boundary ray-map fold/caustic이 생기면 단일값 pupil 경로의 적용 불가를 남긴다. 이를 실제 FFT라고 부르지 않는다.

최신 23개 ×4 field ×3 파장의 **T/S 276쌍**은 `underresolved_quadrature` 200, `conditional` 3, `not_applicable` 73이다. 3쌍의 gate 통과도 해당 scalar 적분 비교에 한정한다. 별도 1024 FFT 진단은 588 nm의 0/0.8F, 총 46개를 계획했고 primary 적용 조건 때문에 일부를 생략하여 저장된 39개는 `underresolved` 17, `not_applicable` 22다. 고NA 파동 계산이 모두 해결됐다는 결과가 아니다.

모든 field에는 회절 없는 geometric ray-intercept OTF/MTF 진단을 따로 저장했다. 보고서의 92개 field 패널 중 64개는 저장된 scalar 파동 곡선, 28개는 명시적인 기하 진단이다. RGB 중 어느 파장의 파동 곡선이 없으면 패널 전체를 기하 진단으로 표시한다. 원래 실패 status를 성공으로 바꾸거나 0 power로 숨기지 않는다. 이 결과 전부를 diffraction-MTF surrogate의 정답 label로 바로 사용하면 안 된다.

## 7. 읽을 파일과 전달 산출물

- [CODEX_START.md](../CODEX_START.md), [HANDOFF.md](../HANDOFF.md), [README.md](../README.md): 진입점과 코어 계약.
- [CT20 candidates.json](../experiments/candidates_ct20_20261022/candidates.json): 23개 처방, 재료, 계산 조건, 출력, status와 이력.
- [CT20 METHOD_KO.md](../experiments/candidates_ct20_20261022/METHOD_KO.md): 생성·재계산·실패·수식·재현 순서.
- [frozen_selection.json](../experiments/candidates_ct20_20261022/results/frozen_selection.json), [final_summary.json](../experiments/candidates_ct20_20261022/results/final_summary.json): 동결 처방 hash와 집계.
- [전체 23개 HTML](../experiments/candidates_ct20_20261022/reports/candidate_report.html), [PDF](../experiments/candidates_ct20_20261022/reports/candidate_layouts_and_mtf.pdf), [전체 단면 모음](../experiments/candidates_ct20_20261022/reports/candidate_contact_sheet.png): 오프라인 탐색·출력용 자료.
- [이전 24개 방법과 결과](../experiments/candidates_20261022/METHOD_KO.md): 이전 세대 기록. CT20로 덮어쓰지 않음.

인계 작성 시점의 로컬 전달 ZIP은 다음과 같다. GitHub 전달 위치는 저장소 `deliverables/20261022/`의 [최신 CT20 23개 ZIP](https://github.com/Dive-deep/Imaging-optics-simulation-engine/raw/refs/heads/main/deliverables/20261022/optics_initial_candidates_ct20_23_20261022.zip)과 [이전 24개 ZIP](https://github.com/Dive-deep/Imaging-optics-simulation-engine/raw/refs/heads/main/deliverables/20261022/optics_initial_candidates_24_20261022.zip)으로 정했다. 업로드 완료 여부는 실제 원격 commit과 ZIP blob에서 확인하며 로컬 생성이나 링크 지정만으로 완료됐다고 말하지 않는다. 전달 목록은 [deliverables/20261022/README.md](../deliverables/20261022/README.md)를 확인한다.

| 묶음 | 로컬 파일 | 크기 | SHA-256 |
|---|---|---:|---|
| 최신 CT20 23개 | `/workspace/outputs/optics_initial_candidates_ct20_23_20261022.zip` | 40,131,585 bytes | `db9e66c27f5eb0b507a3131396c9f64cb3dfea68be72ad8540f0b83bb397dc18` |
| 이전 24개 | `/workspace/outputs/optics_initial_candidates_24_20261022.zip` | 27,149,380 bytes | `7a8c0dc208489b84e5653944f9fc22e9d2642c7783d80aafe0a40665d620f6c7` |

ZIP manifest에는 CRC 검사 통과가 기록되어 있다. 이후 인계 문서 수정분이 앞서 생성된 ZIP에 자동 반영되는 것은 아니므로 최신 작업 상태는 저장소 문서를 우선한다. 재생성·추가 실험은 새 경로와 새 hash로 남기며 기존 결과 파일이나 증거를 덮어쓰지 않는다.

두 실험의 code/data/results/reviews/METHOD와 HTML/PDF/전체 단면은 Git 전달 대상으로 정했다. 컴파일된 `.so`, 캐시, 최종 `reports/plots/`의 중복 개별 그림은 Git에서 제외한다. HTML에는 개별 그림이 내장되고 PDF와 전체 단면도 별도 제공하며 ZIP 원본은 보존한다. 새 환경의 CT20 native 커널은 저장소 루트에서 `python experiments/candidates_ct20_20261022/native/build.py`로 빌드한다. `native/README.md`의 이전 실험 경로 흔적을 최신 CT20 명령으로 오해하지 않는다.

위 `/workspace/outputs`와 기존 JSON·manifest에 남은 절대 경로는 **당시 환경의 경로**다. 새 대화 환경에 자동으로 존재하지 않는다. 저장된 후보 처방 읽기·추적과, 과거 부모 원자료까지 찾아 전체 탐색 이력을 재실행하는 일은 구분한다. 후자에는 이전 `/workspace/...` prefix를 현재 체크아웃으로 대응시키는 provenance resolver가 필요할 수 있다. 기존 METHOD·JSON·manifest의 문자열이나 hash를 고쳐 이 환경 차이를 숨기지 않는다.

## 8. 다음 대화의 진행 기준

다음 실행 범위는 새 대화의 사용자 지시에 따른다. 이어갈 후보 작업은 **23개 주변의 설계 스캔**과 **검증된 성능 label의 정의**다. 하나를 임의 winner로 줄이거나 새 지시 없이 대규모 DB·최적화·학습을 완료했다고 가정하지 않는다.

주변 스캔을 요청하면 재료 순서·면 형식·부모 ID를 유지하며 curvature, 정규화 비구면 계수, CT, air gap, STOP, 공통 상면의 범위와 sampling을 정한다. 보간 처방도 geometry → 광선 상태 → NA/상고/BFL/CT → 성능 품질을 다시 검사해야 한다. 합격 부모 사이의 보간은 자동 합격이 아니다.

Label 작업에서는 지표 정의·단위·센서/상면·파장·field·주파수·수렴 조건을 먼저 고정한다. 처방과 모델 hash, 근사 가정, 실패 이유를 함께 저장하고 feasibility와 성능을 구분한다. `underresolved`, `not_applicable`, 수치 실패를 성공 label로 채택하지 않는다. 현재 고NA 파동 경로의 적용성과 수렴을 해결할지, 별도 검증한 기하 지표만 우선 사용할지는 사용자 목표에 맞춰 정한다.

새 환경 검증은 기존 증거를 보존하는 별도 출력 경로에서 수행한다. `run_tests.py`는 `results/test_results.json`과 `results/test_run.log`를 덮어쓰므로 그대로 실행해 과거 기록을 새 결과로 바꾸지 않는다. 코어 회귀를 다시 확인할 때는 `PYTHONPATH=src:examples python -m unittest discover -s tests -v`로 실행하고 새 로그를 별도 경로에 저장하거나 임시 checkout에서 원래 runner를 사용한다. 설치·native build·검증을 실제로 실행하지 않았다면 완료했다고 보고하지 않는다.
