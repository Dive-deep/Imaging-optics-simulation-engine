# CT20 복구·선정 경로 독립 소프트웨어 검토

검토일: 2026-10-22. 검토자: `candidate_native_kernel` 서브에이전트.
범위는 새 `experiments/candidates_ct20_20261022/`의 복구 최적화, 최종 입력 계약, 보고서 연결이다. 원본 실험, 원본 Python 코어, native 커널 및 현재 실행 중인 `repair_ct20.py`/`rescue_ct20.py`는 수정하지 않았다.

## 결론

최적화 중의 낮은 샘플 수, 완화된 aperture 및 geometric MTF proxy는 후보 이동을 위한 휴리스틱이다. 통과 판정은 실제 clear aperture에서 새로 계산한 13개 field, RGB, STOP 16×48 및 chief/rim 광선과 실제 CT 값으로 이루어진다. 이 구분을 유지하는 현재 경로에서 실패 광선을 제거하여 후보를 통과시키는 우회 경로는 확인하지 못했다.

다만 기존 `rescue_ct20.py`의 일괄 1.25배 aperture 확장은 정상 설계의 conic domain을 벗어나 전체 광선을 수치 실패로 만들 수 있다. 별도 `rescue_ct20_guarded.py`에서 이 문제를 완화했으며 기존 실행 파일은 유지했다. 개선본의 검증은 domain 제어 회귀검증이고, 최적화 수렴이나 최종 후보 합격을 보증하지 않는다.

## 확인된 문제와 대응

### 1. 유효한 물리 aperture가 가상 확대 때문에 무효가 되는 문제

native 교차 계산은 해당 면 전체 clear radius에서 `1-(1+k)c²r² > 0`을 요구한다. `k>-1`, `c!=0`인 면의 한계 반경은 `1/(|c|sqrt(1+k))`이다. 기존 일괄 확대는 이 한계를 넘을 수 있다. 검토 시점에 repair/qualified_repair의 6개 서로 다른 설계에서 이를 관찰했다: seed0011, seed0036, seed0063, seed0072, seed0075, seed0098. 이 상태에서 invalid-ray residual이 거의 평탄해져 복구를 방해할 수 있다.

개선본은 trace에만 쓰는 deepcopy의 반경을 `max(actual, min(1.25*actual, domain*(1-1e-8)))`로 설정한다. 실제 반경이 이미 domain 밖이면 그대로 유지하여 기존 실패 판정을 보존한다. plane과 `k<=-1`인 면에는 유한 domain cap을 적용하지 않는다. geometry, STOP/image 여유 및 footprint penalty는 원래 물리 설계를 참조하고, 최종 `fit_apertures`/quick metrics/proxy는 실제 clear로 재추적한다.

`validate_rescue_domain.py`와 `results/rescue_guarded_domain_controls.json`의 검증 결과:

- 물리 반경을 축소하지 않고 원본 설계 및 곡률/conic/비구면 계수를 유지했다.
- domain에 가까운 정상 구면의 두 광선은 실제 aperture와 guarded aperture에서 성공했으며, 무조건 1.25배에서는 status 3이었다.
- 실제 반경이 domain 밖인 면의 실패를 그대로 유지했다.
- 실제 seed0075의 300개 표본 광선은 물리 aperture에서 300개 성공, 기존 일괄 확대에서 300개 status 3, guarded에서 300개 성공이었다.
- 실제 설계의 guarded residual은 유한했고 입력 설계는 변하지 않았다.

### 2. 평균 clipping penalty만으로 개별 외곽 광선이 남는 문제

guarded 개선본은 기존 광선별 penalty 외에 각 면의 RGB×field×ray 전체에서 최대인 `(footprint-actual_clear+0.10mm)`의 양수 부분에 20을 곱한 residual 6개를 추가한다. 이 변경은 모든 면의 가장 큰 초과를 직접 억제한다. 기본 merit는 8×32, `max_nfev=150`이며 CLI에서 변경할 수 있다. NA 목표는 0.706이고 축상 rim 전체의 최솟값에 큰 penalty를 적용한다. 유한한 penalty는 제약 충족의 증거가 아니므로 실제 13-field 전 광선 통과를 반드시 후속 판정으로 사용한다.

### 3. 작은 provenance 해석상 주의점

`repair_ct20.py`의 `common_image_plane_delta_mm`는 계산된 이동량이며 허용 범위 밖이라 적용하지 않았어도 기록될 수 있다. 최종 `design.axial_distances_mm[5]`가 실제 적용값의 근거이다. 현재 실행 파일은 수정하지 않았다.

`stop_scan_ct20.py`는 STOP만 변경한 trial을 만들지만 후속 `qualify()`가 clear radii를 추가 축소할 수 있다. 따라서 `fixed_CT_and_surface_prescription`은 곡률, conic, A4–A12, 축상 간격 및 image plane 고정의 의미로 읽어야 한다. 최종 실제 clear radii는 qualification 결과 design에 존재한다. 보고서 담당이 최종 generation history에 `ct20_stop_scan` 보존을 반영했으며 원래 STOP 반경과 변경량도 유지한다.

## 물리·광선 계약 검토

- 3개 렌즈, 6면의 순서와 재료는 유지된다. PMMA/PC 면만 conic과 A4–A12를 최적화하고 glass는 standard 구면이다.
- Parameterization의 최종 CT 상한은 20 mm, 하한은 2.02 mm이고 fresh 판정은 실제 값에 대해 정확히 `>2`와 `<=20`을 확인한다. 기존 label을 재사용하지 않는다.
- trace 배열은 failed ray를 제거하지 않는다. status와 image_forward가 모두 유효한 ray만 transmitted centroid/OTF에 기여하며 invalid residual과 footprint residual이 별도로 남는다.
- RGB image centroid는 색마다 재정렬하지 않고 field별 하나를 사용한다. OTF도 RGB 복소값을 합성한다. 이 geometric proxy는 회절 FFT MTF 자체가 아니다.
- 축상 NA 수락값은 3파장 모두 complete-rim의 최솟값이다. `quick_metrics`의 all-sample maximum NA는 진단값이며 최종 수락값으로 쓰지 않는다.
- glass 등 물성은 기존 사용자 승인 nd/Vd normal-line 모델을 그대로 사용한다. 신규 열해석 속성이나 제조사 분산 정밀도를 만들어 넣지 않는다.
- `fit_apertures(coupled=False)`는 원래 실제 추적과 aperture 변경 후 추적의 status, 방향, 위치, OPL, image 좌표를 비교한다. 기존에 막힌 ray를 제거해 성공으로 바꾸는 단계가 아니며, 원래 전 광선 통과 여부도 최종 eligible에서 검사한다.
- STOP 스캔은 사용자 허용 `front_diameter50` 분기에서만 수행하고 EPD50 분기는 건드리지 않는다. ±0.60 mm scan은 새 STOP에 대해 4-field 사전검사 뒤 13-field fresh qualification을 시행한다. STOP 크기와 실제 EPD는 새 결과에 명시된다.

## 전체 후보 보존과 최종 감사

`prepare_candidates.py`는 점수로 전체 후보를 잘라내지 않고 같은 부모의 checkpoint끼리만 선택한다. 모든 eligible distinct parent를 보존하며 최소 20개를 요구한다. persistent `CT20_XXXX` ID를 사용하고 glass 위치/재료/source ID 순서는 품질 rank가 아니다. actual CT, fresh CT gates, 전체 metric gates, RGB NA 및 original all-ray 통과를 확인한다.

보고서 담당이 `finalize_candidates.py`의 모든 geometric fallback curve에도 주파수 배열, T/S 길이, finite, DC 정규화 및 20 lp/mm 범위를 검증하도록 보강했다. 수정 담당에게서 정상 계산 schema 통과와 NaN/길이/DC/주파수 오류 거부 결과를 전달받았다. 이는 직접 실행한 optical 검증과 구분한다.

새 `native/audit_final_candidates.py`는 기존 24개 전용 감사의 별도 사본이며 N≥20에 동작한다. frozen selection의 실제 CT, fresh gate, distinct parent, 전체 eligible parent 보존 및 design SHA를 읽는 단계에서 재검사한다. 설계마다 chief 1개, 중간 pupil 16개, rim 8개를 3파장×4field에서 비교하므로 총 300N rays이다. 원본 Python 엔진의 연속 radial-interval geometry 검사와 native의 매면 위치/방향/OPL 및 image-plane 좌표를 대조한다. 한 ray라도 실패하면 해당 후보와 전체 audit는 실패한다. input과 관련 source SHA의 전후 변경도 검사한다.

입력 계약 회귀검증 `results/native_audit_input_controls.json`에서는 현재 working pool 21개를 전부 받아들였고, CT 20초과/CT 2와 동일/fresh upper gate 누락/부모 중복/pool 일부 누락/19개 입력을 모두 거부했다. 이 입력 control 이후 23개 selection을 동결했고, 아래 최종 광선 감사를 별도로 실행했다.

## 동결된 최종 23개 원본 엔진 대조 결과

`results/frozen_selection.json`의 후보 23개를 하나도 제외하지 않고 `native/audit_final_candidates.py`로 검사했다. CPU 프로세스 1개에 BLAS/OMP/MKL 스레드 1개 제한을 적용했다. 검사 소요 시간은 21.769473초이며, 이는 후보 생성이나 wave 계산 benchmark가 아닌 독립 ray/geometry 감사 시간이다.

- 원본 Python 엔진의 연속 radial-interval geometry 검사: **23/23 통과**.
- RGB 3파장 × field 0/0.5/0.8/1.0 × chief·중간·rim 25 rays × 23개: **6,900/6,900 성공**.
- 상태, 수치 한계 또는 image projection 실패: **0건**.
- 매면 XYZ 최대 절대 차이: **2.795275122e-11 mm** (허용 1e-8 mm).
- 매면 방향 성분 최대 절대 차이: **1.639494096e-12** (허용 2e-10).
- 매면 OPL 최대 절대 차이: **6.993161605e-11 mm** (허용 1e-8 mm).
- 최종 image XY 최대 절대 차이: **3.983213759e-11 mm** (허용 1e-8 mm).
- frozen selection 및 원본/core/model/native/audit 관련 source SHA 전후 변경: **0건**. 모든 후보의 canonical design SHA도 동결 입력과 일치했다.

근거 파일은 `results/native_final_candidates/plan.json`, `summary.json`과 후보별 `CT20_XXXX.json` 23개이다. 입력 frozen selection SHA-256은 `19a6f2e21fff64b2a36947426ffb42a75683617fb3fd41971f7fcf396464497d`이다. `summary.json` SHA-256은 `974e81d06d5218a028ba1f1122ece9d9e2ed5e92f8d2f230fa7165138925cb82`이다.

이 결과는 해당 표본의 원본·native 광선 교차/굴절/OPL 구현 일치와 원본 geometry 판정을 확인한다. 13-field 전체 pupil의 fresh qualification, FFT MTF sampling 수렴 및 제조 검증은 각각 별도의 검증이며 본 감사 결과로 대체하지 않는다.

## 남은 한계

4-field merit와 유한한 pupil sampling은 연속 design space/field/pupil 전체를 보증하지 않는다. geometric metric 합격은 제조 가능성, wave sampling 수렴, 고 NA 벡터 회절의 정확도, 최적해 또는 Zemax 수치 일치를 의미하지 않는다. candidate별 FFT/기하 proxy 구분과 wave diagnostic 상태는 별도 보고서대로 유지해야 한다.

## 검토 시점 파일 SHA-256

- `repair_ct20.py`: `24cb1e073ed22c9486d8810cc479806e3581a00f40ae969bd490fc08e07fd3a6`
- `rescue_ct20.py`: `901de3e48b2ef8ce29fa83f16e2ff65cf53a651fe55cb218e81a8bcd4c676d3d`
- `rescue_ct20_guarded.py`: `5c0964fadfde979db0fcec9cd39711ddbc0687c2e53694599d937e45d01e7175`
- `stop_scan_ct20.py`: `83f3736f679158b30295f47a5a300e4f9061d9a31fc27b64f60f043f48bd4b39`
- `native/audit_final_candidates.py`: `85b2327806fbb0e72f96fc8f13591d6ed77973db1b244f144131023363058c66`
