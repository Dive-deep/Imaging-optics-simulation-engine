# 구현·검증 기록

이 패키지는 기술 명세를 실제 Python 계산 경로로 연결한 첫 구현이다. 테스트는 독립 해석해·보존 관계·수치 수렴·구체적 반례를 기준으로 한다. **실제 광학 처방의 Zemax 대조나 측정 비교는 수행하지 않았다.**

## 실행 결과

최종 전체 회귀 결과는 **183개 통과, 실패 0, 오류 0, skip 0**이다. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0에서 실행했다. 아래 기록은 해당 환경과 합성 표본 범위에 대한 검증이다.

최종 자동 테스트의 권위 있는 기록은 `results/test_results.json`과 `results/test_run.log`다. 개별 test ID, 통과/실패, 환경과 실행시간을 보존한다. `run_tests.py`로 재생성할 수 있다. 테스트 수만으로 정확도를 판단하지 않도록 아래 검증 범위와 반례를 함께 기록했다.

`pip install --no-deps --no-build-isolation -e .` 설치 및 `python -m sequential_optics ...` 명령행 해석을 실제로 실행했다. `compileall`도 통과했다. 설치 검증에 외부 광학 프로그램은 사용하지 않았다.

## 샘플 선정과 검증 범위

| 표본 | 선택한 이유와 확인 대상 |
|---|---|
| 평면·평행 평판 | Snell 각, 정확한 lateral shift/OPL, incident phase, Fresnel·흡수, 역방향 상호성을 직접 확인 가능 |
| 구면·A₂만 있는 면 | 유한 광선의 근축 극한, vertex curvature가 c+2A₂임을 확인; c=0이 항상 plane이 아님 |
| A₂…A₁₂ 비구면 | sag/gradient/Hessian, 전 항의 thermal similarity, 실제 6면 경로의 역추적·회전·길이/파장 scaling |
| 두 양의 교점·접선 | 최초 hit 선택; 앞선 접선·미해결 구간을 건너뛰어 뒤 root를 답으로 취급하지 않음 |
| 내부에서 교차하는 면 | center와 rim 두께만 양수여도 내부가 겹칠 수 있는 기존 반례 |
| 원형·직사각·annular pupil | Airy/circular MTF, exact discrete triangular OTF, 실제 주파수 축, clipping support |
| tilt·defocus·이중상 | OTF phase 부호, PSF 이동, negative OTF, MTF 크기 평균의 잘못, 서로 다른 sampling/padding |
| 표본 격자와 동기인 진동 phase | 중첩 격자의 거짓 수렴을 비정합 표본과 비교; 격자 증가 하나를 보증으로 삼지 않음 |
| Lambertian aperture | 독립 projected solid angle·해석적 수집 Q, fill·면적·n·cos 전력/진폭 관계 |
| 합성 3매 렌즈 | 실제 6면→STOP→OPL→FFT→spot/RI/영상/열을 이어 실행하는 회귀 fixture |
| 유한 물체 collimator | 출력 angular 분석이어도 광학계 EFL은 유한할 수 있음; C=0 오판 금지 |

기본 3매 fixture의 별도 reduced-ABCD 계산은 EFL **84.1490162331 mm**, 마지막 면 기준 BFL **79.9243093735 mm**다. 실제 유한 aperture 광선의 spot이 0이라고 요구하지 않는다. fixture는 unoptimized synthetic glass이며 사용자의 최종 FOV나 화질 목표를 달성한 설계가 아니다. 상세 선정 논의와 에이전트 간 질문·수용·반론은 `reviews/sample_plan.md`에 있다.

## 독립 리뷰에서 발견하여 반영한 사항

| 발견 내용 | 수정·검증 |
|---|---|
| 정확한 근축 초점에서 position Jacobian이 거의 0이면 field curvature 상대오차가 발산 | 물리 길이 scale을 포함한 오차·coupling 검사, 실제 triplet 초점 회귀 |
| STOP aiming 중 STOP 자체의 TIR를 aiming 실패로 분류 | STOP 교점 잔차와 STOP 굴절을 분리; 마지막 real trace가 TIR 기록 |
| 접선 부근의 작은 residual을 실제 root로 채택 | 첫 구간·bracket·단순근 조건 강화; tangent와 뒤 root 반례 |
| fractional 반복 횟수/STOP index, 극단 ray vector, caller dictionary 변경 | 생성 시 검증, overflow-safe 정규화, read-only ray·고정 설정 snapshot |
| JSON thermal lens group이 dict 그대로 전달 | typed group 변환과 입력 key 검사 |
| 사용하지 않는 'air' material 추가가 reported EFL을 바꿈 | `report_reference_material`을 명시적 설정으로 분리 |
| CTE 평가가 물성 유효 온도 범위를 무시 | expansion에서도 범위 검사; 범위 밖 extrapolation 금지 |
| focal snapshot을 afocal에 재사용하면서 조명 재설정 | mode/reference 재사용 조건 검사 |
| FFT padding 부족인데 phase 검사만 통과 | linear autocorrelation 지원 여부를 validity 상태에 반영 |
| 유한각 T/S에 단순 atan2 좌표 회전 사용 | 실제 3D tangent-frame projection으로 angular fan 계산 |
| Gauss 내부 표본의 hull로 actual STOP을 축소 | rim/chief를 topology용 0 가중치 표본으로 추가; 독립 Lambertian Q 회귀로 확인 |
| aiming과 aperture 허용오차 차이로 실제 STOP rim 광선이 잘못 차단 | aiming 종료에도 최종 trace의 반경 조건 적용; 구경 변경 없이 독립 collimator 사례 397/397 광선 통과 |
| 광원 PSF에 uniform lens pupil을 재사용할 위험 | source projected solid angle·mapped area·Fresnel/absorption으로 physical E-field 구성 |
| Q·h 정규화만으로 광량 맞음을 주장할 위험 | 정규화 전 ray/boundary/spectral power와 원래 phase sampling 상태를 따로 보존 |
| 흡수 매질을 무한 전파하는 far-field | 출력 매질 흡수가 있으면 지원되지 않는 모델로 거부; lossless 선택은 명시 |
| 유한 OBJ에서 angle field를 무한 물체로 처리 | 지원하지 않는 공액·field 조합을 입력 오류로 처리 |

리뷰 문서는 `software_review.md`, `diffraction_physics_review.md`, `radiometry_physics_review.md`, `metrics_camera_review.md`, `geometry_thermal_physics_review.md`, `integration_camera_review.md`다. 작성자와 다른 에이전트가 코드를 읽고 별도의 반례 테스트를 작성했다. 이는 전문가의 외부 인증이나 사람의 실제 광학 시험을 뜻하지 않는다.

## 실제 통합 실행 기록

- `results/demo/summary.json`: G 렌즈 MTF·spot, 20/60/85°C 고정 IMA와 공통 refocus, nominal mapping 변화, uniform-scene RI, finite-source far-field, 기하 bitmap 영상.
- `results/demo/wave_*.npz`: 복소 pupil·OPD·PSF·OTF·MTF와 동일 pupil diffraction reference.
- `results/source_wave_demo.json`: 유한 source cell들의 물리 조명 회절 kernel 및 합성 영상, afocal source cell의 global angular intensity. 비축 detector 투영 근사를 명시한다.
- `reviews/integration_referee_results.json`: 별도 검토자의 ring/arm/grid 독립 refinement, 강한 비축/열 defocus, pure-asphere, 직접 RS/local FFT 비교. `*_pre_rim_fix.json`이 있으면 **수정 전 기록**이며 최종 예측 데이터가 아니다.
- `reviews/source_wave_camera_check.json`: actual triplet의 source Q와 **정규화 전** 전력 비교.

물성·thermal stress 크기가 다른 통합 fixture가 있으므로 두 보고서의 온도 수치를 동일 처방 결과처럼 비교하면 안 된다. 실행시간은 병렬 작업과 같은 환경의 영향을 받으며 제품 성능 benchmark가 아니다.

## 데이터로 채택할 때 남는 조건

강한 field에서 실제 추적이 모두 성공해도 wave phase가 `underresolved`일 수 있다. 더 큰 field에서는 `aim_unresolved`가 발생하며 이 경우 통계/회절 결과를 정상 vignetting으로 계산하지 않는다. 고정 IMA의 큰 thermal defocus에서도 같은 구분을 유지한다. 실제 상태는 기록된 JSON을 따른다.

Source-wave focal 예제에서 normalized PSF가 계산되더라도, 별도 physical-plane raw angular-spectrum의 phase 표본은 `underresolved`일 수 있다. 그 power ratio를 검증 완료 값으로 쓰지 않는다. 작은 phase step 통과도 충분조건은 아니다.

현재 FFT는 명세의 독립 scalar profile과 sampled triangulation을 구현한다. 모든 field·NA·pupil distortion에서 정확한 전파, 실제 Zemax 내부 cosine-space와 동일한 수치 결과를 보증하지 않는다. 정확한 비축 global detector diffraction, 임의 단차 rim, 전체 AGF thermal catalog 호환, vector diffraction는 지원 범위 밖이다. 코드가 반환하는 모델·가정·unit·상태를 데이터와 함께 보존해야 한다.
