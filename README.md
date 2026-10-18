# Sequential Optics Core

3매 렌즈 기술 명세를 바탕으로 작성한 **Python 독립 광학 계산 엔진**이다. 실제 비구면 교점·벡터 Snell 추적과 명시적인 scalar 회절 profile을 사용한다. 렌즈 수를 3으로 하드코딩하지 않아 같은 모듈을 다른 면 수에도 재사용할 수 있다.

실행 가능한 코드, 자동 테스트, 합성 처방, 실제 실행 결과, 다중 에이전트의 독립 코드 리뷰가 함께 있다. **Zemax 결과와의 수치 일치 또는 측정 정확도를 검증한 제품은 아니다.** `conditional`, `underresolved`, 수치 실패와 물리적 차단은 구분해서 반환한다.

## 현재 수정본 v0.2.1

후속 리뷰의 C1–C4를 수정하고 field별 T/S FFT MTF의 축 metadata·예제·회귀 검증을 추가했다. 수정 및 독립 후속 리뷰는 `FIX_C_REPORT.md`, `VALIDATION.md`, `reviews/postfix_*_c.md`를 참고한다. 이전 R1–R9의 235개 회귀는 유지하며 그 수정 이력은 `FIX9_REPORT.md`에 있다. GitHub/Codex로 이전했으며 시작 방법은 `CODEX_START.md`를 따른다. D1의 성분별 F-θ 정의는 이번에 변경하지 않았으며 `docs/D1_explanation_ko.md`에서 그림과 함께 설명한다.

핵심 원리와 수식은 `docs/core_principles_ko.md`를 먼저 읽는다. 같은 내용의 `docs/core_principles_ko.tex`도 제공한다. 긴 기존 계산 명세는 상세 알고리즘·예외·논점을 확인할 때 사용한다.

## Codex에서 이어서 작업

`CODEX_START.md`에 설치·검증·시작 프롬프트와 미해결 사항을 정리했다. 소스와 테스트는 저장소 루트 기준이다. 큰 재생성 결과 23개는 원본 v0.2.1 ZIP에 보존하며 Git에서는 생략했다. 대상 목록·해시는 `MIGRATION_MANIFEST.json`에 있다.

## 실행

Python 3.10 이상, NumPy 1.24 이상, SciPy 1.10 이상이 필요하다. 프로젝트 폴더에서:

```bash
python -m pip install -e .
python run_tests.py
python examples/run_demo.py
python examples/run_source_wave.py
python examples/run_field_mtf.py
python -m sequential_optics examples/triplet_G_only.json --settings examples/settings.json --output result.json
```

설치하지 않을 때는 `src`를 `PYTHONPATH`에 추가한다. 테스트는 표준 `unittest`를 사용하므로 pytest가 필요 없다. `run_tests.py`는 실제 실행 환경·개별 테스트·오류·시간을 `results/test_results.json`에 기록한다. 예제는 `results/demo`에 JSON 및 복소 OTF·PSF·pupil 등을 담은 NPZ를 만든다.

```python
from sequential_optics import EngineSettings, RunRequest, Field, evaluate
from sequential_optics.io import load_system

system = load_system('examples/triplet_G_only.json')
settings = EngineSettings.from_json('examples/settings.json')
request = RunRequest(fields=(Field(0, 0), Field(0, 0.03)),
                     metrics=('first_order', 'spot', 'mtf', 'common_focus'))
result = evaluate(system, request, settings)
```

`Field(kind='angle')`는 무한 물체의 rad 단위 X/Y 각도다. 유한 물체는 `object_z`와 `Field(x, y, kind='height')`를 사용한다. `height`의 단위는 mm다. 유한 OBJ와 angle field를 함께 넣어 무한 물체로 조용히 바꾸지 않는다.

## 설정과 선택적 데이터

기본값은 `settings.py`의 dataclass와 별도 `examples/settings.json`에서 관리한다. 사용자 설정을 바꾸려면 JSON 또는 `dataclasses.replace`로 새 설정을 만든다. 입력 dictionary는 방어적으로 복사·고정한다.

| 항목 | 동작 |
|---|---|
| 기본 파장 | `primary_wavelength='G'`; 배포 예제의 G는 **550 nm 진공 파장**이다. 사용자의 실제 G 값으로 변경해야 한다. R650/B450 역시 변경 가능한 예제 기본값이다. |
| 단색 n 하나 | `Material(..., wavelengths_nm=(550,), indices=(1.5,))`. G만 계산하면 R/B, Abbe, dn/dT, CTE, 투과율이 필요 없다. |
| 분산 없는 모델 | `n_constant=1.5`는 사용자가 **명시적으로 선택한 무분산 모델**이다. 단일 측정값을 전 파장 데이터로 자동 확장하는 기능이 아니다. |
| 누락 파장 요청 | 해당 파장만 `MissingDataError`/`missing_data`. 기본적으로 다른 파장 결과는 유지한다. 외삽하지 않는다. |
| 열 데이터 없음 | 기준 자료의 nominal 추적·MTF·spot은 가능하다. 온도를 실제로 바꾸면 필요한 n(T)/CTE를 요구한다. |
| 열 형상만 계산 | `thermal_system`은 형상 팽창용 CTE만 요구한다. 이후 광선 계산 시 n(T)가 필요하다. |
| 흡수 데이터 없음 | 렌즈 단독 spot/MTF는 영향 없다. 광량 API에서는 `missing_absorption='transparent'` 가정을 결과에 기록하며, `'error'`로 엄격하게 요구할 수 있다. |
| 절대 광원 밝기 없음 | RI·정규화 광도·상대 영상 계산 가능. 결과 단위를 relative로 표시한다. 절대 W/sr로 표시하지 않는다. |
| 온도 미지정 | `temperature_c=None`은 각 재료가 제공한 기준 n을 그대로 쓰는 nominal 상태다. 재료 기준온도가 서로 다르면 공통 물리 온도라는 뜻이 아니다. 열 비교는 `temperatures_c=(T0,T1,...)`로 공통 온도를 명시한다. |

길이·OPL·OPD는 mm, 입력 파장은 vacuum nm, 내부 n은 absolute다. mm/μm/nm를 혼용하지 않는다. 렌즈 물성 이름만으로 실제 재료 데이터를 추정하지 않는다. 합성 예제의 crown/flint 수치는 실제 제품 물성표가 아니다.

`SamplingSettings`에는 pupil rings/arms, FFT grid/padding, 경계면 clearance, focal 0–120 lp/mm 범위, 별도 afocal cycles/rad 범위가 있다. `RaySettings`에는 교점·aiming·미분 허용오차와 반복 제한이 있다. Afocal 주파수의 기본 상한 10000 cycles/rad는 독립 설정이며 120 lp/mm의 단위만 바꾼 값이 아니다.

## 모듈과 공개 API

| 모듈 | 주요 기능 |
|---|---|
| `models`, `settings`, `errors` | 공통 자료형, 설정, 누락 데이터/형상/수치/적용성 오류 |
| `materials` | 단일/다파장 n, Sellmeier, T×λ 표, absolute dn/dT, CTE 정의, 흡수·내부 투과율, 데이터 범위 |
| `geometry` | sag/gradient/Hessian, 최초 유효 교점, Snell/TIR, STOP aiming, actual ray bundle, 전체 clear 구간의 면 교차 검사 |
| `paraxial` | A₂를 포함한 vertex curvature, reduced ABCD, EFL/BFL, 공액, EPD/명목 paraxial F/#→물리 STOP |
| `metrics` | spot/angular spot, chief map, distortion, grid·각 점의 spot, focal/angular ray fan, LSA, coma descriptor, 미분적 T/S field curvature, Zernike, 공통 geometric refocus |
| `wave` | canonical 경계면·OPD, FFT PSF/complex OTF/MTF, 동일 pupil의 diffraction reference, angular FFT/MTF, scalar RS 직접 적분 |
| `radiometry` | Fresnel/Jones power·Beer, reverse-ray 일반 RI, 연속 Lambertian source의 far-field, 입체각·전력 ledger |
| `imaging` | 실제 광선 기반 source별 Q·kernel, bitmap 영상, 공통 좌표의 RGB 합성, complex OTF 합성 |
| `source_wave` | 광원 방향의 구면 삼각형 측도→수집 전력→경계장의 전기장 진폭→회절 PSF, raw power 진단, 보존적 image kernel·global angular bin 변환 |
| `thermal` | 기준 처방에서 온도별 독립 생성, 모든 A₂…A₁₂와 크기 팽창, anchor/gap/STOP/OBJ/IMA 정책 |
| `api`, `io` | 선택한 해석만 실행, JSON 처방/결과, Zemax에 익숙한 radius/thickness/diameter 행 입력 변환 |

핵심 API `evaluate`는 first_order/spot/mtf/ray_fans/field_curvature/mapping/common_focus를 묶는다. RI·광원·영상은 서로 다른 물리 입력이 필요하므로 `radiometry`와 `imaging`의 명시적인 API로 호출한다. `examples/run_demo.py`에서 연결 예제를 확인할 수 있다. 제조 공차는 사용자 범위에서 제외되었다.

## 처방 입력

`Surface`의 `z`는 절대 vertex 위치, `curvature=1/R`, `clear_radius`는 active semi-diameter다. `coefficients`는 `(A2,A4,A6,A8,A10,A12)` 순서이며 생략된 뒤 항은 0이다. Standard는 A항을 받지 않는다. 각 면의 `medium_after`가 그 다음 구간의 물성을 정한다.

`from_sequential_rows`는 `radius_mm`, `thickness_mm`, `active_diameter_mm`, `mechanical_diameter_mm`, `conic`, `coefficients`, `surface_type`, `medium_after` 행을 받는다. 앞면 thickness는 중심 두께, 뒷면 thickness는 다음 면까지 거리라는 순차 의미를 그대로 vertex 위치로 변환한다. 마지막 thickness와 별도 `image_z`를 동시에 강제하면 오류다. 실제 .ZMX 파일 parser나 GUI는 포함하지 않는다.

STOP은 surface index 및 물리 반경으로 지정한다. `resolve_aperture`에 EPD 또는 명목 paraxial F/# 하나를 주어 기준 온도에서 변환할 수 있다. 제품의 모든 working-F/#·자동 pupil solve 옵션을 복제하지 않는다. 열 변환 전에 STOP을 resolve하고, 열 정책에서 유지·팽창 여부를 지정한다.

`mechanical_radius` 외의 실제 무효부 단차는 단일 높이만으로 형상을 정할 수 없다. 현재 `inactive_rim_policy='analytic_sag_opaque'`는 sag를 mechanical 반경까지 연장한 **명시적 불투명 경계 모델**이다. `inactive_height_mm`는 입력 보존과 균일 팽창용 metadata이며 단차·edge wall의 실제 교점을 만들지 않는다. 임의 무효부 형상/기구 간섭을 지원한다고 해석하면 안 된다.

## 결과를 해석하는 기준

1. `status='ok'`는 API 실행 성공이다. `analysis_status`와 내부 상태를 함께 본다. FFT의 `conditional`은 표본·scalar 적용성·local propagation 검증이 남는 독립 모델 결과이며 검증 완료 표시가 아니다. `underresolved`는 분석용 정답 데이터로 채택하지 않는다.
2. `aperture_blocked`/`tir`와 `aim_unresolved`/수치 실패를 같은 vignetting으로 처리하지 않는다. 미해결 광선을 0으로 넣어 그럴듯한 MTF를 만들지 않는다.
3. pupil ray 수, wave grid, FFT padding을 각각 바꿔 검토한다. Padding만 늘리면 원래 파면 표본 부족은 해결되지 않는다. 샘플 수렴은 Zemax 호환이나 Maxwell 정확도의 증명이 아니다.
4. 렌즈 단독 FFT는 별도 uniform-Fourier-pupil 규약이다. LED fill/pixel MTF를 곱하지 않는다. 실제 source 영상에는 source에 맞는 Q와 pupil illumination이 필요하다.
5. Focal FFT PSF의 국소 출력은 chief-normal plane이다. 실제 global-z 센서의 기하 spot과 좌표면이 다르다. 비축 detector 회절을 정확히 옮기려면 별도 전파를 해야 한다.
6. FFT의 `T/S`는 fixed chief-normal Y/X 주파수 절단이며, field별 `mtf_t`, `mtf_s` 및 diffraction reference를 반환한다. `axis_convention_id`와 `mtf_frequency_axes_global`, `mtf_target_line_axes_global`에 실제 방향을 기록한다. Zemax FFT의 고정 object-X/Y 이름 규약을 따르는 방향이지만 제품의 임의 2D mapping과 수치 동일성은 미검증이다. 기하 수차의 field-meridional T/S와 구별한다. `docs/field_mtf_usage_ko.md`에 접근 예제와 단위를 설명한다.
7. 고정 영상면과 **모든 field가 공유하는 하나의 refocus**를 구분한다. 파면 focus sweep은 같은 BoundaryState를 재사용해 조명을 유지한다. 매 field 최적값의 envelope를 한 센서 성능으로 저장하지 않는다.
8. MTF magnitude 외에 complex OTF·전체 곡선·광량을 보존한다. RGB는 등록된 intensity PSF 또는 complex OTF를 합친다. MTF 절댓값의 평균은 일반적으로 틀리다.

## 검증 자료와 남는 한계

`docs/calculation_spec_ko.md`에는 식·계산 절차·44개 출처와 기존 미결 쟁점이 있다. `reviews`에는 샘플 선정 논의, 작성자와 다른 에이전트들의 물리/카메라/소프트웨어 코드 리뷰, 실제 통합 실행 결과가 있다. `VALIDATION.md`는 최종 테스트와 수정 사항을 요약한다.

검증에는 plane/plate/Snell/analytic diffraction의 독립 해, first-root 반례, A₂와 모든 고차항, 전체 aperture 내부 교차, reciprocal ray/eikonal, raw power·Jacobian, 공통 초점·RGB 반례, 누락 데이터와 JSON, 실제 6면 처방을 사용한다. 예제는 최적화된 카메라/헤드램프나 실제 N-BK7 제품 검증이 아니다.

주요 제한은 scalar 회절·단일 경로, 회전대칭 동축 면, local regular aiming, sampled pupil topology, 실제 코팅 위상/벡터 회절/ghost·scatter 미지원이다. AGF/제조사 thermal catalog adapter 전체, T-dependent absorption 표, decenter/tilt, 곡면 센서, stress/GRIN/기구 FEA도 포함하지 않는다. 제공한 absolute n(T) 표/선형 모델로 현재 범위의 균일 열 해석은 실행할 수 있다.

기하 far-field와 기본 `geometric_image`는 회절을 포함하지 않는다고 결과에 표시한다. 이 결과를 FFT image simulation으로 명명하지 않는다. 실제 시스템에 사용할 데이터는 해당 처방에서 수렴·모델 적용성을 확인한 뒤 선택해야 한다.

`source_wave.source_diffraction_kernel`은 실제 finite source에 맞는 별도 회절 경로다. STOP 삼각형이 source에서 차지하는 projected solid angle과 투과율로 Q를 구하고, 같은 삼각형이 출사 경계면에서 차지하는 면적으로 전기장 진폭을 만든다. `sqrt(ray.power)`를 바로 FFT pupil에 넣지 않는다. `raw_power_diagnostics`는 PSF 정규화 전 ray/boundary/spectral power 차이와 raw phase underresolution을 보존한다. 정규화 후 Q·h의 전력 일치만으로 진폭 정확도를 주장하지 않는다.

이 경로의 축상 focal PSF는 실제 z 영상면에 보존적으로 resample할 수 있다. 비축에서는 `allow_local_detector_projection=True`를 명시해야 chief-normal PSF의 국소 평행 투영 근사를 사용한다. 정확한 큰 tilt/NA의 센서 회절 전파는 미구현이다. Afocal response는 local 방향을 global 방향으로 3D 회전한 뒤 실제 sr bin으로 변환한다. 여러 source cell의 kernel은 `synthesize_image`, angular bin 전력은 비간섭 합산으로 연결한다. Macro-cell midpoint·삼각형·FFT grid 모두 각각 수렴을 확인해야 한다.

기본 ray bundle에는 Gauss 적분 광선 외에 **STOP rim과 chief를 topology용 0 가중치 광선으로 추가**한다. 기하 quadrature의 가중치 합은 그대로 1이며, 실제 pupil 경계를 Gauss 최외곽 ring으로 축소하지 않는다. `include_pupil_boundary=False`는 연구·비교용 선택이며 그 경우 wave support가 내접 polygon으로 줄어드는 것을 감안해야 한다.

## v0.2.0의 광선·경계장 재사용과 각주파수 판정

`trace_bundle`은 추적 당시 처방·STOP·OBJ·field·파장·온도·실제 사용 n·광선 설정과 전체 광선 내용의 provenance를 저장한다. 변경된 처방에 이전 bundle을 전달하면 소비 API가 거부한다. IMA만 이동하는 refocus와 FFT grid/padding/보고 주파수 변경은 같은 광선으로 계산할 수 있다. 실제 bundle의 pupil 좌표·가중치는 내용 검증에 포함되며, 다음 추적에 쓸 기본 rings/arms 설정과 혼동하지 않는다. `BoundaryState`는 불변 snapshot이며 전체 내용을 검증한다.

직접 구성한 analytic bundle은 `analyze_wave(..., allow_unverified_bundle=True)`와 같이 명시해야 하며 `trace_provenance_status='unverified_external_bundle'`로 기록한다. 물리 광원 회절은 검증된 engine trace만 허용한다. 이 표시는 출처 일치 판정이며 물리 정확성 인증이 아니다. 흡수 자료만 바뀌어 기하 광선이 그대로 유효하더라도 source 진폭과 전력은 현재 자료로 다시 계산한다.

Afocal angular OTF는 실제 angular cell span과 요청 주파수의 적분 위상을 검사한다. `SamplingSettings.angular_otf_phase_step_limit_rad`와 `angular_otf_unresolved_power_tolerance`에서 기준을 관리한다. 저수준 `angular_otf`는 부족한 표본을 거부하며 `return_diagnostics=True`는 해당 주파수에 NaN과 resolved mask를 반환한다. 고수준 lens/source 분석은 `underresolved`를 전달한다. 이 검사를 통과해도 수렴은 별도로 확인해야 한다. Focal MTF의 lp/mm와 angular cycles/rad를 혼용하지 않는다.
