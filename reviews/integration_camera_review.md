# 독립 통합 코드 검토 — 카메라·프로젝터 광학 / 시험 설계

검토일: 2026-10-18. 검토자는 광선·회절·재료·metrics·radiometry production 모듈을 작성하지 않았다. 직접 작성한 것은 표본 검토, 독립 시험과 이 통합 실행 기록이다. 기준은 `docs/calculation_spec_ko.md`, 기존 여러 역할의 반론, 실제 현재 코드다. 실제 Zemax·Code V 결과나 실측 정답은 사용하지 않았다.

**판정:** 실행한 범위에서 실제 여섯 굴절면의 ray→spot→OPD→FFT, finite-source angular 경로, 균일 thermal→fixed/common-focus의 연결이 작동한다. Nominal G-only 계산은 RGB·열·광원·흡수 자료 없이 가능하다. 실패와 미수렴을 수치 정상 결과와 구분하는 방향도 적절하다. 다만 수치 결과가 존재한다는 사실과 물리적 적용성/수렴/제품 일치를 구분해야 한다. 아래의 고시야·큰 defocus 예는 실제로 `aim_unresolved` 또는 `underresolved`를 반환하며 이 상태를 제거해서는 안 된다.

원시 요약: `integration_referee_results.json`. 재현: project root에서 `PYTHONPATH=src python reviews/run_integration_referee.py`. rim/chief 수정 후 전체18개 case 재실행은269.52 s이며 direct-RS·pure-asphere와 chief-map 진단을 포함한다. 그 다음 발견된 정확한 STOP 경계의 종료 predicate 수정은 별도16개 독립 회귀 시험과 source x=3.2 mm 케이스 재실행으로 확인했고 해당 행을 갱신했다. 나머지 diagnostic snapshot의 원시 수치는 rim 수정 후 전체 실행 기록으로 보존한다. 실행 시간은 성능 benchmark가 아니다. 독립 회귀 시험은 `PYTHONPATH=src python -m unittest tests/test_adversarial.py tests/test_wave_adversarial.py -v`: **16건 통과, 8.814 s**. 최종 전체 suite의 권위 있는 수와 실행 시점은 root 통합 보고서를 따른다.

## 1. 검토 도중 발견해 수정된 문제

| 문제 | 영향 | 조치·확인 |
|---|---|---|
| Afocal batch의 `mtf` 요청이 angular frequencies 없이 wave를 호출 | 반환은 성공인데 angular PSF만 있고 MTF가 없을 수 있었음 | 독립 angular 주파수 설정을 추가하고 API가 전달. 사용자 설정 400 cycles/rad,9표본이 실제 `mtf_s`에 반영되는 독립 회귀 시험 통과. 120 lp/mm를 전용하지 않음. |
| Afocal field-curvature dispatch가 image_z=None인 focal 경로 사용 | finite projector의 미분 vergence 진단이 불필요하게 실패 | 마지막 면 뒤 유한 진단 reference plane을 API가 구성. `converging_positive_independent` 결과를 반환하는 회귀 시험 통과. 제품의 diopter convention 인증은 아님. |
| Focal boundary state를 afocal로 재사용하면 조명이 조용히 바뀔 수 있음 | 고정 AΠ를 보존한다는 계약 위반 | mode 변경 재사용을 거부. chief와 material index·object plane도 fingerprint에 포함하도록 수정됨. wave 작성자의 해당 회귀 시험과 코드 변경을 확인. |
| Gauss STOP quadrature에 실제 rim 표본이 없어 wave/source triangle hull이 구경을 작게 근사 | 낮은 radial order에서 pupil support와 수집 power의 체계적 감소; source analytic Q 시험에서 발견 | 경계와 chief ray를0 weight로 추가하는 설정을 기본 활성화. 기하 quadrature는 유지하고 실제 pupil edge를 표현. 담당 작성자 수정 후 이 문서의 통합 수치를 다시 실행하여 갱신. |
| 정확한 STOP-rim target인데 aim tolerance가 aperture predicate보다 느슨하여 일부 rim ray가 차단됨 | finite collimator source x=3.2 mm의 rim5개가 false clipping되어 wave support를 왜곡 | aiming 종료 시 잔차뿐 아니라 실제 STOP hit가 tracing과 동일한 radial predicate를 통과하도록 수정. aperture 확대·target 이동 없음. 독립 실패 회귀 시험이 통과했고 해당397-ray 케이스를 재계산하여397/397을 확인. |
| 기하 T/S와 FFT T/S 이름이 서로 다른 축을 뜻할 수 있음 | field가 X쪽일 때 같은 T label을 같은 물리 단면으로 비교할 위험 | FFT 결과에 S=e₁(투영된 X),T=e₂라는 별도 `axis_convention`을 추가. 기하 T는 field의 radial direction인 점을 유지. |

별도 software review에서 수행한 입력 불변성·settings 검증 등은 중복 검토인 것처럼 세지 않는다. 이 검토의 production 수정은 각 소유 작성자가 수행했다.

## 2. 합성 triplet의 실제 계산과 수렴

처방과 독립 ABCD 정답은 `sample_plan.md` §2 및 `tests/test_adversarial.py`에 있다. n=1.50/1.62 합성 상수 매질이며 실제 glass grade를 가장하지 않는다. EFL=84.1490162331 mm, BFL=79.9243093735 mm의 **first-order** 값이 독립 행렬과 일치한다. 유한 구경의 real-ray focus를 이 숫자와 같다고 강제하지 않는다.

G=550 nm, field=(0.03,0) rad, fixed global IMA z=91.9243093735 mm의 예. 아래는 **rim/chief0-weight 표본 추가 후 최종 재실행 값**이다.

| STOP ray sample | FFT grid | S MTF@30 lp/mm | T MTF@30 lp/mm | S MTF@120 lp/mm | max traced triangle phase span |
|---|---:|---:|---:|---:|---:|
| 6 rings ×24 arms + rim/chief | 64² | 0.516478 | 0.637534 | 0.057042 | 2.3056 rad |
| 10 rings ×36 arms + rim/chief | 64² | 0.514085 | 0.635052 | 0.053997 | 1.4724 rad |
| 14 rings ×44 arms + rim/chief | 64² | 0.514081 | 0.634622 | 0.054911 | 1.2001 rad |
| 14 rings ×44 arms + rim/chief | 80² | 0.516202 | 0.636966 | 0.055806 | 1.2001 rad |

앞 세 행은 ray/OPD 보간 표본 변화, 마지막은 같은 ray에서 FFT pupil raster만 바꾼 것이다. N=64→80의 비정수 배율 변화도 사용했다. boundary 표본은 모든 해상도에서 pupil rim을 명시적으로 나타내며,0 quadrature weight이므로 geometric moments를 바꾸지 않는다. **완전 수렴 또는 Zemax 오차율을 인증하지 않는다**. 원시 pre-fix 값은 `integration_referee_results_pre_rim_fix.json`에만 보존했다. 이는 구경을 축소하는 sampling defect 수정 전의 기록이므로 현재 성능 값으로 사용하면 안 된다.

동일 두 field(0 및0.03 rad)에 하나의 공통 RMS 초점을 계산하면 z=91.7551612496 mm다. 축상 RMS spot은4.5563→2.2917 µm,0.03 field는8.5512→3.9365 µm로 감소했다. 이 공통 위치의 세 ray 해상도 간 범위는2.56e-12 mm였다. 이는 이 온건한 표본의 geometric quadrature 안정성이다. 각 field의 독립 best focus를 동시에 사용한 결과가 아니며, MTF 기준 best focus와 같다고 주장하지 않는다.

### FFT와 직접 scalar 적분의 국소 비교

같은14×44+rim/chief ray,80²의 **동일 경계장**에서 focal direct RS와 local Fourier integral을 비교했다.

| chief-normal 출력 좌표 x,y mm | `(I_local−I_RS)/I_RS` |
|---|---:|
| 0.000,0.000 | +4.589e-10 |
| 0.005,0.000 | +5.177e-06 |
| 0.020,0.000 | +1.888e-06 |
| 0.040,0.000 | -3.207e-05 |

이 범위에서는 local propagation 근사가 잘 맞는다. 두 계산은 ray/OPD/경계장·pupil raster를 공유하므로 그것들의 독립 검증은 아니다. 전체 FFT window 모서리의 최대 phase expansion 잔차는3.96–6.19 rad이다. finer pupil grid에서는 같은 support에 대해 PSF의 표현 window가 넓어지므로 full-window 최대값이 오히려 커질 수 있다. 이를 전체 MTF 오차 한계나 중심 PSF 오차로 바꾸어 해석하지 않는다. 모든 위 wave 결과는 `conditional`이며 local 비교와 별도로 수렴·edge energy를 검토해야 한다.

## 3. 비구면·큰 field·실패 상태

- 제3렌즈 앞면의 A₂…A₁₂를 모두 nonzero로 둔 controlled 표본은289/289 ray가 통과했다. 같은 fixed IMA에서 S MTF@30 lp/mm는0.158555로 바뀌었다. 고차항이 전체 경로에 반영됨을 보여 주는 연결 증거이며 최적화된 처방이 아니다.
- 첫 면을 **c=0,A₂=0.01 mm⁻¹**, 나머지 A₄…A₁₂도 nonzero로 바꾼 pure-asphere 표본은661/661 ray가 통과했다. vertex curvature=0.02 mm⁻¹로 원 구면의 first-order power를 유지한다. RMS spot=8.2746 µm,S MTF@30 lp/mm=0.524077,wave=`conditional`다. c=0을 평면으로 잘못 해석하지 않는다.
- field0.15 rad의 ray 상태: `{'ok': 289}`. wave=`underresolved`,triangle phase span=13.535 rad. 광선 추적 성공은 회절 표본 충분성의 증거가 아니다.
- field0.45 rad의 ray 상태: `{'escaped': 30, 'aim_unresolved': 16, 'ok': 243}`. aggregate 계산을 `NumericalError`로 거부했다. unresolved를 aperture loss로 넣지 않는다.
- field0.6 rad의 ray 상태: `{'aperture_blocked': 7, 'escaped': 62, 'aim_unresolved': 49, 'ok': 171}`. aggregate 계산을 `NumericalError`로 거부했다. unresolved를 aperture loss로 넣지 않는다.
- 큰 field의 단일-branch aiming이 임의의 유효 ray를 모두 찾는다고 보장할 수 없다. 실제 실패 상태를 보존하는 것이 적절하다.
- 독립 반례 `z=r²`의 첫 양의 교점 선택, 중심·rim은 양수지만 중간에서 겹치는 gap=`1−8ρ²+8ρ⁴` 거부 시험이 통과했다.

## 4. Finite-source projector와 angular 출력

OBJ z=−78.5212427950 mm의 paraxial collimator 공액을 사용했다. 시스템 자체의 C는0이 아니며 EFL은 유한하다. "finite object의 출력이 평행"과 "zero-power afocal lens train"을 혼동하지 않는 독립 시험도 통과했다.

| Source x mm | 성공 ray(rim/chief 포함) | Angular RMS spot rad | Wave 상태 | S angular MTF@200 cycles/rad |
|---:|---:|---:|---|---:|
| 0 | 397 | 7.603005e-05 | conditional | 0.981692 |
| 3.2 | 397 | 1.599714e-04 | conditional | 0.971134 |
| 6.4 | 397 | 3.934212e-04 | underresolved | 0.899757* |

`*` underresolved 수치는 진단용이며 확정 성능 데이터로 채택하지 않는다. 6.4 mm chief는 `(−0.07582434,0,0.99712119)` 방향이다. 이 표본은 수평 FOV 약8.70°이고 사용자의 최종24° 목표 렌즈가 아니다.

12.8×3.2 mm, fill0.64의 **연속 Lambertian** source를1024개의4D Sobol ray로 entrance-plane importance sampling한 geometric far-field도 실행했다. emitted relative power=82.3549665, exit relative power=0.1635172, direct efficiency=0.00198552였다. 절대 radiance를 주지 않았으므로 이 전력을 실측 W라고 표시하면 안 된다. status는378 ok,209 aperture_blocked,437 escaped이며 보고 angular ROI 밖의 전달량은0이었다. 한 sample count의 집계값이지 Imax 또는 효율의 완전 수렴 인증은 아니다.

launch plane−0.01 mm는 보수적인 첫 표면 sag lower bound 앞이라는 guard를 통과하지 못했다. −0.5 mm로 바꾸어 실행했다. 실제 표면에 대한 가장 타이트한 bound는 아니지만 importance PDF가 optical surface를 건너뛰는 위험을 막는 보수적 조건이다. 검토 중 이 입력 거부를 solver failure나 광속 손실로 오인하지 않았다.

## 5. 온도, 공통 초점과 고정 sensor

시험용 양의 재료 CTE=70×10⁻⁶/K,dn/dT=−10⁻⁴/K; 음의 재료 CTE=8×10⁻⁶/K,dn/dT=3×10⁻⁶/K를 **명시적으로 추가**했다. 기준 매질의 실제 grade 또는 실측 thermal law가 아니다. 공기는 명시적 temperature-invariant index 모델이다. default fixed anchors, fixed global IMA, 표면을 따르는 STOP을 사용했다.

| T °C | 공통 RMS focus z mm | field0.03 fixed RMS µm | common-focus RMS µm | fixed wave | common wave |
|---:|---:|---:|---:|---|---|
| 20 | 91.755161 | 8.5512 | 3.9365 | conditional | conditional |
| 80 | 96.573252 | 155.5433 | 3.6296 | underresolved | conditional |
| 0 | 90.256440 | 62.5910 | 4.0523 | underresolved | conditional |
| 20 재호출 | 91.755161 | 8.5512 | 3.9365 | conditional | conditional |

모든 T는 동일 T₀ 처방에서 독립 생성되었고 T₀ 재호출은 수치상 원래 결과와 같았다. fixed IMA는 항상91.9243093735 mm다. hot/cold fixed-plane MTF가 underresolved이므로 그 절댓값으로 열성능을 확정하지 않는다. common refocus는 명시적 보상량이며 actual fixed detector가 저절로 이동한 것이 아니다. 같은 boundary state의 물리 AΠ를 유지하고 reference만 바꾸어 refocused FFT를 계산했다.

온도별 chief coordinates도 JSON에 별도로 남긴다. 각 온도의 best-fit focal length로 왜곡을 다시 정규화하여 nominal scale/FOV drift를 지우지 않는다. 길이와 파장을 동시에 scale하는 독립 시험은 단위·Fourier 일관성용이고 이 실제 열 시나리오와 구분했다.

## 6. 의존성과 재사용성

- G를532 nm로 바꾼 settings와 **단 한 개의532 nm n 표본**만으로 actual triplet spot와 FFT batch가 실행되었다. 자료에 CTE,dn/dT,absorption,R/B 또는 source pixels는 없다. 분석에 필요 없는 정보를 요구하지 않는 계약을 실제 API 경로에서 확인했다.
- 같은 재료로650 nm 추적을 요청하면 MissingDataError다. single-G 표를 임의의 achromatic material로 승격하지 않는다. 의도적인 constant-n 자료는 별도 모델로 허용한다.
- T=T₀의 geometry no-op에는 CTE가 필요 없다. 다른 T의 geometry에는 CTE, 그 T의 굴절률 계산에는 n(T) 자료가 각각 필요하다. geometry-only와 index-only 의존성이 분리되어 있다.
- Source pixel 개수·fill factor는 lens-only MTF의 함수 인수가 아니다. 별도 source 분석에만 들어간다. 이 구조와 G-only FFT 호출은 적절하다.
- 전역 z translation, 회전된 field+pupil, 전체 길이+λ similarity를 실제 여섯 면으로 시험했다. 좌표 원점 변화의 공통 piston을 OPL 차이에서 구분했다.
- low-level FFT는 이미 wrapped된 복소 pupil만으로 숨은 고주파 phase를 알아낼 수 없다.64cycles ramp는 N=32와64에서 모두 잘못된0 shift로 보이며 N=70에서 불안정성, N=160에서 올바른 shift가 드러났다. 이 시험은 alias를 자동 해결했다는 뜻이 아니라 "두 해상도 일치만으로 충분하지 않다"는 반례를 보존한다.

## 7. 남는 범위·적용성 제한

1. **기본 focal FFT 출력면:** 현재 chief-normal virtual plane의 독립 scalar profile이다. 기하 spot의 global-Z IMA와 비축에서는 다른 평면이다. 축별 label과 output basis를 확인해야 한다. 정확한 tilted global detector diffraction 또는 Zemax의 모든 cosine-space 옵션과 동일하다고 할 수 없다.
2. **Ray→wave 표본:** 현재 conservative transmitted STOP triangles와 piecewise-linear unwrapped OPD를 사용한다. 명세가 설명한 uniform exit-grid 직접 aiming의 완전한 구현은 아니다. sampled topology/phase 진단은 전역 caustic 부재의 수학적 증명이 아니다.
3. **광원 회절 연결:** 완성된 `source_wave.py`의 actual triplet finite-source afocal 경로를 별도로 실행했다. on-axis source cell0.001 mm²,시험 radiance1 W/(mm² sr),fill0.64에서 Q=0.003976102466 W/mm², cell power=3.976102466×10⁻⁶ W였다. global angular bin 합은 같은 cell power이며 이는 Q/h 정규화·3D 등록 검사다. **정규화 전** boundary grid/ray power=0.99158445, wave/boundary=0.999999984였다. 전자의 약0.84% 차이는 지워지지 않고 raster 오차 진단으로 남는다. raw phase step=0.282 rad,triangle span=0.709 rad,status=`conditional`이다. 원시 기록은 `source_wave_camera_check.json`; 독립 actual-triplet 회귀 시험도 추가했다. on-axis detector 변환은 지원하지만 비축 global-Z detector 변환은 local chief-parallel projection을 명시적으로 opt-in해야 한다. 이것은 정확한 비축 detector diffraction의 구현이 아니다. 상세 독립 수식 검토는 geometry/wave 담당 검토와 최종 suite를 함께 따른다.
4. **기계 rim:** active aperture는 광선 전달을 제한한다. mechanical diameter와 무효부 높이는 전체 rim CAD 형상·접촉·응력을 유일하게 정하지 않는다. metadata를 보존하더라도 mechanical collision/thermal mount 해석을 구현했다고 표시하지 않는다.
5. **물성·공정:** 합성 물성 검증은 실제 N-BK7/PMMA/PC grade, 분산/열/투과 측정값의 검증을 대신하지 않는다. 제조 공차, stress/gradient/GRIN, scattering/ghost는 사용자 제외 범위다.
6. **수렴 완료 판정:** 이번 표본은 representative integration이다. 개별 설계를 데이터셋으로 채택할 때는 중요한 metric마다 별도의 수렴·aperture·throughput·모델 적용성 기록이 필요하다. MTF 한 숫자나 성공한 ray 개수만으로 물리 정확도를 주장하면 안 된다.

수식/과정의 주요 위험을 실제 반례로 검사하고 고쳤다는 점은 긍정적이다. 검토한 실행 범위에서 새로 확인된 미수정의 확정 수식 오류는 없으나, 이 판정은 위 조건부 표본·구현 profile에 한정된다.
