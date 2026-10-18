# Sequential Optics Core v0.1.0 — 2차 독립 코드 리뷰

검토일: 2026-10-18. 대상 ZIP의 SHA-256은 `d854cbaba525581088fe7b890ee54e2f1fb663c8ef8c6fb562b06a66a65284ef`다. 배포본을 다시 확보하고 manifest와 대조한 고정 사본을 사용했다. **엔진 소스는 수정하지 않았다.**

## 판정

**기본 광학 계산의 기반은 있으나, 현재 버전을 자동 데이터 생성용으로 포괄 승인할 수 없다. 추가 수정이 필요하다.** 실광선 Snell, 두꺼운 렌즈 EFL/BFL, 모든 A₂–A₁₂ 항의 균일 팽창, 선택적 G/열/분광 데이터 처리에 대한 새로운 독립 대조는 통과했다. 반면 원거리 MTF, 광선 캐시 재사용, 비네팅 처리에서 잘못된 값 또는 불필요한 전체 분석 실패를 재현했다.

기존 **183개 테스트는 이번에도 모두 통과**했다. 그 테스트가 검사하지 않은 조건에서 **확인된 결함 9개**가 추가로 발견되었다. 아래 우선순위는 이 리뷰의 데이터 생성 용도 기준이다.

- **P1 4개:** 해당 조건에서 잘못된 해석 데이터가 저장되거나 정상 해석이 막힌다. 해당 경로의 데이터 채택 전에 수정해야 한다.
- **P2 4개:** 특정 지표·사용자 설정·객체 사용 방식에서 계약과 결과가 어긋난다.
- **P3 1개:** 잘못된 입력의 검증/오류 처리 문제다.

이는 모든 기본 G 결과가 틀렸다는 판정도, 새로 나열하지 않은 코드에 결함이 없다는 보증도 아니다. 기본 공개 `evaluate`는 매번 광선을 추적하므로 R2의 stale-bundle 사례는 주로 저수준 API와 캐시 재사용 경로에 해당한다. R1 역시 요청 각주파수에 비해 출력 각도 표본이 부족한 조건에서 발생한다.

## 독립성 및 검증 방법

새 에이전트 3명이 각각 광학 물리학자, 카메라/투사 렌즈 설계 엔지니어, 과학 계산 소프트웨어 개발자 역할을 맡았다. 코드를 먼저 검토하고 기존 리뷰를 정답으로 취급하지 않았다. 이는 역할을 부여한 AI 검토이며 사람 전문가의 외부 인증은 아니다.

각 검토자는 배포 소스를 수정하지 않고 별도 재현 스크립트를 작성했다. 기존 editable 설치를 잘못 읽지 않도록 고정 소스의 import 경로를 강제했다. 주 검토자는 각 역할의 최종 재현 스크립트를 다시 실행하고, 원인·영향·중복을 검토했다. 물리 검토자와 렌즈 검토자는 두 aperture 처방의 물리적 타당성과 가짜 cylinder miss도 상호 확인했다.

추가 재현 스크립트의 assertion은 **현재 버그가 실제로 발생함**을 확인하는 것이다. 이를 올바른 구현에 대한 통과 테스트 수로 합산하지 않았다. 수정 후에는 기대 동작을 검사하는 회귀 테스트로 전환해야 한다.

## R1 — 각주파수 샘플링 부족을 감지하지 못하는 angular MTF

**P1.** `wave.py:204–219 angular_otf`, `analyze_wave`의 afocal 분기 및 상태 판정. `source_wave`의 afocal MTF도 같은 함수에 의존한다.

실제 공기 aperture를 통과한 광선과 매끄러운 Gaussian 경계 진폭에서, 경계장·광선·파장을 유지하고 FFT padding만 4→8→16으로 바꿨다. 약 **4581.818 cycles/rad**에서 MTF가 약 **0.999999466 → 4.19×10⁻⁸ → 2.97×10⁻⁸**로 바뀌었다. 세 계산 모두 경계 위상 변화 0이며 `conditional`을 반환했다. 별도 스크립트에는 기본 0–10000 cycles/rad 주파수 격자에서도 발생하는 사례를 포함했다.

원인은 `Σ pᵢ exp(-2πifαᵢ)` 적분의 요청 주파수에 대한 각도 샘플링 검사가 없다는 것이다. 경계 파면의 위상 샘플링을 검사하는 것만으로 이 출력 적분의 aliasing을 검출할 수 없다. 경계 에너지가 충분히 작다는 검사도 해결책이 아니다.

수정 방향: 요청 최대 주파수에 맞춰 실제 비균일 각도 좌표에서 적분 위상 변화와 수렴을 검사하고, 부족하면 refinement 또는 `underresolved`/명시적 거부를 반환한다. `atan` 변환 뒤 각주파수에 선형 direction-cosine 공간의 hard cutoff를 그대로 적용해서 0으로 덮어쓰면 안 된다.

증거: `physics/repro_angular_otf_alias.py`, 대응 JSON. Gaussian 해석값은 **소각 근사 참고값**이고 exact-angle 정답이라고 주장하지 않는다. 동일 입력의 padding 수렴 비교와 출력 적분 구조가 결함의 직접 근거다.

## R2 — 변경 전 광선을 변경 후 처방의 결과로 승인

**P1.** `models.py RayBundle`, `wave.py:320–332 _fingerprint`, `make_boundary_state`, `analyze_wave`, `source_wave.py:135–161`.

합성 triplet 마지막 면 곡률을 **0.5% 변경**한 뒤 이전 bundle을 넘겨도 거부하지 않는다. 60 lp/mm MTF는 이전 값 **0.508147**을 그대로 반환하지만 새 처방으로 재추적하면 **0.053244**다. 둘 다 `conditional`이다.

기존 `BoundaryState` 재사용 검증에도 별도 누락이 있다. STOP 반경을 4→2 mm로 바꿔도 fingerprint에 STOP index/radius가 없어서 캐시가 승인된다. 캐시 사용 시 MTF60은 **0.508147**, 올바른 재추적은 **0.22660**이다. 물리 광원 회절 API의 첫 면 Snell/OPL 검사는 이후 면의 변경을 검증하지 못한다.

원인: trace 당시 처방 identity가 bundle에 없다. 새 상태를 만들 때 현재 처방과 전달받은 오래된 광선을 함께 hash하여 불일치 자체를 정상 상태로 기록할 수 있다.

수정 방향: 추적 시점의 처방·굴절률·온도·STOP·공액·field·파장 identity를 bundle에 보존하고 모든 소비 API에서 비교한다. IMA만 이동할 때 허용되는 재사용과 재추적이 필요한 변경을 명시적으로 나눈다. 광원 진폭 캐시는 광량 모델의 identity까지 구분한다.

증거: `physics/repro_stale_bundle.py`, 대응 JSON. 이 결함은 회절 근사에 관한 논쟁과 무관한 데이터 일관성 문제다.

## R3 — STOP 앞의 정상 비네팅이 aiming 실패로 바뀜

**P1.** `geometry.py aim_ray.residual`의 `ignore_apertures=True` 추적과 `_bound_interval`; `metrics.py _valid_bundle` 및 wave 진입 검사.

공기 aperture 두 개를 z=0,10 mm에 두고 둘 다 반경 2 mm, 두 번째 면을 STOP으로 지정한다. 무한 물체 field는 X=0.1 rad, IMA는 z=20 mm다. 서로 어긋난 두 aperture의 통과 영역은 독립 기하로 직접 계산할 수 있다.

rim 표본을 꺼서 R4를 분리하면 192개 중 134개가 정상 통과하고 58개는 물리적으로 잘려야 한다. 현재 코드는 그 58개를 `aim_unresolved`로 반환하므로 spot 전체가 실패한다. 첫 aperture의 clear radius는 그대로 두고 mechanical radius만 4 mm로 확장한 대조군에서는 58개가 `aperture_blocked`로 분류된다. 대조군 RMS **1.185733581 mm**, 통과 quadrature fraction **0.688257821**은 같은 표본의 독립 해와 일치한다. 연속 aperture의 해석적 겹침 비율 **0.684006291**과 유한 quadrature 값을 혼동하지 않았다.

수정 방향: aiming을 위한 가상 연장/잔차 계산과 실제 물리 경계 차단을 분리한다. 다만 모든 `escaped`를 무조건 물리 차단으로 바꿔 수치 실패를 숨겨서는 안 된다. 도달 가능한 STOP 영역과 실패 원인을 판별할 수 있어야 한다.

증거: `lens_design/reproduce_lens_findings.py`의 LD1.

## R4 — bounding cylinder 반올림으로 정상 STOP rim 교점 소실

**P1.** `geometry.py:175–190 _bound_interval`.

R3와 다른 원인이다. front mechanical radius를 충분히 크게 하고 기본 rim 광선을 포함하면, 정규화 pupil `(약 0,-1)`의 광선 하나가 실패한다. 해석적인 실제 평면 교점은 `(약 0,-2,10) mm`, 반경은 정확히 2 mm이고 `normal·direction≈0.995`로 표면 접선도 아니다.

bounding cylinder 이차식의 판별식이 반올림 때문에 **−1.38778×10⁻¹⁷**이 되어 `disc<0`에서 miss로 처리된다. 결과는 146 정상/70 차단/1 `aim_unresolved`이며 이 한 개의 0 가중치 경계 광선 때문에 spot 전체가 실패한다.

수정 방향: 실제 표면의 접선과 보조 bounding cylinder의 접선을 구분하고, 크기에 맞춘 반올림 오차 한계와 안정적인 이차식/평면 교점 처리를 적용한다. 모든 음수 판별식을 0으로 만드는 무조건적 보정은 잘못된 교점을 만들 수 있다.

증거: 같은 재현 스크립트의 LD2. 이전 리뷰에서 수정한 STOP aiming 종료 허용오차와는 별도 결함이다.

## R5 — 사용자 기준 초점이 LSA의 실제 축 교점을 덮어씀

**P2.** `metrics.py longitudinal_spherical_aberration`, `qi==0` 처리.

실제 paraxial focus **91.924309 mm**에서 기준 초점만 10 mm 뒤인 **101.924309 mm**로 지정하면 q=0의 `axis_crossing_z_mm`까지 101.924309로 바뀌고 LSA가 0이 된다. q→0의 인접 실제 광선은 원래 초점으로 수렴하므로 여기서 기대하는 LSA는 **−10 mm**다.

수정 방향: 실제 paraxial limit와 값을 뺄 기준면을 별도 변수로 둔다. 기준면 변경은 LSA offset만 바꿔야 한다. 증거: LD3.

## R6 — 최대 spot 반경에서 이미 추적한 rim이 제외됨

**P2.** `metrics.py:45–56 _valid_bundle`, `_moments`, `spot_diagram`.

`w>0` mask를 통계·좌표·최대 반경에 공용으로 적용한다. 적분용 가중치가 0인 rim/chief는 성공적으로 추적되었어도 `geo_radius_sampled`와 plot 좌표에서 빠진다.

기본 8 rings/24 arms triplet에서 반환 반경은 **0.008772395 mm**, 동일 기준점으로 모든 성공 표본을 사용한 반경은 **0.009032284 mm**로 약 **2.877%** 작다. 1 ring의 큰 차이는 민감도 대조일 뿐 기본 성능 오차로 확대 해석하지 않는다.

수정 방향: centroid/RMS의 quadrature mask와 plot/최대 반경의 기하 표본 mask를 분리한다. **현재 weighted RMS를 rim에 임의 양의 가중치를 줘서 바꾸면 안 된다.** 연속 pupil의 진짜 최대값까지 보증하려면 별도 refinement/극값 탐색이 필요하다. 증거: LD4.

## R7 — 고수준 API에서 설정 전달 누락

**P2.** `api.py:117,139`.

spot과 common refocus 호출에 `settings=settings`가 없다. 사용자 `direction_tol=10⁻¹⁵` 설정에서 16/16 광선 추적과 직접 spot 호출이 성공하는 grazing 수치 경계 사례가, `evaluate`에서는 기본 `10⁻¹²`로 되돌아가 `Ray is parallel`로 실패한다. 결과에 기록된 설정과 실제 사용한 판정 설정이 다르다.

일반 광축 부근 처방의 모든 spot이 잘못되었다는 주장은 아니다. 그러나 사용자 설정을 별도로 관리한다는 핵심 API 계약의 위반이다. 두 호출에 설정을 전달하고 저수준/고수준 결과 동등성 테스트를 추가한다. 증거: `software/reproduce_software_findings.py`의 SW-01.

## R8 — 영상 객체의 생성 후 배열 변경이 전력/면적 조건 우회

**P2.** `imaging.py PixelGrid`, `SourceKernel`, `synthesize_image`.

frozen dataclass 내부 NumPy 배열은 writable이다. **호출자가 배열을 in-place 수정하는 경우** 생성 시 검사를 우회한다. kernel mass를 1→2로 바꾸면 수집 전력 1/검출 전력 2를 반환하고, edge 순서를 뒤집으면 음수 irradiance가 나온다. 정상 실행 중 저절로 광량이 늘어난다는 주장은 아니다.

수정 방향: 방어적 copy의 write flag를 잠그고, 외부 커널을 받는 계산 경계에서도 필요한 불변식을 검사한다. `max(0,total-inside)`만으로 전력 위반을 숨기지 않도록 진단한다. 증거: SW-02.

## R9 — thermal group index 입력형 검증 누락

**P3.** `thermal.py:30–36,71`.

`front_index=0.0, back_index=1.0`인 입력이 생성 단계를 통과한 뒤 raw `TypeError`를 일으켜 `fail_fast=False` batch도 중단한다. bool index도 정수처럼 통과한다. 이것은 **잘못된 입력형 검증 문제**이며 열 물성이 생략되었을 때의 선택적 데이터 처리 문제와 구별한다.

수정 방향: 정수 index 및 bool 배제 정책을 적용하거나 허용할 숫자형을 명시적으로 정규화한다. 모든 예외를 catch해 물리 실패로 바꾸는 방식은 피한다. 증거: SW-03.

## 새 검토에서 확인한 정상 동작

| 대조 항목 | 관찰 |
|---|---|
| 두꺼운 렌즈 독립 해 | EFL 50.8474576271 mm, BFL 49.1525423729 mm 일치 |
| Vector Snell | 예상/실제 transmitted sin이 약 0.197013471108로 일치 |
| 굴절률 고정 균일 팽창 | 길이 배율 1.01에 대해 EFL 비율 1.01 |
| 모든 A₂–A₁₂의 thermal similarity | sag 항등식 오차 약 2.8×10⁻¹⁷ mm |
| nominal G만 있는 물성 | CTE·dn/dT·흡수 없이 first-order/spot 정상 |
| B 및 고온 물성 누락 | G/기준 상태 결과 보존, 해당 요청만 missing_data |
| 기존 검증 | 183개 전체 통과; 이번 새 반례를 커버하지 못한 사실과 양립 |

검토한 파면/광량 식의 범위와 아직 보증하지 않는 사항은 `physics/review.md`에 별도로 있다. 기본식에서 새 결함을 찾지 못했다는 것과 모든 실제 처방의 물리 정확도를 검증했다는 것은 다르다.

## 구현 결함과 구별할 현재의 모델 범위

현재 FFT는 자체 scalar 경계장 모델이다. Zemax와의 수치 일치, 정확한 비축 global detector 전파, vector diffraction, 임의 무효부 단차를 구현했다고 해석할 수 없다. README에 이미 선언된 범위를 이번에 새로 발견한 버그로 세지 않았다.

Ansys의 공개 문서도 FFT의 scalar/Fraunhofer 적용 조건과 pupil distortion·계산면 문제를 별도로 다룬다. 따라서 이 리뷰에서 발견한 버그를 수정하는 일과, 사용자 처방에서 해당 회절 모델이 충분히 정확한지 평가하는 일은 별도 검증 단계다. 특히 `conditional`은 데이터 정확성 승인 표시가 아니다.

- [Ansys — Methods for analyzing MTF in OpticStudio](https://optics.ansys.com/hc/en-us/articles/42661772706835-Methods-for-analyzing-MTF-in-OpticStudio)
- [Ansys — FFT와 Huygens PSF의 차이 및 적용 조건](https://optics.ansys.com/hc/en-us/articles/42661795535251-What-is-the-difference-between-the-FFT-and-Huygens-PSF)
- [Ansys — Tilted image surfaces에서 FFT/Huygens MTF 차이](https://optics.ansys.com/hc/en-us/articles/42661833682707-Why-are-FFT-and-Huygens-MTF-results-different-on-tilted-image-surfaces)

## 수정 순서와 재검증 조건

1. **추적 provenance와 angular MTF 주파수 검증:** 틀린 숫자가 정상 실행 결과로 저장되는 경로를 먼저 차단한다.
2. **비네팅/교점 분류:** 광선의 물리적 차단, bounding 계산의 반올림, 실제 수치 실패를 분리한다. `aim_unresolved`를 일괄 무시하지 않는다.
3. **LSA·최대 spot·설정 전달:** 기준면 불변성과 고수준/저수준 동등성을 검사한다.
4. **객체·입력 불변식:** 생성 후 수정과 malformed thermal 입력 경로를 명확히 한다.
5. 위 반례들을 기대 동작의 회귀 테스트로 전환하고 기존 183개를 함께 실행한다. 이후 대표 실제 처방에서 ray density, FFT grid, padding, 요청 주파수 범위를 독립적으로 바꾸어 수렴을 판단한다.

이 산출물은 **리뷰와 재현 증거**다. 수정된 엔진 또는 새 검증 완료 버전으로 표시하지 않았다.
