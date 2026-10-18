# 독립 표본·반례·시험 설계 검토

검토 역할: 광학 물리 / 카메라·프로젝터 시험 설계. 작성일: 2026-10-18. 이 문서는 실제 처방의 설계 최적화나 Zemax 수치 정답을 제공하지 않는다. `docs/calculation_spec_ko.md`, `CONTRACT.md`, 기존 광학·회절·radiometry 감사와 세 관점 검토를 바탕으로, 구현 작성자와 구분된 검토자가 시험 표본과 반례를 제안한다. 시험의 실행 상태는 최종 실행 기록과 대조해야 한다. 아래의 수치 기준은 사용자 설계 합격 기준이 아니다.

## 1. 증거를 네 가지로 구분한다

1. **해석식 검증:** 단순화한 문제에 닫힌 해가 있는 경우다. Snell, 평행판, Fresnel, 원형 pupil overlap, affine 영상 power 등은 실제 정답과 비교할 수 있다. "정답지가 없다"는 것은 이 독립 정답까지 없다는 뜻은 아니다.
2. **불변량·상호관계 검증:** 회전 대칭, piston 불변, amplitude² power, 선형 광량, 단위 변환, T₀ 복귀, 길이+파장 similarity 등이다. 공통된 잘못된 구현이 양쪽 결과를 동시에 오염시킬 수 있으므로 독립 해석식과 함께 사용한다.
3. **수치·연결 검증:** 실제 3매 처방의 여섯 굴절면, STOP aiming, OPL, 경계장, FFT, RI, 온도와 공통 초점을 함께 실행한다. 유한 결과, 상태 일관성, 해상도·step 변화의 안정성은 구현이 작동한다는 증거다. 외부 정답 또는 실제 물리 정확도의 증명은 아니다.
4. **모델·제품 대조:** scalar 회절 적용 범위, 실제 재료·실장, Zemax의 특정 release와 결과 일치는 별도다. 이번 synthetic 표본으로 인증하지 않는다.

검사 통과 개수만으로 정확도를 백분율로 제시하지 않는다. 서로 다른 solver의 차이를 자동으로 통계적 신뢰구간으로 부르지 않는다.

## 2. 실제 3매 연결 표본

### 2.1 처방

이름: `synthetic_triplet_v1`. 단위 mm, curvature mm⁻¹. 무코팅·비흡수의 **명시적 상수 n 합성 매질**이며 BK7, crown, flint의 실측 자료를 가장하지 않는다. n이 상수라는 모델 선언과 G 한 점만 알고 있는 측정 표는 별개의 입력이다.

| 면 | z | curvature | conic | clear radius | medium after |
|---|---:|---:|---:|---:|---|
| 0 | 0 | +1/50 | 0 | 6 | synthetic positive n=1.50 |
| 1 | 3 | −1/50 | 0 | 6 | air n=1 |
| 2 | 5 | −1/40 | 0 | 5 | synthetic negative n=1.62 |
| 3 | 7 | +1/40 | 0 | 5 | air n=1 |
| 4 | 9 | +1/45 | 0 | 6 | synthetic positive n=1.50 |
| 5 | 12 | −1/45 | 0 | 6 | air n=1 |

STOP은 면 2, 물리 반경 4 mm. 첫 표본은 Standard sphere만 사용한다. 기계 반경은 clear radius 이상으로 선언하되 기계 rim에서 투과를 허용하지 않는다. 구면 radicand의 여유가 충분하고 각 렌즈의 전·후면과 인접 렌즈 사이 간격은 유효 반경 전체에서 양수다. 양의 단일렌즈 두 장과 음의 한 장으로 굴절 부호, 재질 변경, 공기 gap, STOP 전·후 경로가 모두 포함된다. 그러나 이것은 Cooke triplet의 검증된 처방도 색수차가 보정된 렌즈도 아니다.

독립 reduced-state `(y,nθ)` 행렬 곱:

```python
import numpy as np
z = [0, 3, 5, 7, 9, 12]
c = [1/50, -1/50, -1/40, 1/40, 1/45, -1/45]
n_after = [1.50, 1, 1.62, 1, 1.50, 1]
M = np.eye(2)
n_before = 1.0
for j in range(6):
    if j:
        M = np.array([[1, (z[j]-z[j-1])/n_before], [0, 1]]) @ M
    M = np.array([[1, 0], [-(n_after[j]-n_before)*c[j], 1]]) @ M
    n_before = n_after[j]
A, B, C, D = M.ravel()
print(M, -1/C, -A/C, D/C)
```

결과는 다음과 같다. 이것은 **first-order 정답**이며 유한 구경 real-ray 초점을 같은 숫자라고 강제하지 않는다.

- `M ≈ [[0.94979493, 9.56993756], [-0.01188368, 0.93312134]]`, det M=1.
- 공기 image-space EFL = 84.14901623305201 mm.
- 마지막 vertex로부터 BFL = 79.92430937348901 mm.
- 무한 물체의 paraxial IMA z = 91.92430937348901 mm.
- 유한 물체에서 paraxial collimation을 주는 OBJ z = −78.52124279501913 mm.

### 2.2 카메라·프로젝터 두 공액

카메라 경로는 infinite OBJ의 field `(0,0)`, `(0.03,0)`, `(0,0.03)` rad에서 시작한다. 먼저 공통 paraxial IMA를 사용하고, 모든 field에 공통으로 적용되는 real-ray RMS best plane과 그 이동량을 별도 출력한다. spot와 FFT가 나빠도 정상 상태와 근사 조건이 정직하면 실패한 렌즈 설계일 수 있으며 엔진 버그라고 단정하지 않는다.

프로젝터 경로는 위 finite OBJ z에 12.8×3.2 mm의 균일 Lambertian source를 놓는다. 320×80, fill=0.64는 연속 근사의 metadata/상수 광량 배율이다. paraxial 수평 FOV 추정은 `2 atan(6.4/EFL)≈8.70°`이며 실제 chief map으로 다시 측정한다. 사용자의 최종 ±12° 목표를 만족하도록 만든 표본이 아니다. 중심·가장자리·모서리의 **출력 방향**, angular spot, far-field W/sr 또는 정규화 광도를 검사한다. 120 lp/mm를 angular 주파수로 이름만 바꾸지 않는다.

역추적/상호성 시험은 같은 개별 ray의 경로와 Snell을 대상으로 한다. 새 camera plane-wave bundle과 finite Lambertian projector의 RI 또는 MTF가 무조건 같다고 검사하지 않는다.

### 2.3 변형 표본

- STOP 반경을 4에서 2로 줄여 pupil clipping과 raw flux를 관찰한다. MTF가 올라가더라도 광량 감소를 숨기지 않는다. MTF가 반드시 올라간다는 보편 명제는 시험하지 않는다.
- 마지막 렌즈의 한 면에 A₂…A₁₂를 모두 넣는다. 반경 a=6에서 각 항의 sag 기여가 ±0.003 mm가 되도록 `A₂m = (−1)^(m+1)·0.003/a^(2m)`로 정하면 계수 overflow 없이 고차항 전 경로를 자극할 수 있다. 유효면 전체의 radicand, 경사, 겹침을 다시 확인한 뒤 채택한다.
- 동일 합성 재질에 명시적 시험용 dn/dT와 CTE를 주고 T₀→T₁, T₀→T₂, T₀ 재실행을 각각 독립 계산한다. 실제 grade의 열거동을 대신하지 않는다.
- 비축 field는 원점·가장자리만 보지 말고 중간값도 사용한다. optimizer가 샘플 사이의 나쁜 거동을 피해서 점수만 개선할 수 있다.

## 3. 독립 단위 표본과 반례

| ID | 표본 / 독립 기대 관계 | 잡으려는 오류 | 한계 |
|---|---|---|---|
| G1 | 단일 평면 경계의 알려진 Snell angle·TIR; n₁=n₂ identity | normal 부호, 분기, 정규화 | 복합 비구면 교점 검증 아님 |
| G2 | 평행판의 각도 복귀와 입사·출사 평행 ray 사이 수직거리 `t(sinθ−cosθ tanθt)` | 축방향 두께와 실제 광로 길이 혼동 | 다중 반사·etalon 제외 |
| G3 | c=0, A₂≠0의 sag와 vertex power c+2A₂ | c=0을 무조건 plane으로 처리 | 큰 aperture에서 paraxial 등식 강제 안 함 |
| G4 | z=r², origin=(0,0,−2), d=(1,0,3)/√10: t=√10와2√10 중 첫 hit | Newton의 뒤쪽 root 선택 | 임의 다중근의 완전성 보장 아님 |
| G5 | 면 간 gap `1−8ρ²+8ρ⁴`: 중심·rim=1, ρ=1/√2에서−1 | center/rim만 보는 렌즈 충돌 검사 | 기계 mount 교차 제외 |
| G6 | STOP이 curved surface에 있을 때 실제 sag 위 target residual | vertex-plane aiming 오인 | 모든 multi-branch 해 탐색 보장 아님 |
| W1 | 직사각 pupil의 separable triangular OTF + 원형 pupil overlap | FFT 축·normalization·주파수 단위 | high-NA vector diffraction 정답 아님 |
| W2 | piston 불변, 알려진 linear phase의 PSF shift와 complex OTF phase | 부호 오류를 MTF 절댓값이 숨김 | phase ramp alias 여부를 별도 확인 |
| W3 | 물리 A→2A에서 raw flux→4flux, normalized PSF/MTF 불변 | 정규화가 진폭 오류를 숨김 | Q·h closure만으로 증명하지 않음 |
| W4 | off-axis 평면파가 flat identical-medium stack을 통과 | launch phase 누락, raw OPL을 OPD로 오인 | focal reference는 별도 구성 |
| W5 | 선언한 실제 pupil mask의 diffraction reference | clipped pupil을 full circle limit로 비교 | stop 변경을 같은 설계라고 취급 안 함 |
| W6 | FFT와 direct scalar RS의 제한된 점/범위 비교 | local expansion·cosine-space scaling | 두 scalar 모델 일치가 Maxwell 검증 아님 |
| W7 | 비정수 배율 N 및 offset sample; phase slope 점검 | N→2N에서 반복되는 commensurate alias | refinement만으로 caustic 해결 안 됨 |
| W8 | afocal hΩ 적분은 dΩ=dl dm/γ; αβ chart Jacobian | 방향코사인을 radians로 오인 | horizon에 가까운 cutoff는 별도 적용성 검사 |
| R1 | 무광학 Lambertian source: I(θ)=L·A·fill·cosθ, Φ=πL·A·fill | 두 cosine·면적·fill 중복 | finite detector bin에서는 bin 평균과 비교 |
| R2 | cone 반각60°: projected Ω=3π/4, Ω=π | irradiance와 intensity 측도 혼동 | 임의 lens pupil은 실제 mapping 필요 |
| R3 | 정상 n=1↔1.5 여섯 면: τ=0.96⁶ | source/ray/pupil transmission 중복 | direct path만; ghost 제외 |
| R4 | 60° 평행 여섯 면의 `(Ts⁶+Tp⁶)/2`와 `((Ts+Tp)/2)⁶`의 불일치 | 매 계면 비편광 reset | full vector diffraction 검증 아님 |
| R5 | g(q)=2q의 zero-blur image E=Q/4, 총전력 불변 | warp Jacobian 누락/중복 | finite pixel integration 포함 필요 |
| R6 | uniform-scene RI와 유한 source edge/corner image 비교 | 서로 다른 두 출력의 혼동 | 대칭 작은 PSF의 직선edge~1/2,corner~1/4는 제한 예 |
| M1 | 비대칭 spot의 chief와 centroid, 왜곡 axis 분모0 | chief mapping 삭제·축상0/0 | 별도 정의를 같은 값으로 강제 안 함 |
| M2 | 서로 다른 field focus와 unequal weights의 공통 optimum | field별 best focus를 한 IMA로 보고 | RMS optimum≠MTF optimum |
| M3 | clipped nonuniform pupil의 알려진 완전 Zernike fit, rank-defect | disk orthogonality 투영 오용 | 계수는 선언한 basis/measure 의존 |
| T1 | 모든 A₂m의 `s^(1−2m)`와 zT(r)=s z₀(r/s) | 고차항 CTE scale·c=0 처리 | thermal drift 전체 증명 아님 |
| T2 | 기준T에서 thermal data 없이 G 계산; 다른λ/T 요청만 적절한 failure | 필요 없는 입력을 전역 필수로 강제 | 누락값을 임의0으로 채워도 안 됨 |
| T3 | n 기준과 λ 기준의 독립 네 조합 | relative index⇒air wavelength 잘못된 결합 | 실제 AGF build와 동일성 별도 |
| T4 | 동일 처방을 T₀로부터 독립 생성; 순서 바꿔 동일 결과 | 누적 geometry·state mutation | mount/stress 영향 미지원 |
| T5 | 모든 길이와λ를 s배: spatial scale만 s배, normalized angular invariant | 단위·Fourier scale 연결 | **실제 heating 시험이 아님** |

## 4. 실행 순서와 실패 판정

1. 재질·설정·입력의 최소 의존성부터 검사한다. 열자료와 RGB 전체가 없어도 G nominal lens-only analysis는 가능해야 한다. 반대로 요청한 계산에 필요한 자료가 없으면 이름 있는 오류가 나야 한다.
2. 분석식 unit specimens로 부호·측도·정규화·단위를 확인한다.
3. 합성 triplet에서 실제 finite/infinite object, chief/marginal, intermediate fields를 연결한다.
4. 한 번에 한 numerical knob를 바꾸어 root/aim residual, pupil support·phase sampling, FFT padding/window, quadrature/bin, thermal finite differences를 각각 평가한다. padding은 pupil phase sampling을 개선하지 않는다.
5. 실패 ray는 `physical clip/TIR`, `solver failure`, `missing data`, `unsupported branch`를 구분한다. solver failure를 단순 zero pupil로 누적하면 PSF/RI가 그럴듯한 수치로 바뀌므로 이 경로를 검사한다.
6. raw power·chief map·complex OTF·focus compensation·상태를 함께 남긴다. zero flux는 normalized PSF/RI를 정의할 수 없으며 정상 zero RI로 취급하지 않는다.
7. tolerance는 analytic roundoff, quadrature error, approximation residual별로 다르게 정한다. 실패한 예를 통과시키려고 의미 없이 완화하지 않는다.

## 5. 작성 에이전트와의 반론·응답 기록

### Materials / thermal

제안: single-G table과 explicit constant n를 구분하고, missing dn/dT가 nominal 계산을 막지 않도록 한다. 각 coefficient와 anchor를 독립 시험한다.

검토 반론: T₀ no-op에서 CTE/dn/dT가 둘 다 없어도 작동해야 한다. 지원하는 경우 geometry-only는 dn/dT, index-only는 CTE가 필요하지 않다. c=0,A≠0는 평면이 아니다. length+λ similarity를 thermal truth로 오인해서는 안 된다. 응답: geometry-only 변환은 CTE만 요청하고 `Material.index`는 CTE를 요청하지 않는다. T=T₀는 thermal data가 없는 no-op을 허용한다. fixed anchors/fixed gaps, front/mid/back anchor와 STOP radius·OBJ·IMA 정책을 분리했다. hot state의 누적 재팽창은 거부하며 기준 처방으로부터 다시 생성한다. 독립 G532 단일 자료 시험에서 nominal trace와 T₀ 변환이 통과했다.

### Metrics

제안: 두 field focus 10/20의 common15, clipped pupil Zernike, cubic distortion, T/S 회전.

검토 반론: 대칭 weight만으로는 weighting bug를 잡기 어려워 unequal weight를 추가한다. chief와 centroid를 구분하는 비대칭 표본과 nominal scale drift를 반드시 남긴다. axis distortion과 rank defect는 명시적 상태로 처리한다. 응답: 비대칭 weighted spot에서 chief/centroid를 모두 남기고, pupil·field의 unequal weight를 적용한 common-focus 시험, clipped/rank-deficient Zernike, nominal 대비1% map drift, arbitrary field T/S 회전 및 real-triplet curvature 시험을 구현했다고 보고했다. 실제 최종 suite 기록으로 실행 여부를 대조한다.

### Wave

제안: circle overlap, ideal sphere, tilt/defocus, clipped pupil, direct scalar RS, angular closure.

검토 반론: rectangle은 축과 cutoff를 더 명확히 검사한다. A²의 **정규화 전** power와 off-axis flat-stack launch phase를 추가한다. high-NA에서 FFT와 exact RS의 차이를 무조건 버그로 취급하지 않는다. canonical plane은 sparse ray maximum으로 정할 수 없다. focus sweep에서 같은 AΠ 유지와 B=1 재설정은 서로 다른 metric이다. 정수 배율 phase alias를 피하려면 noncommensurate 또는 offset sample가 필요하다. 응답: 전체 마지막 clear disk의 보수적 sag/tilt bound로 Π를 정하고, 초점 이동에서 AΠ를 유지한다. trace triangle·grid의 unwrapped phase span을 기록하며 conditional/underresolved 상태를 구분한다. 독립 검토자가 추가한 rectangle/Parseval/정확한 angular chart/alias 시험4건이 통과했다. boundary 재사용은 focal↔afocal 전환을 거부하고 chief·material n·object plane 변경도 fingerprint에 포함하도록 수정되었다.

### Radiometry / imaging

제안과 응답: 무광학 source의 deterministic cosine-weighted QMC에 대해 총 Lambertian flux와 **bin-averaged** polar intensity를 비교한다. normal/oblique plate의 Jones와 scalar product, projected solid-angle aperture와 제한된 slow cos⁴ fixture, raw RI·fill scaling을 검사하겠다고 응답했다. `lens_metric` provenance를 physical image kernel로 받지 않고 physical-source Q/h를 요구한다. 검토자는 A→2A 검사는 wave 모듈에 맡기고, radiometry는 Q→2Q image power 선형성을 담당하는 분업에 동의했다. uniform RI와 finite-edge source는 같은 출력이라고 강제하지 않는다.

### Root / 통합

Root는 위 triplet을 **unoptimized synthetic regression fixture**로 채택했다. 추가 반론으로 center/rim은 통과하지만 중간 반경에서 겹치는 G5 표본과 W7 commensurate alias 표본을 요구했다. 이 두 표본은 단순한 정상 렌즈 성공 사례만으로는 검출되지 않는 실패를 겨냥한다.

## 6. 최종 결과 해석에서 보존할 제한

- production 재료의 신뢰할 수 있는 grade별 n(λ,T), transmission, CTE가 자동으로 검증된 것은 아니다.
- 수렴한 독립 scalar FFT profile은 Zemax cosine-space 내부 구현과의 동일성을 증명하지 않는다.
- clipping·aberration·source weighting은 서로 다르다. raw total Q 하나로 pupil 내 조명분포 차이를 고칠 수 없다.
- 기하광학과 scalar 회절의 적용 한계는 실행 성공과 다르다. 의도적으로 실패시키는 표본이 적절한 상태를 반환하는 것도 성공한 시험이다.
- 개별 테스트·통합 예제의 실제 실행 결과, 남은 오류와 unsupported case는 별도 최종 코드 검토/실행 보고서가 권위 있는 기록이다.

## 7. 독립 검토자 직접 실행 기록

초기 독립 시험 `tests/test_adversarial.py`8건과 `tests/test_wave_adversarial.py`4건을 실제 실행했고 모두 통과했다(0.202 s, 코드 작성 중 snapshot). 이후 nominal G-only FFT/API와 afocal angular-frequency/vergence dispatch 회귀 시험을 추가했다. 전체 통합 수치·최종 실행 상태는 `reviews/integration_camera_review.md`와 `reviews/integration_referee_results.json`에 기록한다. 해석식이 없는 triplet의 MTF 수치를 정답값으로 고정하지 않는다.

최종 독립 suite에는 source_wave actual-triplet global angular power 등록 시험까지 포함되어15건이며, rim/chief zero-weight 표본 추가 후15/15 통과(3.511 s)를 확인했다. Gaussian STOP 적분 표본과 실제 pupil edge 표본을 분리해야 한다는 반론이 source analytic Q 시험에서 발견되어, 현재 기본 설정은 둘을 함께 제공한다.

마지막 exact STOP-rim 경계 반례를 추가한 독립 suite는16건이며16/16 통과(8.814 s). source x=3.2 mm의5개 rim이 aim/aperture tolerance 불일치로 차단되는 실제 결함을 검출했고, 종료 시 실제 STOP predicate를 확인하는 수정 뒤397/397 및 해당 angular MTF 재실행을 확인했다.
