# 초기 후보의 광학 물리 및 상대조도 검토

이 문서는 새로운 후보 생성 실험의 독립 물리 검토다. 원본 `src/`와 과거 검증 기록을 바꾸지 않았다. 사람 전문가 인증, Zemax 수치 동치 또는 완성 렌즈의 성능 보증은 아니다.

## 조건의 의미와 물리적 난점

상고6.4mm와 반화각12°는 rectilinear 기준 EFL30.1096327mm에 해당한다. EPD50mm이면 약f/0.602다. 일반적인 근축 NA≈1/(2F/#)를 실제 고NA 값으로 사용하면 안 된다. `radiometry.py`는 실제 image-space 방향을 사용하며, 588nm 축상에서 정상적으로 전달된 명목 pupil rim48개 전체의 최소 `n_image*sin(u)`를 기준값으로 사용한다. 일부 내부 광선 한 개의 큰 각도를 NA 합격 근거로 사용하지 않는다. 양·음 T/S rim 사이 cone NA와 RGB의 값도 함께 기록한다.

전체 원형±12° 시야와 EPD50mm를 모두 손실 없이 지름12.8mm image circle에 전달한다면, étendue 용량으로 필요한 image-side projected acceptance의 등가 NA는 `25*sin(12°)/6.4=0.812155` 이상이다. 이것은 시야와 pupil 전반의 광량 가정에 따른 조건이며, 모든 상황의 축상 NA가 반드시0.812155라는 법칙은 아니다. 사용자의 대안인 첫 렌즈 opening50mm는 EPD50mm와 구별한다. object-side STOP 지름과 첫 surface clear 지름을 따로 저장하고, 축상 pupil에서 일부 광선이 차단되면 명목 STOP 지름을 검증된 EPD로 표시하지 않는다.

양의 plastic 두 장과 음의 고굴절 glass만으로 색수차를 보정하면 positive Petzval curvature가 남기 쉽다. 양의 PMMA / 음의 PC / 양의 구면 glass 같은 배열을 포함하고 두께·간격·비구면을 함께 최적화하는 것이 합리적이다. Thin/contact 근사에서 색수차·Petzval·총 파워를 동시에 맞춘 해는 구조의 난도를 진단하는 도구이며, 그 자체를 실제 제작 가능한 처방으로 쓰지 않았다.

## 기하 성능과 합격 조건

`evaluate_design(design, trace_result=None, rings=16, arms=48)`가 JSON 직렬화 가능한 결과를 반환한다. 기본값은 RGB486/588/656nm, field0/.5/.8/1이며 모든 계산은 하나의 고정 영상면을 사용한다.

- 각 파장·field에서 원래 entrance STOP 면적 가중치를 사용해 centroid, RMS radius와 diameter, covariance를 계산한다. Spot 통계 자체는 통과 광선에 대한 조건부 통계이며, vignetting과 실패 수를 같이 기록한다.
- Polychromatic spot은 RGB의 개별 중심을 이동시키지 않고 같은 sensor 좌표에서 계산한다. RGB incident radiance의 가중치는 동일1/3이다.
- Chief height와 rectilinear distortion을 보고한다. 왜곡 기준은 해당 파장의 ABCD EFL이며, 실제 고정 영상면에서 미소각 chief slope를 사용한 별도 왜곡도 보존한다.
- 마지막 surface vertex부터 고정 영상면까지의 거리와 근축 BFL은 서로 다르다. 두 값을 각각 보고하며 보수적으로 둘 다3mm 이상인지 별도로 검사한다. 후자는 BFL 정의의 모호함을 드러내기 위한 추가 검사다.
- CT는 엄격히2mm 초과, 축상 air gap은1mm 이상이어야 한다. 각 clear aperture에서1025개 반경 표본으로 인접 surface의 두께/간격이 양수인지 검사한다. 이는 조밀한 표본 검사이며 연속구간의 수학적 보증은 아니다. Conic radicand는 전 구경에서 양수여야 한다.
- Glass의 conic과A4…A12는 모두0이어야 한다. Polynomial 슬롯은A4/A6/A8/A10/A12 다섯 개만 허용한다.
- 상고 target은6.4mm다. 기본 수치 근접 판정0.032mm는 선언된 계산상 근접 기준이며 사용자가 허용한 설계 공차가 아니다. 최종 호출에서 더 작은 값을 지정할 수 있고 실제 tolerance와 residual을 결과에 기록한다.
- `all_declared_numeric_constraints_and_sampled_rays_passed`는 모든 기준과 모든 chief/rim/내부 표본의 전송을 함께 요구한다. 수치 실패·TIR·clear aperture 이탈·backward projection을 삭제하거나 성공 광선으로 재분류하지 않는다. 이 상태가 최적 MTF 또는 실제 제작 가능성을 뜻하지는 않는다.

## 상대조도 정의

무한 물체가 균일 radiance를 가지고, 모든 powered surface 앞의 고정 STOP 평면에서 면적을 적분한다고 가정한다. 회전대칭 chief map을 `r(theta)`라 하면, object solid-angle annulus와 image annulus의 면적 변환으로

`E(theta)/L = A_eff(theta)*cos(theta)*sin(theta)/(r(theta)*dr/dtheta)`

를 사용한다. 축상 극한은 `A_eff(0)/f0^2`이며 `f0=dr/dtheta|0`다. 유효 면적은 원래 STOP 면적에 유효 광선의 원래 quadrature weight와 선택적인 Fresnel 전송을 곱한 값이다. **손실 광선을 제거한 뒤 RI 면적을1로 다시 정규화하지 않는다.**

Rectilinear map `r=f*tan(theta)`와 동일 pupil이면 cos⁴ 법칙을 회복한다.12°에서0.915414다. Sine map `r=f*sin(theta)`는 같은 가정에서 균일 기하 RI를 준다. 실제 처방에서는 수치 chief derivative를 h=2e−5rad와h/2로 계산하고 차이를 기록한다. 각 RGB의 irradiance를 먼저 동일 incident-radiance 가중치로 합하고 그 합의 축상값으로 정규화한다. RGB RI 값을 단순 평균하지 않는다.

RI는 다음 두 가지다.

1. `RI_geometric`: aperture/TIR/trace 전송만 반영한 lossless 기하값.
2. `RI_uncoated_fresnel`: 실제 surface normal과 입출사 방향에서 구한 각 경계의 평균 s/p Fresnel power transmission을 곱한 값. 무코팅·흡수 없음·각 경계의 scalar unpolarized power 근사다. Polarization/Jones vector는 운반하지 않는다.

이 RI는 smooth chief mapping에 기반한 국소 기하 방사측정 근사다. 큰 수차의 PSF를 sensor 전체에 쌓아 계산한 검출기 조도나 OpticStudio의 exit-pupil projected-solid-angle 적분과 동일한 계산을 수행했다고 주장하지 않는다. Mapping이 정의되지 않으면 RI를 null과 상태로 반환하며0으로 숨기지 않는다.

공식 정의 참고: [Ansys Relative Illumination](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Relative_Illumination.html).

## 고NA scalar MTF와 좌표 보간 검토

NA0.7–0.8에서 scalar FFT와 chief-normal 면의 국소 sensor 투영은 조건부 모델이다. Vector diffraction·polarization·실제 측정/상용툴 상관을 인증하지 않는다. [Ansys FFT와 Huygens 설명](https://optics.ansys.com/hc/en-us/articles/42661795535251-What-is-the-difference-between-the-FFT-and-Huygens-PSF)은 FFT 적용에서 느린F/#, 작은 chief angle과 pupil aberration 조건을 설명한다. 본 목표는 그 전형적인 적용 범위보다 빠르다.

128 pupil에서 고NA 방향 코사인 간격은 작은 image window와 거친 OTF 주파수 간격을 만든다. FFT padding만 늘리면 PSF pixel 간격은 작아져도 OTF 주파수 간격이 개선되지 않는다. 따라서0–20lp/mm의 유효 곡선은 직접 continuous-frequency pupil autocorrelation 또는 충분한 pupil 격자를 요구한다. 직접 autocorrelation을 실행한 곡선을 실제 FFT 실행 결과라고 부르지 않는다. 기하 ray-intercept OTF는 회절을 포함하지 않는 별도 진단이다.

Boundary UV와 direction cosine의 사상 `chi=UV/sqrt(L²+|UV|²)`는L>0에서 일대일이다. 정점만 angular 좌표로 옮겨 직선 triangle을 연결하면 강한 비선형 rim 영역에서 실제 사상에는 없는 chord overlap이 생길 수 있다. UV 삼각망을 보존하고 각 angular query를 `UV=L*chi/sqrt(1-|chi|²)`로 역변환해 support와 unwrapped OPD를 조회하는 방법은 합리적인 조건부 수치 대안이다.

이때 적분 측도는 여전히 `dchi_x*dchi_y`이고 기존 선언된 Fourier pupil 진폭이B(chi)=1이므로 진폭에 UV Jacobian을 추가하지 않는다. UV에서 면적을 독립 적분할 때만 `L²/(L²+u²+v²)²`를 사용한다. OPD 보간 chart는 원래 linear-chi에서 linear-UV로 바뀌므로 원본과의 수치 동일성은 주장하지 않는다. Constant OPD, linear-chi tilt, analytic circular support, traced mesh와 quadrature의 별도 수렴 제어가 필요하다. 이 좌표 변경이 물리적 Debye/aplanatic apodization을 추가한 것은 아니다.

Unused clear aperture를 줄여도 모든 RGB/field의 전체 ray footprint와 margin을 유지하면 실제 전송 광선은 같을 수 있다. 그러나 이에 따라 canonical boundary가 바뀌고 scalar MTF/triangulation이 바뀐다면 boundary/model 민감도다. 물리 성능이 개선됐다고 해석하지 않는다.

## 완료한 독립 검증

`validate_radiometry.py`의15개 제어가 모두 통과했다. 원자료는 `results/radiometry_controls.json`에 있다.

- 6개 평면 경계의 폐형식 Fresnel 전송과 최대 절대차1.11e−16.
- nd와Vd의 네 재료 anchor 확인. Ri를550nm로 오해하지 않음.
- 손실 ray를 포함한 비균일 가중 centroid/RMS와 전송 면적 보존.
- 실제 고NA pilot0083에서 RGB×4field×6chief/marginal/중간 ray의72개를 원본 Python 교점·gradient·Snell 코드와 독립 비교. Status 불일치0, 최대 위치차2.14e−11mm, 방향차2.15e−13, OPL차3.54e−11mm, 상면좌표차9.34e−12mm.
- 실제 original infinitesimal ray와 ABCD EFL 차이1.47e−12mm, BFL차1.06e−13mm.
- Pilot0083의 축상588nm rim48개 NA는0.7092790163905844–0.7092790163905850. 한 광선의 극값을 사용하지 않음.
- 별도 더 작은 chief derivative step으로 RI를 확인했으며 green edge RI 차이는5.19e−10.
- 강제로 aperture를 줄인 반례에서 nominal EPD와NA 합격을 거부하고 모든 실패 수를 유지.

Pilot 수치를 전체 최종 후보의 성적으로 일반화하지 않는다. 최종 후보마다 반환된 constraint 결과와 sampling/MTF 상태를 읽어야 한다.

추가 afocal 평면 stack 대조에서도 무한 EFL/BFL을 양의 유한 결상 해로 합격 처리하지 않음을 확인했다.

## 최종 동결된24개 후보의 독립 감사

`results/frozen_selection.json`의 SHA256은 `ec008e48d0b5156c793f9d04fb8224e8e91c28c3029479dcef9017c9caab3435`다. `validate_selection_physics.py`가 저장된 합격 boolean만 읽지 않고24개를 새로 추적해 검사했다. 결과는 `results/final_selection_physics_audit.json`, 범위 집계는 `results/final_selection_physics_statistics.json`이다.

**24개 모두 검사 통과**했다. 서로 다른 parent seed24개와 서로 다른 처방 hash24개이며, 13field×RGB×817ray×24 = **764,712개 새 광선에서 실패0**이다. 전체 rim의 RGB 최소NA≥0.7,588nm 상고 오차≤0.001mm, 재료·유리 standard·CT·axis gap·두 BFL 정의를 확인했다. Surface 간격은 독립적으로4097개 반경 표본에서 다시 계산했다.

공개 처방의 radius를 curvature로 되돌리고 A4…A12·conic·clear radius·material-after·vertex z·thickness chain·STOP·image plane을 재구성해 native 입력과 일치하는지 검사했다. 기존 결과와 새 계산의 spot·RI·chief·왜곡을 비교하고, 저장 MTF proxy 곡선에서 AUC와 선택 score를 독립 합산해 순위와 일치함을 확인했다. 원래 aperture-fit의 전체 ray 수·처방 hash·전송 이력 보존도 검사했다.

| 최종24개 범위 | 결과 |
|---|---:|
| 선택된 구경 해석 | **모두 첫 렌즈 전면 opening50mm 대안** |
| 실제 지원되는 축상 EPD | 43.200000–43.959376mm |
| 588nm 축상 complete-rim 최소NA | 0.707111–0.723077 |
| 각 처방 RGB 중 가장 작은 축상NA | 0.701772–0.719056 |
| 실제 상면 vertex BFL | 3.000670–9.027529mm |
| 588nm 근축 BFL | 3.014000–9.617752mm |
| RGB 합성 edge RI, 무코팅 근사 | 0.856233–0.988907 |
| RGB 합성 edge RI, 기하 근사 | 0.860721–0.994833 |
| RGB 합성 축상 spot RMS radius | 70.98–169.95µm |
| RGB 합성 full-field spot RMS radius | 114.63–428.09µm |
| 588nm full-field distortion | −2.092566–+1.537250% |
| 588nm target6.4mm 최대 오차 | 0.000506802mm |

선택된 후보에 **EPD50mm 조건을 만족한 처방은 없다.** 사용자가 허용한 첫 렌즈 opening50mm 대안을 만족하는 후보들이다. 두 조건을 합쳐 EPD50mm라고 보고해서는 안 된다.

상고6.4mm는588nm 기준이다.486nm의 full-field chief는6.292902–6.359159mm,656nm는6.417300–6.443117mm다. 모든 색이6.4mm에 맺히는 조건을 달성한 것이 아니며, red chief 일부가 이상적인6.4mm sensor edge 밖에 있다. 이 결과의 RI는 detector crop와 큰 spot의 active-area 손실을 계산한 것이 아니다. 후속 실제 sensor 기준 DB에서는 이 색별 위치와 유한 detector 적분을 별도로 고려해야 한다.

NA≥0.7 판정은 **축상**이다. 최종 후보의 full-field tangential opposite-rim cone NA는 RGB 전체에서0.659662까지 내려간다. 전 field에서NA≥0.7이라고 주장하지 않는다.

선택 score는 `0.6*mean normalized geometric-MTF AUC + 0.4*minimum normalized geometric-MTF AUC −0.25*|RI_edge−1|`라는 명시적인 초기 탐색 heuristic이다. Scalar 파동 MTF의 정확도 순위나 전역 최적해 순위가 아니다.24개 중23개에 centroid 기준 geometric OTF의 음의 실수부 표본이 존재한다. 이것만으로 복소 OTF의 절댓값이0을 통과했다고 단정할 수는 없지만, magnitude-AUC가 단조롭게 높은 MTF 곡선을 보장하지 않는다는 점을 드러낸다. 각 후보의 실제 곡선·phase 관련 진단을 함께 읽어야 한다.

## 최종 파동 좌표·boundary 상계의 추가 코드 검토

`InverseMappedBoundaryPupil`의 보수적 angular bounding rectangle, exact inverse query, uniform-angular 측도, complex OTF 부호와 T방향 주파수의 `1/cos(chief)` sensor 투영을 확인했다. RGB 공통 sensor 원점의 phase를 유지하는 구현과 수식 부호가 일치하도록 방법 문서의 수식을 정정하도록 요청했고 반영되었다.

`conservative_surface_projection_upper`는 반경r에서의 정확한 azimuth 최대값 `axis_z*(z+sag(r))+hypot(axis_x,axis_y)*r`를 사용한다. 양의 axis_z 조건 아래 conic의 곡률 부호에 따른 단조 끝점, 각 signed polynomial 항의 구간 끝점 범위, 전체 clear aperture의 양의 conic domain 검사를 확인했다.32768개 구간의 합집합 상계이며 단순히 조밀한 점의 최대값을 취한 것이 아니다. 기본 연산의 outward rounding과 명시적 공학적 여유를 쓰며, 모든 초월함수 구현에 대한 형식적 반올림 증명이라고 주장하지 않는다.82개 analytic/무작위/경계 제어는 `results/wave_projection_bound_controls.json`에 있다.

이 검토에서 선언된 조건부 모델 안의 차단 결함을 확인하지 않았다. 그래도 고NA 벡터 회절의 물리 정확도, 실제 pupil apodization, ray mesh와 적분 격자 수렴은 별도 문제다. 최종 파동 결과의 `underresolved`·실패 상태를 합격으로 바꾸지 않는다.
