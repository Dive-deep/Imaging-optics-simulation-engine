# 독립 2차 코드 리뷰 — 카메라·투사 렌즈 설계 엔지니어

- 대상: 동결된 `sequential_optics_core_v0.1.0` / `input/optical_engine`.
- 관점: 실제 sequential 처방에서 aperture, STOP, real ray, 수차, spot, 공액, 온도 변화의 결과를 설계 판단에 사용할 수 있는가.
- 방법: 구현을 먼저 읽고 독립 해석식과 최소 처방으로 재현했다. 이전 작성 에이전트들의 리뷰 결론을 전제로 삼지 않았다. 엔진 소스는 수정하지 않았다.
- 실행: `python lens_design/reproduce_lens_findings.py`. 다른 위치에서는 `ENGINE_ROOT=/path/to/optical_engine` 지정. 실행 수치: `findings.json`.
- 판정: 주요 기본식은 합리적으로 구현되어 있으나, **현재 버전을 임의 처방의 학습 데이터 생성기에 바로 투입하기에는 아래 재현 결함이 남아 있다.** 특히 평범한 비축 비네팅과 STOP 경계에서도 전체 해석이 실패한다. 정상 fixture의 테스트 통과가 이 조건을 보증하지 않는다.

## 확인된 결함 요약

| ID | 우선순위 | 확인된 문제 | 주된 영향 |
|---|---|---|---|
| LD-1 | P1 | STOP 앞의 물리적 유한 면 이탈을 `aim_unresolved`로 바꾼다 | 정상 비네팅 처방의 spot/MTF/refocus 해석을 중단 |
| LD-2 | P1 | 유한 원통 경계의 거의 0인 판별식에서 부동소수점 음수를 실제 miss로 취급 | 정상 STOP rim 광선 하나가 전체 해석을 실패시킴 |
| LD-3 | P2 | LSA 사용자 기준면을 q=0 실제 초점으로 잘못 사용 | 축상 limit 값과 LSA 곡선이 불연속/오류 |
| LD-4 | P2 | `geo_radius_sampled`가 추적 완료된 zero-weight rim을 제외 | Zemax식 GEO radius/spot 외곽으로 사용하면 과소평가 |

P1은 일반 입력 범위의 분석을 막거나 데이터 수집의 유효/무효 분류를 바꾸는 문제, P2는 특정 결과·옵션의 정확성/의미 문제다. LD-4는 weighted RMS 오류와 구분한다.

## LD-1 — 물리적 비네팅이 수치 미해결로 변환됨

**위치**

- `src/sequential_optics/geometry.py:419–426`, `aim_ray.residual`.
- `geometry.py:175–188`, `_bound_interval`은 `ignore_apertures=True` 상태에서도 physical cylinder를 적용한다.
- `geometry.py:484–487`, `trace_bundle.one`은 예외를 `aim_unresolved` trace로 전달한다.
- `src/sequential_optics/metrics.py:45–55`, `_valid_bundle`이 이 상태 때문에 bundle 전체를 거부한다.

**최소 처방**

공기 중 평면 aperture 두 개: `(z=0 mm, clear radius=2 mm)`, `(z=10 mm, clear radius=2 mm, STOP)`. IMA는 z=20 mm, 물체 무한대, field `(0.1,0)` rad. 중심 시야는 모든 ray가 통과한다. 비축에서는 두 원형 개구의 투영이 `d=10 tan(0.1)`만큼 상대 이동하므로 정상적인 부분 비네팅이다. 렌즈/굴절/caustic이 없어 ray 방향과 차단 여부를 해석적으로 확정할 수 있다.

LD-2와 분리하기 위해 이 사례는 `include_boundary=False, rings=8, arms=24`로 실행했다.

| 항목 | 기대 | 실제 |
|---|---:|---:|
| 동일 샘플의 정상 ray | 134 | 134 |
| 물리적으로 차단/이탈하는 ray | 58 | 58개 모두 `aim_unresolved` |
| 통과 quadrature 비율 | 0.6882578210 | 전체 분석 예외 |
| 동일 샘플의 centroid RMS | 1.1857335813 mm | 전체 분석 예외 |

독립 계산은 STOP 좌표 `p`에서 첫 면 좌표가 `2p−(10tan0.1,0)`라는 관계를 사용했다. 투과 여부는 이 벡터의 길이가 2 mm 이하인지로 결정된다. 연속 면적 기준 두 disk의 overlap 비율은 0.6840062909이며, 위 quadrature 값과의 약 0.00425 차이는 유한 샘플 경계 오차이지 본 결함의 증거로 쓰지 않았다.

첫 면의 **mechanical radius만 4 mm**로 바꾸면 clear radius와 통과 ray는 그대로인데, 58개가 `aperture_blocked`로 바뀌고 RMS가 독립 기대값에 정확히 일치한다. 즉 physically unavailable path와 numerical unknown이 aperture bookkeeping 방식에 따라 뒤바뀐다.

**영향과 수정 방향**

미비네팅 처방만 통과시키는 데이터 편향, 정상 off-axis field의 누락, 최적화 목적함수의 불필요한 중단을 만들 수 있다. Aiming용 수학적 표면 연장과 실제 유한 면 clipping의 분리를 설계해야 한다. 다만 **아무 Newton 시도에서 발생한 escape나 TIR을 일괄 `blocked`로 치환하면 안 된다**. 현재 시도는 잘못된 초기값일 수도 있기 때문이다. 해석적으로 정해진 ray 경로의 physical miss와 실제 미해결을 구분하고, 일반 refracting 처방에는 continuation/검증된 branch 전략이 필요하다. 위 평면 사례는 실제 경로가 유일하므로 미해결로 볼 근거가 없다.

독립 광학 물리 리뷰어도 이 처방이 비축 평면파의 정상적인 disk overlap 문제이며, 미지의 aiming 해가 아니라고 교차 확인했다.

## LD-2 — 실제 광학면과 무관한 bounding-cylinder 접선 수치 오류

**위치**: `src/sequential_optics/geometry.py:180–188`, 특히 `disc<0` 즉시 escape. 후속 `aim_ray.residual:424–426`과 `metrics._valid_bundle:50–52`에서 전체 해석 실패로 확대된다.

**재현**: LD-1 처방에서 첫 면의 mechanical radius=4 mm, 기본 `rings=8, arms=24, include_boundary=True`.

- 217개 trace 중 정상 146, 정상적인 aperture 차단 70, 잘못된 `aim_unresolved` 1.
- 실패한 normalized STOP 표본은 약 `(-1.83697e-16,-1)`, quadrature weight 0.
- 직접 평면 교점 계산: `(-4.44089e-16,-2,10) mm`, 반경은 정확히 2 mm.
- 면 법선·광선 방향 내적: 0.9950041653. **광학면에 대한 grazing이 아니다.**
- bounding cylinder의 계산값: `b²=0.040134137246281794`, `4ac=0.04013413724628181`, 차이 `−1.3877787807814457e−17`.

이 음수는 원통 경계 접선의 반올림 오차다. 실제 aperture와 plane 교점 조건은 만족한다. 현 구현은 평면의 해석적 교점을 구하기 **전에** 원통 discriminant에서 종료한다. 기본으로 추가된 zero-weight rim 하나의 실패가 정상 positive-weight spot의 계산까지 막는다.

**수정 방향**: 평면은 실제 plane hit를 계산한 다음 동일한 radial tolerance로 범위를 판정하는 경로를 우선 사용할 수 있다. 일반 곡면 bounding interval에서는 scale-aware 오차 경계 및 보수적인 interval 계산이 필요하다. 단순히 물리 개구를 임의로 늘리거나 모든 음수 판별식을 0으로 만들면 안 된다. zero-weight topology ray 실패가 각 분석에 어느 정도 필수적인지도 별도 계약으로 정하되, wave pupil을 조용히 잘라서 성공 처리하지 말아야 한다.

LD-1은 물리 loss 분류, LD-2는 수치 bounding 오류다. 같은 예외 메시지로 나타나지만 수정 대상이 다르다.

## LD-3 — 사용자 지정 LSA 기준에서 q=0 물리 초점을 덮어씀

**위치**: `src/sequential_optics/metrics.py:432–444`, `longitudinal_spherical_aberration`.

함수는 사용자 `reference_focus_z`를 받아 `crossing−reference_focus_z`를 보고한다. 그런데 q=0에서는 실제 paraxial image를 계산하지 않고 `crossing=reference_focus_z`로 대입한다. 좌표 원점을 바꾸면 실제 초점 위치까지 바뀌는 오류다.

제공된 synthetic triplet, pupil `[0,0.001,0.1]`, 실제 근축 초점 91.9243093735 mm, 사용자 기준 101.9243093735 mm:

| q | 예상/실제 물리 교점 | 현재 보고 교점 | 현재 LSA |
|---:|---:|---:|---:|
| 0 (근축 limit) | 91.9243093735 | 101.9243093735 | 0 |
| 0.001 | 약 91.9243091870 | 91.9243091870 | −10.0000001865 |

q=0의 LSA는 −10 mm이어야 한다. 작은 q의 실제 ray 결과가 같은 물리 limit를 독립적으로 뒷받침한다. 사용자 기준면은 subtraction 기준으로만 사용하고, q=0 crossing은 별도로 구한 물리 paraxial conjugate로 채워야 한다. 참조를 기본값으로 둘 때에는 드러나지 않는 결함이다.

Ansys의 LONA 정의도 zone=0을 paraxial image 위치로 계산하며, 현재 image에 대한 defocus로 보고한다 [R2]. 엔진이 LONA의 모든 옵션을 복제하지 않더라도 기준면 변경이 물리 초점을 바꿀 수 없다는 점은 동일하다.

## LD-4 — GEO radius에서 이미 추적한 rim을 제외

**위치**: `src/sequential_optics/metrics.py:53`, `_valid_bundle`의 `w>0` mask; `spot_diagram:121–124`; `_moments:107`.

기본 `trace_bundle`은 chief/rim을 zero quadrature weight로 추가한다 (`geometry.py:495–500`). 이 표본은 aperture topology를 위한 것으로, **weighted centroid/RMS에서 0을 주는 것은 옳다**. 그러나 같은 mask로 최대 반경과 표시 좌표까지 제한한다.

제공된 on-axis synthetic triplet, 모두 정상 trace:

| sampling | 보고 GEO radius | 추적된 모든 transmitted ray의 최대 반경 | 과소값 |
|---|---:|---:|---:|
| 기본 8 rings × 24 arms | 0.0087723954 mm | 0.0090322837 mm | 2.877% |
| 보조: 1 ring × 24 arms | 0.0032638407 mm | 0.0090322837 mm | 63.865% |

1 ring 결과는 coarse quadrature의 예시이며 일반적인 오차율로 주장하지 않는다. 기본 설정에서도 차이가 있으며, 이미 추적한 ray를 최대값에서 제외한다는 점이 핵심이다.

Ansys Standard Spot Diagram은 GEO를 기준점에서 가장 먼 transmitted ray까지 거리, 즉 표시되는 ray 전체를 감싸는 반경으로 정의한다 [R1]. 현재 이름 `geo_radius_sampled`를 **양의 quadrature 표본만의 최대값**이라고 매우 좁게 해석하면 내부적으로 일관적이다. 그러나 Zemax와 유사한 GEO/spot 외곽으로 사용자에게 제공하기에는 의미가 다르다. 이것은 RMS 적분식의 오류가 아닌 결과 계약·호환성 결함이다.

수정은 moments의 가중치는 그대로 두고 성공한 모든 표본의 geometric maximum을 따로 구하는 방식이 적절하다. 혹은 기존 값을 `geo_radius_positive_quadrature_only`처럼 명확히 구분하고 all-traced maximum과 표시 좌표를 별도로 반환해야 한다. 최대값조차 연속 aperture 전체 최대값의 보증이 아니라 표본 최대값임은 유지한다.

## 독립 양성 검증과 수식 평가

재현 스크립트에 아래 해석적 control을 포함했다. 앞선 결함들 때문에 모든 계산이 잘못됐다고 판단하지 않는다.

1. 공기 중 R1=50 mm, R2=−50 mm, n=1.5, 중심두께 5 mm인 thick lens에서 독립식 `Φ=Φ1+Φ2−(t/n)Φ1Φ2`, `EFL=1/Φ`, `BFL=(1−tΦ1/n)/Φ`를 사용했다. 코드 EFL 50.8474576271 mm, BFL 49.1525423729 mm가 일치한다.
2. index가 온도 불변인 해당 lens를 CTE=1e−4/K, ΔT=100 K로 자유 팽창시켰을 때 EFL 비율 1.01을 정확히 반환한다.
3. A2–A12를 모두 넣은 asphere에서 `sag_T(sr)=s sag_0(r)`의 차이는 최대 2.78e−17 mm다. `A2m(T)=A2m(0)s^(1−2m)` 구현은 타당하다.
4. 입사각 0.3 rad, n1=1→n2=1.5에서 vector Snell의 transmitted sinθ가 `sin(0.3)/1.5`와 double precision 오차 범위에서 일치한다.
5. code path상 real ray 교점과 vector Snell을 실제 사용한다. 근축식은 aiming 초기값과 EFL/BFL 같은 정의된 일차 특성에 사용한다. 실광선 전체를 ABCD로 대체한 구현은 아니다.
6. `c_vertex=c+2A2`, reduced-coordinate ABCD 및 finite-conjugate 식은 구조상 합리적이다. Field curvature는 chief 주변 real rays의 미분 공액이며, Seidel 값/최소 RMS spot면과 구분한다. Fan의 even/odd 분해를 Seidel coma 계수라고 부르지 않은 것도 적절하다.

## 이미 명시된 제한 — 이번 결함 목록과 구분

다음은 README와 API에서 선언되어 있어 새 수식 오류로 세지 않았다. 다만 사용자 목표와의 적합성은 별도로 판단해야 한다.

- `resolve_aperture`의 F/#는 명목 infinite-conjugate EFL/EPD다. finite working F/#, 모든 Zemax aperture solve를 지원하지 않는다.
- uniform STOP-area quadrature가 렌즈 지표의 측도다. Lambertian source의 power weighting과 동일하지 않다. 물리 source RI/영상은 별도 radiometry/source-wave 경로를 써야 한다.
- `inactive_height_mm`는 실제 단차/측벽 geometry가 아닌 metadata다. mechanical radius까지 analytic sag를 연장한 opaque 모델만 있으므로 실제 기계 무효부 영향은 계산했다고 볼 수 없다.
- thermal 분석은 균일 온도, 자유 등방 팽창, 명시적인 anchor/gap 정책이다. mount/응력 해석을 제외한 사용자 범위에는 합리적이나, anchor 정책을 지정하지 않은 실제 조립품의 유일한 초점 이동을 예측할 수는 없다.
- scalar FFT, chief-normal focal 면과 실제 detector의 비축 투영, afocal angular chart의 정의는 상용툴의 모든 관례를 그대로 복제한 것이 아니다. 이 부분은 다른 독립 물리 리뷰와 함께 읽어야 한다.
- 여기서 확인한 analytic control은 Zemax correlation/실측 검증이 아니다. catalog glass 온도 모델의 정확도를 synthetic fixture로 검증했다고 해석하지 않는다.

## 출처와 보존 범위

[R1] Ansys OpticStudio Help, **Standard Spot Diagram**. GEO는 기준점으로부터 가장 먼 ray, RMS는 전체 ray의 분산 척도이며 vignetted ray를 제외한다는 정의 확인. 직접 문장을 복제하지 않고 요약했다.
https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Standard_Spot_Diagram.html

[R2] Ansys OpticStudio Help, **Aberrations (optimization operands by category)**, LONA. zone=0에 paraxial image, 양의 zone에 real marginal ray image를 사용하고 current image 대비 defocus를 보고한다는 정의 확인.
https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25204/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Aberrations_optimization_operands_by_category.html

[R3] Ansys Optics, **How to use Ray Aiming**. 실제 STOP 위치로 ray를 맞추는 동기와 sequential surface termination, paraxial/real pupil 구분 확인. 이번 LD-1의 실제 기대값은 문서 추정이 아니라 독립적인 평면/disk 해석으로 구했다.
https://optics.ansys.com/hc/en-us/articles/42661778056083-How-to-use-Ray-Aiming

[R4] Ansys OpticStudio Help, **Vignetting Factors**. 비네팅 설계가 정상적인 sequential 해석 대상이고 ray aiming 유무에 관계없이 처리된다는 설명 확인. 이번 엔진에 vendor vignetting-factor 편집 기능 자체가 있다고 주장하지 않는다.
https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v261/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Vignetting_Factors.html

조회: 2026-10-18. 위 출처는 지표 의미 확인용이며, 결함의 핵심 수치 증거는 함께 제공한 실행 스크립트와 JSON이다.
