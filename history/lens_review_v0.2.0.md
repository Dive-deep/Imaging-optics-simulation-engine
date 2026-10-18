# v0.2.0 독립 코드 리뷰 — 결상·프로젝터 렌즈 설계 관점

검토일: 2026-10-18. 역할을 부여한 새 AI 리뷰이며 사람 전문가 인증이 아니다. 대상은 `../input/engine/`의 고정 v0.2.0이다. 엔진 소스와 기존 테스트를 수정하지 않았다. 이 보고서의 재현은 `reproduce_lens_review.py`, 실제 값은 `results.json`에 있다. **재현 스크립트의 성공은 발견한 문제가 현재 존재함을 확인한 것이며 제품 합격 판정이 아니다.**

## 판정

추가 **구현 결함 1개와 F-theta의 2D 정의·호환성 문제 1개**를 확인했다. 둘 다 P2로 제안한다. 전체 분석이 근축 광선으로 바뀌었다거나 Snell/기본 T/S 초점 식이 전반적으로 틀렸다는 증거는 찾지 못했다. 오히려 독립적으로 미분한 단일 구면의 비축 T/S 초점은 세 방위각에서 최대 약 `5.02e-8 mm` 오차로 일치했다.

이 결과는 기존 235개 테스트를 부정하지 않는다. 기존 테스트가 다루지 않은 방위각·동일 곡면 접속 조건을 추가한 것이다. 전체 235개 재실행은 부모 리뷰어의 통합 기록과 구별한다.

## LD02 — 동일 곡면 dummy 삽입이 정상 광선을 `escaped`로 바꾼다

**분류:** 확정 구현 결함, P2. **대상:** `src/sequential_optics/geometry.py:267–292`, 특히 `284`, `291`. 허용 계약은 같은 파일 `149–155`에 있다. 해당 docstring은 동일 위치 dummy/refracting 면을 허용한다고 명시한다.

### 처방과 독립 기대값

공기 1.0, 유리 1.5, 앞면 `z=0, R=50 mm`, 뒷면 `z=3, R=-50 mm`, clear radius 10 mm, 앞면 STOP radius 5 mm인 일반 singlet을 사용했다. IMA는 `52.4949494949 mm`다. 첫 면과 **위치·sag가 완전히 같으며 양쪽 매질도 유리인 dummy 면**을 앞면 바로 뒤에 하나 삽입했다.

이 삽입은 추가 거리 0, 굴절률 차이 0이므로 광선 방향·OPL·aperture 통과·spot을 바꿀 수 없다. 물리적 접촉 모델이나 제조 오차의 가정이 필요한 주장이 아니다. 같은 광학 경계를 기록용으로 한 번 더 관찰하는 것이다. 두 처방 모두 `validate_geometry`가 `valid`로 받아들인다.

### 실제 결과

| Field X/Y [rad] | 원래 처방 | 동일 곡면 dummy 추가 | 결과 영향 |
|---|---|---|---|
| 0 / 0 | 217 ok | 217 ok | RMS 동일: 0.04059858628047327 mm |
| 0.01 / 0.02 | 217 ok | **217 escaped** | 원래 RMS 0.04312536723603562 mm인데 추가 후 spot 계산 불가 |
| 0 / 0.1 | 217 ok | **115 ok, 102 escaped** | 통과 quadrature fraction이 1 → **0.6082827156417662**. Spot은 `ok`이며 RMS가 0.09044853455719173 → **0.08936669979161653 mm**로 바뀜 |

모든 광선이 사라지는 사례뿐 아니라, 수치적인 광선 소실을 실제 vignetting으로 처리하여 정상 상태의 잘못된 수치를 내는 사례도 있다. 이 때문에 단순히 “특정 면에서 수치 오류가 발생한다”로만 취급하기 어렵다. 실제 연속 광학계의 throughput 손실은 아니며 위 fraction은 광량이 아니라 STOP quadrature 통과 비율이다.

### 원인

첫 곡면 교점은 허용오차 안에서 계산되지만 sag 방정식 잔차가 정확한 floating-point 0은 아닐 수 있다. 재현의 한 광선에서 첫 hit은 `(4.000000000000016, 1.1061e-14, 0.16025682249210343) mm`, 잔차는 **+5.512812428776215e-13 mm**였다. 이는 기본 교점 허용오차 `1e-10 mm`보다 작다.

같은 곡면을 다시 만나면 그 미소 잔차 때문에 수학적 root가 미소 음의 t 쪽에 있다. 그런데 곡면 경로는 `t>=0` 구간에서 **정확한 0 또는 부호 변경**만 찾으므로 해당 접속을 `escaped`로 반환한다. 평면 경로는 `t>=-tol`을 별도로 허용하고 같은 삽입 대조 실험에서 완전히 동일한 결과를 냈다.

### 수정 방향

직전 면에서 검증한 교점과 동일한 곡면을 연결하는 상황을 명시적으로 인식해 t=0 접속을 보존하거나, 거리 오차·법선 입사·root의 단순성까지 검사하는 원점 접속 규약을 추가한다. **모든 작은 sag 잔차를 root로 인정하는 완화는 권하지 않는다.** 코드가 방지하려는 실제 tangency/다중 root의 조기 오인 문제가 다시 생길 수 있다.

후속 회귀는 축상·비축, 동일 매질 dummy·같은 위치의 명시적 추가 굴절, 평면·곡면, 미소 양의 간격, 실제 tangent 및 더 앞선 다른 root를 함께 구별해야 한다. `escaped`는 실제 물리적 miss일 때만 사용해야 한다.

## LD01 — `f_theta`가 방위각에 따라 달라지는 2D 좌표 정의

**분류:** 확정 정의·Zemax 유사성 문제, P2. **대상:** `src/sequential_optics/metrics.py:298–316`; 원리 문서 `docs/core_principles_ko.md:194`.

이 문제는 문서와 코드를 서로 다르게 구현한 버그가 아니다. **문서도 현재의 componentwise chart를 그대로 정의하고 있다.** 문제는 그것을 회전대칭 광학계에서 통상 사용하는 radial F-theta 왜곡과 같은 이름으로 노출한다는 점이다. X 또는 Y 단일 meridian에서는 두 정의가 같으므로 기존 단축 테스트만으로 발견되지 않는다.

### 방향과 F-theta의 관계

현재 Field는 광선 방향을 다음처럼 만든다.

\[
\mathbf d\ \parallel\ (\tan\theta_x,\tan\theta_y,1).
\]

따라서 실제 광축에 대한 polar angle과 방사 단위벡터는

\[
\theta=\arctan\sqrt{\tan^2\theta_x+\tan^2\theta_y},\qquad
\hat{\mathbf e}=\frac{(\tan\theta_x,\tan\theta_y)}{\sqrt{\tan^2\theta_x+\tan^2\theta_y}}.
\]

축상을 기준으로 한 radial F-theta의 2D chart는 `q = theta * e`다. 현재는 `(theta_x, theta_y)`를 사용한다. 두 축 성분이 모두 있으면 `norm(q)`가 실제 theta보다 크다. 이로 인해 좌표계가 임의로 만드는 왜곡이 렌즈 왜곡에 더해진다.

Ansys의 공식 Field Curvature and Distortion 설명은 F-theta 기준을 물체 공간 chief의 각도와 영상 반경의 관계로 제시한다.[1] 본 프로젝트는 임의 2D field 입력을 지원하므로 이 차이가 직접 노출된다. 공식 도구의 내부 다차원 처리까지 수치적으로 검증했다는 뜻은 아니다.

### 회전 불변성 재현

제공된 `synthetic_triplet()`에서 실제 polar angle을 `0.15 rad`로 고정하고 방위각만 90° → 45°로 바꾸었다. 면·STOP·IMA가 모두 동축 회전대칭이므로 chief 반경과 radial 왜곡은 같아야 한다.

| 항목 | Y meridian | 45° 방위 |
|---|---:|---:|
| Field X [rad] | 약 0 | 0.10646465627572703 |
| Field Y [rad] | 0.15 | 0.10646465627572702 |
| 실제 polar angle [rad] | 0.15 | 0.15 |
| 실제 chief 반경 [mm] | 12.708843304909514 | 12.70884330490951 |
| 현재 `f_theta` [%] | **0.6852198935656716** | **0.3082209286634375** |
| polar F-theta 기대값 [%] | 0.6852198935656716 | 0.6852198935656435 |
| rectilinear 대조 [%] | -0.07105439685840156 | -0.07105439685842625 |

원인은 chief trace가 아니다. 실제 chief 반경과 rectilinear 결과는 회전 불변이며 F-theta reference chart만 잘못된 방사 기준을 만든다.

### 수정 방향

1. `f_theta`를 위 polar-angle chart로 바꾸고 역변환 `tan(theta) * e`를 이용하여 derivative perturbation도 같은 chart에서 정의한다. 축상 극한은 연속적으로 0으로 둔다.
2. 기존 componentwise chart가 필요하다면 별도 이름으로 남기고 radial F-theta와 다르다고 결과 metadata·원리 문서에 명시한다.

비축 `reference_field` 주위의 local Jacobian 정의도 어느 chart의 미분인지 유지해야 한다. 같은 polar angle에서 여러 방위각·큰 field·축상 극한을 비교하는 회귀가 필요하다.

## 독립 대조에서 확인된 기본 계산

### 비축 구면의 T/S 초점

추가 입사각 theta=0.4 rad, R=50 mm, n1=1, n2=1.5의 단일 구면과 vertex STOP을 사용했다. 입사 평면파를 구면 vertex 근처에서 직접 미분하면, `theta_t = asin(n1*sin(theta)/n2)`, `D=n2*cos(theta_t)-n1*cos(theta)`에 대하여

\[
z_T=\frac{n_2R\cos^3\theta_t}{D},\qquad
z_S=\frac{n_2R\cos\theta_t}{D}.
\]

이는 이 리뷰에서 Snell 법칙으로 직접 유도한 **global Z 교점**이며 chief를 따른 전파 거리와 구별했다. 기대값은 `(128.0489504847748, 137.30296565592562) mm`다. 엔진은 방위각 0°, 45°, 90° 모두 `status='ok'`이고 최대 오차 `5.02e-8 mm`였다. 실제 ray-neighbor의 미분으로 T/S focus를 구하는 방향은 합리적이며 Ansys가 공개한 parabasal field-curvature 설명과도 일치한다.[1]

### 추가로 대조·검토한 항목

- **실행 대조:** planar dummy 삽입의 RMS 변화 0, 축상 curved dummy의 RMS 변화 0, rectilinear 왜곡의 방위 불변성. 전부 재현 스크립트에서 assert했다.
- **코드·수식 검토:** A2를 포함한 vertex curvature와 reduced ABCD, EPD/명목 F/# → STOP의 변환, weighted spot와 chief/rim의 분리, LSA 기준면 이동, meridional ray fan과 coma descriptor, 공통 geometric RMS refocus, 모든 A항과 크기의 균일 팽창. 이 항목 전체를 새 독립 수치 검증한 것으로 세지는 않는다.
- EFL/BFL의 근축 정의와 real-ray spot/수차 추적은 서로 다른 역할이다. 현재 분리는 타당하다. 명목 F/#와 finite-conjugate working F/#도 README에서 구별한다.
- 같은 field의 local best focus를 공통 sensor 성능으로 오인하지 않도록 common focus가 한 z를 사용하며, 별도 spot 재평가를 수행하는 구조는 타당하다.

## 개선 권고 — 확인된 결함 수에는 포함하지 않음

1. **Vignetting 수렴 지표:** 기본 polar Gauss quadrature는 smooth circular pupil에 효율적이나 clipping mask의 불연속에는 수렴 속도가 떨어진다. 이미 결과가 quadrature measure와 샘플링 한계를 명시하므로 새 계산식 오류로 세지 않는다. 다만 데이터 생성에는 rings/arms와 별개의 boundary refinement, rectangular 또는 adaptive clipping 적분, 변화량 보고를 추가하는 편이 좋다. Ansys도 강한 비원형 vignetting에서 rectangular sampling을 권한다.[4]
2. **공통 초점의 merit 재평가:** `metrics.py:698`의 큰 2차 모멘트 차감은 큰 defocus와 매우 작은 잔여 spot에서 cancellation 가능성이 있다. 최적 z를 찾은 후 이미 가지고 있는 centered ray coordinates로 RMS를 다시 계산하는 편이 수치적으로 더 안전하다. 이번 리뷰에서는 실용 처방의 명백한 새 실패를 확정하지 않았으므로 버그로 계수하지 않는다.
3. **분석 좌표 계약:** 현재 global-z spot, chief-normal FFT, 고정 X/Y FFT alias, 진짜 meridional T/S 수차는 서로 같지 않다는 문서가 있다. 이들을 소비하는 GUI/데이터셋에도 axis/plane/measure 정보를 필수로 보존해야 한다. 기존에 명시된 scalar·tilt·mount 제한은 새 결함으로 세지 않았다.

## 출처

아래는 제조사 공식 공개 자료이며 2026-10-18 직접 조회했다. 수식 검토의 외부 기준으로 사용했으며 상용툴의 실제 숫자와 대조한 것은 아니다.

1. Ansys OpticStudio, **Field Curvature and Distortion**, 2026 R1 계열 공개 도움말. Parabasal Z-focus, 실제 chief image height와 small-field reference, F-theta의 물체 각도 정의. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Field_Curvature_and_Distortion.html>
2. Ansys OpticStudio, **Standard Spot Diagram**, 2025 R2. Chief/centroid 기준, sampled RMS/GEO, vignetting된 광선의 제외, 기본적인 spot 해석. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Standard_Spot_Diagram.html>
3. Ansys Optics, **How to use Ray Aiming**. Ray aiming으로 물리 STOP을 채우는 과정과 paraxial/real STOP size 차이. <https://optics.ansys.com/hc/en-us/articles/42661778056083-How-to-use-Ray-Aiming>
4. Ansys Optics, **Why is the RMS spot size listed in Spot Diagram different than the values reported by RSCE and RSRE operands?** 표본과 기준점에 따른 RMS 차이, 비원형 vignetted pupil의 GQ/RA 차이. <https://optics.ansys.com/hc/en-us/articles/43071095528979-Why-is-the-RMS-spot-size-listed-in-Spot-Diagram-different-than-the-values-reported-by-RSCE-and-RSRE-operands>

동일 매질·동일 위치 dummy의 광학 불변성과 구면 T/S 검산식은 위 제품 설명을 복제한 것이 아니라 이 리뷰에서 사용한 독립적인 물리 불변량과 Snell 미분 유도다. Zemax/Code V 상관 검증과 실측은 수행하지 않았다.
