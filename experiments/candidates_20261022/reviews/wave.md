# 후보 평가의 파동·MTF 방법 검토

작성 범위는 `wave_evaluation.py`와 `validate_wave.py`다. 루트 v0.2.1 엔진 소스는 수정하지 않았다. 이전 환경의 FP16 최적화 엔진은 이 체크아웃에 없어 사용했다고 주장하지 않는다. 이번 평가의 광선·OPL·파동 계산은 FP64다.

## 계산 정의

새 후보는 무한 물체, Y 방향 field, 공통 평면 IMA, RGB 486/588/656 nm를 사용한다. STOP에 Gauss–Legendre rho² 16 rings × 48 arms와 0-weight chief/rim을 놓아 817개 광선을 추적한다. 새 native kernel의 각 굴절면 위치·방향·OPL을 기존 `TraceResult` 형태로 연결한다. 경사 입사면의 초기 위상은 `d_initial dot origin`으로 포함한다. 이 fresh 외부 번들에 기존 엔진의 검증 seal을 임의로 붙이지 않고 `allow_unverified_bundle=True`와 `unverified_external_bundle_fresh_native_trace`를 명시한다. 교점 실패를 물리적 clipping으로 바꾸지 않는다.

기존 `make_boundary_state`가 렌즈 후방의 chief-normal 물리 평면과 STOP triangle 연결을 만든다. 같은 reference sphere의 OPD와 direction-cosine Fourier pupil을 사용한다. 하지만 고NA 조건에서 grid128은 image window가 대략 0.05 mm에 불과하여 OTF native 주파수 간격이 약20 lp/mm가 될 수 있다. FFT padding은 PSF 픽셀 간격만 줄이고 이 window/주파수 간격을 해결하지 않는다. 따라서 0–20 lp/mm의 촘촘한 곡선을 단순 두 점 보간으로 만들지 않았다.

주 출력은 연속 주파수의 scalar pupil autocorrelation이다.

\[
\mathrm{OTF}(f)=\frac{\int P(\boldsymbol\chi-\boldsymbol\delta/2)
P^*(\boldsymbol\chi+\boldsymbol\delta/2)\,d^2\chi}
{\int|P(\boldsymbol\chi)|^2d^2\chi},\qquad
\boldsymbol\delta=\lambda_0\boldsymbol f/n_\mathrm{image}.
\]

이번 image medium은 air, n=1이다. 이 식은 scalar Fourier PSF의 Wiener–Khinchin 관계이며 **주 곡선을 계산할 때 FFT를 실행했다는 뜻은 아니다.** 0–20 lp/mm를 0.5 간격으로 직접 평가한다. 위상을 먼저 wrap한 복소 배열을 보간하지 않고 unwrapped OPD 차이를 구한 뒤 복소 지수를 계산한다. 심한 OPD가 있어도 phase-field FFT의 grid alias를 그대로 받아들이지 않지만, ray mesh가 OPD 자체를 충분히 표현하는지는 별도 문제다.

## Interpolation chart의 변경과 적용 범위

기존 엔진은 물리 boundary의 triangle 꼭짓점을 direction-cosine 공간으로 옮긴 뒤 그곳에서 straight triangle 보간을 수행한다. 고NA·압축된 rim 표본에서는 원래 물리 UV triangle의 orientation이 정상이더라도 이렇게 옮긴 straight chord가 겹치는 사례가 있었다. 대표 pilot0083은 16×48, 32×96, 64×192에서 이 문제가 남았다.

이번 조건부 평가에서는 물리 UV triangle의 reference-sphere OPD를 선형 보간하고, 요청한 angular 좌표를 정확한 역변환

\[
(U,V)=\frac{L(\chi,\eta)}{\sqrt{1-\chi^2-\eta^2}}
\]

으로 찾아간다. 원래의 곡선 support를 유지하며 새로운 ray path나 새 회절 전파식을 만드는 것이 아니다. 다만 **원래 angular-linear 보간과 수치적으로 동일하다고 주장하지 않는다.** 출력에 interpolation chart 변경을 기록하고 mesh refinement를 별도로 요구한다.

적분 측도는 계속 `dchi deta`이며, 기존의 fixed-image lens-metric Fourier amplitude B=1을 유지한다. 이때 UV Jacobian을 amplitude에 곱하면 이중 계수가 되므로 추가하지 않는다. 독립적인 angular-area 진단만 UV triangle 위에서 `L²/(L²+U²+V²)²` Jacobian을 7-point quadrature로 적분한다. 이 면적 진단은 OTF의 실제 정규화 분모가 아니다. 독립 물리 에이전트가 이 좌표·측도 구분을 검토했다.

## Detector plane, RGB와 한계

기존 파동 profile은 chief-normal 국소 면이다. Y field의 평면 sensor로 국소 투영할 때 sagittal X 주파수는 그대로, tangential Y 주파수는 chief_z_cosine으로 나눈다. 이는 국소 scalar 투영이며 high-NA의 정확한 tilted-plane vector diffraction을 주장하지 않는다.

RGB 합성은 각 파장의 complex OTF에 실제 chief image 위치의 translation phase를 포함하여 같은 물리 sensor 원점에서 더한 뒤 magnitude를 취한다. 각 RGB MTF magnitude를 평균하거나 파장별 PSF를 따로 중앙 정렬하지 않는다. 지정되지 않은 광원 spectrum/검출기 response를 추론하지 않으며, 기본 합성은 각 파장의 정규화된 image energy를 동일하게 가중한다.

MTF에는 코팅·Fresnel·흡수 가정을 새로 넣지 않았다. RI의 기하 및 무코팅 Fresnel 분석과 구별한다. NA 0.7–0.8의 scalar uniform-Fourier-pupil 모델은 vector Debye·편광·aplanatic apodization이나 Zemax 수치 동치 인증이 아니다.

## Sampling과 실패 공개

- 기본 angular quadrature 128→256을 비교한다. 모든 비영 MTF점의 순수 상대 변화 ≤0.5%와 정확한 0의 불변을 확인한 경우만 quadrature 조건을 통과한다. max absolute change ≤0.005 진단은 별도로 기록하며, 이를 0.5% 상대 기준의 대용으로 사용하지 않는다.
- 모든 점이 통과해도 ray mesh refinement와 모델 정확도는 별도이므로 `usable_for_final_accuracy_claim=False`를 유지한다.
- 실패는 `underresolved_quadrature` 또는 `not_applicable`로 반환한다. geometric ray-intercept characteristic function은 별도 진단이며 FFT/회절 MTF로 이름을 바꾸지 않는다.
- 별도 actual FFT1024/padding2 진단은 기존 scalar profile을 사용한다. 전체 grid를 triangle마다 검색하던 구현만 같은 bbox/weight/owner 공식의 local grid 범위로 바꿨다. 이때 `underresolved`와 angular triangle overlap 실패를 그대로 남긴다. actual FFT 실패가 별도로 계산한 inverse-UV 주 곡선을 지우지 않는다.
- 사용하지 않은 clear aperture를 실제 광선 footprint+margin으로 줄이면 canonical boundary와 모델 수치가 바뀔 수 있다. 광선과 OPL이 동일한지 aperture-fit 증거를 확인해야 하며, 이를 실제 광학 성능이 좋아진 것으로 표현하지 않는다. 기계적 blank·mount·rim 형상 검증과도 다르다.

## 제어 실험

`results/wave_controls_v2.json`: 11개 제어 모두 통과.

| 제어 | 결과 |
|---|---:|
| 원형 pupil MTF vs analytic diffraction | 최대 절대차 1.08113e−4 |
| pupil tilt의 complex OTF translation 부호 | 최대 복소차 1.12921e−14 |
| 두 RGB 위치의 complex 합성 | magnitude 평균이 숨기는 20lp/mm contrast0 재현 |
| 매끄러운 defocus vs actual FFT | 최대 복소차 4.28494e−4 |
| local bbox interpolation의 mask/owner | 기존과 동일 |
| local bbox interpolation 값 | 기존과 차이0 |
| inverse-UV constant-OPD 원형 pupil | 최대 절대차 1.08113e−4 |
| inverse-UV angular-area Jacobian | 상대차 6.95074e−5 |
| inverse-UV χ-linear tilt, 6×48→12×96→24×192 | 복소오차 .003193→.001066→.000331 |

이 제어는 수학적 관계와 구현 경계의 검증이다. 최종 후보 전체의 수렴·고NA 물리 정확도·제조 가능성을 대신하지 않는다. 초기 `wave_controls.json` 및 `wave_pilot_*.json`은 개발 증거로 보존하며 최종 후보 결과로 바꾸어 표현하지 않는다.

## 후방 boundary의 보수적 상계 개선

후방 plastic asphere에서 A4/A6 등의 부호가 상쇄될 때 기존 profile의 whole-disk positive-term 상계가 실제 sag보다 매우 클 수 있다. 새 평가 전용 private profile은 회전대칭 clear surface의 정확한 방위각 최대

`f(r) = chief_z * (vertex_z + sag(r)) + hypot(chief_x, chief_y) * r`

를 반경 32,768개 구간에서 보수적으로 감싼다. Conic sag의 단조성은 curvature 부호에 따라 endpoint를 고른다. 각 signed polynomial 항의 구간 endpoint 거듭제곱과 radial 항의 upper bound를 합한다. 구간 전체를 포함하는 envelope이므로 단순한 최대값 샘플링과 다르다. 전체 clear aperture에서 conic domain이 양수임을 먼저 확인한다. 기본 산술마다 `nextafter`로 바깥쪽 반올림하고 `1e−12 mm + 64 eps × magnitude`의 추가 여유를 둔다. 실제 원본 엔진 파일을 바꾸지 않고 primary autocorrelation 및 optional FFT의 private profile에만 같은 정책을 적용한다.

`results/wave_projection_bound_controls.json`: random asphere/tilt 80개, 상쇄 다항식 1개, 음의 curvature에서 undefined rim 거부 1개, 합 **82개 통과**. 각 random case에서 32,769개 dense 실제값이 계산 upper bound를 넘지 않았다. 수학적으로 signed-term interval envelope은 보수적이며, dense 시험은 구현 확인이다. 이를 모든 플랫폼 transcendental library의 형식적 IEEE 증명으로 표현하지 않는다.

## 최종 단계 대표 2처방의 ray-mesh refinement

`results/wave_mesh_probe.json`은 `refine_pilot/seed_0012.json`과 `qualified_round1/seed_0108.json`의 저장된 처방을 그대로 사용했다. 최종 선정 파일과 동일하다는 전제는 두지 않는다. 588 nm, 0F와0.8F에서 16×48(817ray)→32×96(3,169ray)을 비교했고 각 실행의 angular quadrature는128→256이었다. 실제 FFT1024/padding2는16×48 실행에서만 별도 진단했다. 한 worker·수치 thread1개로42.950초가 걸렸고 입력 SHA256 및 evaluator/model/native/script SHA256을 결과 파일에 기록했다. 측정 전후 source 변경0건이다.

| 대표 처방 / field | T MTF20, 16×48→32×96 | S MTF20, 16×48→32×96 | 전체 곡선 최대 절대 변화 | 전체 곡선 최대 상대 변화 |
|---|---:|---:|---:|---:|
| seed0012 /0F | .799912→.797418 | .799912→.797418 | .002495 | .312858% |
| seed0012 /0.8F | .392519→.404751 | .060351→.068732 | .017443 | 163.144% |
| seed0108 /0F | .219884→.218548 | .219884→.218548 | .004662 | 1.79585% |
| seed0108 /0.8F | .062222→.069730 | .029449→.012843 | .022317 | 129.297% |

상대 변화의 분모는 finer32×96 값이며, 작은 비영 값을 제외하거나 절대 tolerance를 더하지 않았다. seed0012/0F만 전체 T/S 곡선의 ray-mesh0.5% 조건을 통과했다. 이 경우 angular quadrature의 최대 상대 변화도 두 ray mesh에서 각각 .149778%/.247158%로 통과했다. 따라서 **이 저장된 처방의588nm 축상 scalar 모델에서는 MTF20≈0.80, 정규화 AUC≈0.926이라는 수치 요약이 시험한 두 분해능에서 안정적**이라고 말할 수 있다. 이것은 전체 RGB·field·최종 후보의 정확도 또는 실제 광학계의 보증 하한이 아니다. 실제 FFT의 상태는 여전히 `underresolved`였다.

나머지 3개 조합은 점별 기준을 통과하지 못했다. 특히 seed0012/0.8F의 S-AUC는 .18675446→.18675544로 거의 같았지만 곡선의 최대 상대 변화는163%였다. **AUC 안정성이나 평균 오차를 근거로 개별 MTF값을 신뢰 가능한 DB 학습 label로 바꿀 수 없다.** seed0108/0.8F의 S-AUC도 .364186→.354447로 변했다. 이들은 표시·초기 탐색용 수치 진단이며, 더 조밀한 ray mesh/적분 및 모델 검토가 필요한 상태다.

Actual FFT 진단은 seed0012 두 field에서 `underresolved`, seed0108 두 field에서는 straight angular-triangle overlap으로 `not_applicable`였다. 이 실패를 숨기지 않았으며 geometric MTF 또는 inverse-UV autocorrelation으로 이름만 바꾸어 FFT 통과로 표시하지 않는다. `underresolved` 곡선에서 실제 성능의 비자명한 엄밀 하한을 얻었다고 주장할 근거는 없다.

## 최종 C01의 후방 공간 ray-map fold 진단

최종 wave 평가 중 C01의 0.8F/1F에서 `Sampled ray map has a fold/caustic or singular Jacobian`이 발생했다. 별도 읽기 전용 `wave_recovery.py`로 `results/frozen_selection.json`의 C01 처방을 조사했다. 588 nm, 16×48→32×96→64×192를 사용했으며, 증거는 `results/wave_recovery_C01_diagnosis.json`에 있다. 입력 선정 파일 SHA256 및 평가 소스의 실행 전후 SHA256을 기록했다. 최종 `wave_evaluation.py`, 처방 및 최종 wave 결과는 바꾸지 않았다.

| field | canonical boundary 뒤 reference 거리 | canonical UV map 음의 Jacobian triangle 수, 16→32→64 rings | 64×192 최소 Jacobian |
|---|---:|---:|---:|
| 0.8F | 2.05619 mm | 8 → 24 → 80 | −12.6915 |
| 1.0F | 2.10721 mm | 44 → 168 → 662 | −41.6685 |

모든 추적 광선이 유효했고 최소 ray/boundary cosine은 약 .688이다. 전체 마지막 렌즈면의 보수적 boundary 대신 **추적한 실제 exit footprint 직후**에 평면을 두는 진단에서도 64×192에서 음의 triangle은 각각 30개와 634개였다. 더 후방으로 평면을 이동해도 공간 map의 regularity는 복원되지 않았다. 실제 S1–S5 교점의 XY map은 양의 Jacobian이지만, S6 자체의 XY map은 0.8F에서 30개, 1.0F에서 622개의 음의 triangle을 보였다. 따라서 이것은 앞서 보정한 angular straight-chord 보간 artifact와 다르며, 사용하지 않는 clear aperture가 정하는 plane만 옮겨 해결할 수 있는 현상도 아니다. 더 조밀한 mesh에서도 실제 후방 공간 ray map의 fold가 남는다.

현재 profile은 단일값 공간 boundary 위상·진폭을 전제로 하므로 해당 실패를 유지하는 것이 맞다. 일부 광선을 삭제하거나 fold guard를 우회하여 scalar MTF를 생성하지 않는다. 각 field의 실제 image-ray 분포로 계산한 **geometric MTF 진단 곡선**은 별도 panel에 제공할 수 있지만, 회절/FFT MTF의 대체 검증값이나 .5% 정확도 보증값은 아니다.

같은 64×192 광선의 output-direction XY map은 두 field에서 양의 Jacobian이었다. 이는 향후 momentum/angular-spectrum 좌표, multi-branch wave propagation 또는 Huygens 적분을 검토할 단서다. 다만 좌표의 regularity만으로 phase generating function, 수송 진폭·apodization, caustic/Maslov 위상과 고NA 편광 처리가 검증된 것은 아니다. 이번 결과에 새 전파 모델을 급히 혼합하지 않았으며 별도 물리 개발·수렴 검증 과제로 남긴다. 이 진단은 C01의 지정된 두 field/파장에 대한 것이며, 다른 후보의 실패 원인이 모두 같다고 일반화하지 않는다.
