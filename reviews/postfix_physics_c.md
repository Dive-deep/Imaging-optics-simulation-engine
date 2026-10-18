# C1–C4 수정 후 독립 광학 물리 리뷰

검토일: 2026-10-18. 검토자: 구현에 참여하지 않은 별도 광학 물리학자 역할의 AI 에이전트. 작업 소스 `optics_core_v0_3/src/`와 변경 전 고정 사본 `optics_v0_2_independent_review/input/engine/src/`의 차이를 직접 읽었다. `AGENTS.md`, `HANDOFF.md`, `README.md`, `FIX9_REPORT.md` 및 관련 계산 명세도 확인했다. 엔진 소스는 수정하지 않았다.

## 판정

검토한 변경에서 새 차단 결함을 확인하지 않았다. **C2는 잘못된 광학 geometry를 정상으로 간주해서 생긴 문제가 아니라, 정상적인 작은 cell 면적을 큰 전역 좌표로 계산하면서 생긴 수치적 상쇄였다.** 이번 수정은 상대 좌표와 전력 보존 관계를 유지한다. C1은 이미 인증된 동일 경계의 순차 재방문을 거리 0으로 처리하며 실제 광선 경로·OPL·손실을 바꾸지 않는다.

## C2: geometry 및 전력 적분

관련 코드: `source_wave.py:277–367`.

1. 원래 반례에는 렌즈 면 자체가 필요 없다. 정상적인 직사각형 PSF cell, 증가하는 detector edges, 공통 평행 이동만으로 충분하다. 따라서 렌즈 면 교차·꼬임·부호 뒤집힘을 실패 원인으로 설명할 수 없다.
2. Chief-normal 평면에서 detector로 옮기는 선언된 선형 변환을 \(J\)라 하면, FFT cell의 detector상 면적은 \(A_d=\Delta x\Delta y/|\det J|\)다. 코드의 cell 면적 및 꼭짓점 변환이 이에 일치한다. 직교정규 basis에서 만들어진 가역 변환은 원래의 볼록 사각형을 여전히 볼록하게 유지한다.
3. Cell 질량 \(m_k\)가 detector pixel \(P_j\)에 기여하는 값은 \(m_k |C_k\cap P_j|/|C_k|\)다. 전체 격자가 cell을 포함하면 합은 \(m_k\), 잘라내면 실제 overlap의 합만 남아야 한다. 수정 코드의 최종 합에는 강제 재정규화가 없다.
4. 큰 원점을 먼저 빼고 각 cell 중심에서 clipping한 뒤 triangle cross products로 면적을 구하는 변경은 그 적분식을 바꾸지 않는다. Clipping 경계 좌표를 정확히 지정하는 것도 다음 half-plane 판정에서 불필요한 반올림을 줄이는 방법이다.

### 구현과 다른 독립 적분 대조

`run_postfix_physics_c.py`에서 polygon clipping이나 shoelace/triangle area를 사용하지 않는 기준을 작성했다. 변환된 cell을 두 개의 선형 strip 부등식으로 표현하고, 각 x에서 허용되는 y 구간의 길이를 detector y 구간과 교차시킨 뒤 적분했다. Piecewise-linear 경계들의 교차점에서 구간을 분할했다.

- 독립 tilt/azimuth 12개, 각각 정상 정규직교 basis와 양의 가역 \(\det J\) 확인.
- 비등방 cell, `(5,-7)` 또는 `(1000,-1200)` mm 전역 위치, 비대칭 부분 crop 사용.
- 각 경우 12×10 detector pixels를 직접 비교.
- 최대 pixel 확률 차이 **1.0408340855860843×10⁻¹⁷**.
- 실제 포착률은 약 **0.36254–0.55767**이고 독립 적분값과 일치. 1로 재정규화되지 않았다.
- 이 검사는 선언된 affine 면적 변환의 정확성을 확인한다. 높은 chief tilt에서 이 affine 모델이 실제 회절 전파를 정확하게 나타낸다는 검증은 아니다.

## C1: 동일 경계의 0 거리 처리

관련 코드: `geometry.py:304–339`, `geometry.py:386–407` 및 aiming의 동일 경계 처리.

- 이전 면을 성공적으로 통과한 직후, 처방상 z·곡률·conic·A2…A12가 정확히 같은 다음 면에서만 교점을 재사용한다. 임의 외부 ray 또는 단지 가까운 면을 자동으로 경계에 붙이는 처리는 아니다.
- 교점은 좌표를 다시 보정하지 않고 그대로 사용하므로 추가 OPL은 \(n\times0=0\)이다. 법선은 현재 면 gradient에서 계산한다.
- 새 mechanical domain, active aperture, STOP 및 매질 변화의 Snell은 여전히 적용한다. 같지 않은 면의 first-root 절차는 그대로다.
- 무굴절 dummy에서 Fresnel 전력은 1이며, 길이 0인 Beer factor는 1이어야 한다.

독립 추가 대조에서는 흡수계수 0.027 mm⁻¹인 유리 singlet의 앞·뒷면을 각각 같은 매질의 dummy로 중복했다. 축상/비축 2개 field, pupil 3개 지점에 대해 수정 후 두 처방을 동일한 실제 입사 ray로 추적했다.

- 6개 ray 모두 성공. 대응 교점과 OPL은 부동소수점 값까지 동일했다.
- `unpolarized_jones_power`와 `scalar_interface_average` 두 손실 모델의 transmission/absorption/reflection을 대조한 12개 조건에서 최대 차이는 **5.551115123125783×10⁻¹⁶**.
- 전력 ledger residual은 모두 **2×10⁻¹⁵ 미만**.

동일 위치에 서로 다른 매질 경계를 연속 선언하는 경우의 기하 Snell 합성과 실제 coherent 초박막의 전력은 다른 문제다. 현재 엔진은 코팅 간섭/다중반사를 계산한다고 선언하지 않으므로, 그런 물리 현상까지 이 수정의 검증 범위로 확대하지 않는다.

## C3/C4 및 MTF 출력

`registered_spectral_sum` 변경에서 전력 배열의 shape·비음수·유한성·공통 물리 격자를 확인하고, irradiance를 authoritative cell power와 pixel area로 다시 계산함을 확인했다. 합산 전후 유한성을 검사하며 negative/NaN를 물리적 0으로 바꾸지 않는다. 이는 비간섭 spectral power 합성의 수식을 보존한다. NumPy 정수의 Python 정수 정규화와 JSON 변환은 광학 수식을 바꾸지 않는 범위다. 소프트웨어 계약에 대한 세부 독립 평가는 해당 역할의 리뷰를 따른다.

각 field별 MTF T/S 곡선은 기존에도 존재했고 이번 변경은 실제 계산축을 metadata로 드러낸다. `wave.py:113–126`에서 S는 첫 chief-normal basis, T는 둘째 basis이며, target line 방향은 주파수 방향과 수직이다. Afocal 결과는 동일 frame의 angular chart 좌표축을 나타낸다.

Ansys 공식 FFT MTF와 Sagittal and Tangential 문서를 직접 열어 확인했다. FFT MTF는 object X와 평행한 줄무늬의 응답을 T, object Y와 평행한 줄무늬의 응답을 S로 정의한다. 임의 field 방위각을 따라 고정 label을 자동으로 radial/azimuthal로 돌리는 것은 해당 규약과 다르다. 현재 metadata의 `field_radial_axes=False`는 실제 구현과 일치한다.

다만 현재 projected-global-X로 만드는 직교 chief-normal frame이 모든 2D field에서 제품의 object/pupil mapping과 정확히 같다는 주장은 할 수 없다. `vendor_equivalence='unverified'`를 유지하는 것이 타당하다. 기하학적 meridional T/S와 FFT 고정 T/S를 문서에서 구별해야 한다.

## D1 설명 그림 확인

`d1_explanation/d1_azimuth_explanation_ko.png`와 대응 Markdown을 확인했다. 그림은 렌즈를 회전시키거나 한 field의 T/S 성능을 비교하는 그림이 아니라, 입사 방향의 광축 polar angle \(\theta\)를 고정하고 광축 주위 azimuth \(\varphi\)만 바꾸는 비교다.

\[
a_x=\arctan(\tan\theta\cos\varphi),\qquad
a_y=\arctan(\tan\theta\sin\varphi).
\]

따라서 실제 polar angle \(\theta\)와 \(\sqrt{a_x^2+a_y^2}\)는 일반적으로 다르다. 그림은 10°의 실제 기울기와 7.107076°씩인 대각선 성분 각도를 올바르게 구별하며, 방향 투영을 실제 렌즈 교점으로 표시하지 않는다. D1은 실제 chief 영상 반경의 회전대칭 실패가 아니라 component-wise F-θ reference의 정의 문제라는 설명이 적절하다. 이 리뷰에서 D1 엔진 동작은 변경하지 않았다.

## 재현 및 범위

실행 명령:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 reviews/run_postfix_physics_c.py
```

실행 결과: `status=passed`, detector controls 12개, power controls 12개. 세부 값은 `postfix_physics_c_evidence.json`에 저장했다. 이 수를 root unittest 전체 통과 수에 합산하지 않는다. 전체 suite는 통합 작업자가 별도로 실행한다.

기존 scalar 회절, 단일 경로, 국소 chief-parallel detector 투영, 샘플 수렴 및 underresolved 상태 전파의 미해결 개선점은 이 수정으로 해결되었다고 주장하지 않는다. 상용툴/측정 상관 검증 역시 별도다.

## 직접 확인한 공식 자료

- Ansys OpticStudio 2025 R2, [FFT MTF](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html): field별 FFT MTF, T/S target orientation, image-space frequency 및 FFT 적용 한계.
- Ansys OpticStudio 2025 R2, [Sagittal and Tangential](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Sagittal_and_Tangential.html): 고정 object-local X/Y 규약과 true meridional plane의 차이.
