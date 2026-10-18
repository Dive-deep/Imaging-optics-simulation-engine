# MTF T/S 규약 확인 — 독립 광학 자문

검토 범위: `wave.py`, `source_wave.py`, `metrics.py`, `api.py`, `docs/core_principles_ko.md`, `docs/calculation_spec_ko.md`. 엔진 소스는 수정하지 않았다.

## 결론

**현재 엔진은 이미 각 field에서 `mtf_t`, `mtf_s` 및 동일 pupil의 `diffraction_mtf_t`, `diffraction_mtf_s`를 계산한다.** `evaluate()`의 `states → wavelengths → fields → mtf → value`로 접근한다. Focal 주파수 배열은 `frequencies_lp_mm`, afocal은 `frequencies_cycles_rad`다.

사용자가 Zemax FFT MTF를 기준으로 요구했으므로, **기존 고정 T/S 규약을 자동 field-radial 규약으로 바꾸면 안 된다.** Ansys 공식 FFT MTF 문서는 T를 object X에 평행한 줄무늬, S를 object Y에 평행한 줄무늬의 응답으로 정의한다. 줄무늬의 연장 방향과 공간주파수 벡터 방향은 수직이다. 통상 field가 YZ면에 있을 때 T의 주파수는 Y 방향, S의 주파수는 X 방향이다.

Ansys의 별도 Sagittal and Tangential 정의는 field 위치와 관계없이 object-local Y/X의 고정 규약을 쓴다고 명시한다. 따라서 X/Y 양쪽 field까지 임의로 회전하면서 T/S 이름을 자동으로 교환하는 동작은 제품 관례의 보완이 아니라 **다른 규약을 추가**하는 것이다.

## 현재 구현의 정확한 뜻과 한계

최종 chief 방향을 단위벡터 $\mathbf d_c$라 하면, 일반적인 유효 분기에서

$$
\mathbf e_1=
\frac{\mathbf e_X-(\mathbf e_X\cdot\mathbf d_c)\mathbf d_c}
 {\|\mathbf e_X-(\mathbf e_X\cdot\mathbf d_c)\mathbf d_c\|},
\qquad \mathbf e_2=\mathbf d_c\times\mathbf e_1.
$$

특이한 경우 global Y를 투영하는 기존 fallback이 있다. `output_basis=[e1,e2,dc]`이며

$$
\mathrm{MTF}_S(f)=|H(f,0)|,\qquad
\mathrm{MTF}_T(f)=|H(0,f)|.
$$

각 주파수 벡터는 각각 $f\mathbf e_1$, $f\mathbf e_2$이다. 해당 계산면에서 S target 줄무늬는 $\mathbf e_2$, T target 줄무늬는 $\mathbf e_1$와 평행하다. `source_wave.py`도 같은 축을 쓴다. Afocal 역시 동일 frame의 정확한 각도 $\alpha=\operatorname{atan2}(l,\gamma)$, $\beta=\operatorname{atan2}(m,\gamma)$에서 S/T를 평가하며, rad 단위와 기존 sampling guard를 보존해야 한다.

이 구현은 독립 scalar chief-normal profile이다. **모든 2D field에서 제품의 pupil orientation, working F-number scaling, object-space 축 mapping과 수치적으로 동일하다는 보증은 없다.** 특히 $d_{c,x}d_{c,y}\ne0$이면 global X와 Y를 chief-normal 평면에 각각 투영한 벡터는 일반적으로 직교하지 않는다. 현재 코드는 projected X 하나로 직교 frame을 만든다. 따라서 `vendor_equivalence='unverified'`를 유지하고, fixed-axis 이름의 일치와 전체 제품 수치 일치를 구별한다. 현재 계산 명세도 이 제한을 이미 선언한다.

## 최소 변경 권고

1. 현재 `mtf_s/t` 값을 유지한다. field별 결과 사용 예시와 metadata를 명확하게 추가한다.
2. `axis_convention_id='fixed_chief_normal_xy'`처럼 식별자를 제공한다.
3. `frequency_axis_S_xyz=output_basis[0]`, `frequency_axis_T_xyz=output_basis[1]`를 제공하거나 이에 상응하는 묶음을 저장한다. 이 벡터는 **주파수 증가 방향**이라고 명시한다.
4. `target_line_axis_S_xyz=output_basis[1]`, `target_line_axis_T_xyz=output_basis[0]`는 선택적 편의 정보다. 물리적으로는 계산면의 국소 target 방향이다. 임의 2D field에서 실제 object plane 방향과 정확히 같다고 이름 붙이지 않는다.
5. `metrics.py`의 spot/ray-fan/field-curvature T/S는 field meridian을 따라가는 별도 규약임을 사용자 문서에서 구분한다. 같은 글자 T/S가 분석 종류 사이의 동일한 3-D 축을 보장하지 않는다.

## 회귀 검증 권고

- `RunRequest`에 on-axis, Y-field, X-field를 함께 넣고 **각 field 결과에** 길이가 같은 frequency/MTF T/S/diffraction T/S 배열이 존재함을 확인한다. DC는 1이다.
- 원형 무수차 pupil에서는 T/S가 일치하고 해석적 원형 MTF와 수렴한다.
- 축 정렬 비등방 pupil(직사각 aperture 또는 한 방향 quadratic phase)을 사용해 T/S가 서로 달라지는 사례를 확인한다. $H(f,0)$와 S, $H(0,f)$와 T가 일치함을 비교하면 잘못된 label swap을 잡는다.
- 90° 회전시킨 동일 비등방 pupil의 두 곡선이 교환됨을 확인한다. 이는 **고정 축의 의도된 행동**이다. 모든 azimuth에서 같은 label의 MTF가 불변이어야 한다는 검사는 잘못이다.
- 출력 `output_basis`의 정규직교성, 오른손성 및 chief-normal 성질을 검사한다. 축 메타데이터가 실제 slice 방향과 같은지 확인한다.
- Afocal의 높은 주파수에 기존 R1 sampling guard가 계속 적용되어 `underresolved`/NaN을 보존하는지 확인한다. field별 표를 만들면서 NaN을 0으로 바꾸면 안 된다.

## D1과의 구별

D1의 회전은 입사 방향의 polar angle $\theta$를 유지한 채 **광축 z 주위의 azimuth $\phi$를 변경**한 것이다.

$$
\mathbf d=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta),
\quad \theta_x=\arctan(\tan\theta\cos\phi),
\quad \theta_y=\arctan(\tan\theta\sin\phi).
$$

같은 $\theta$를 갖는 field의 chief 방향이 원뿔 둘레를 도는 비교다. 하나의 field에서 주파수/target 방향을 T와 S로 바꾸는 비교가 아니다. 회전대칭 계에서 **방사 방향의 F-θ 왜곡 정의**라면 같은 polar angle에서 방위각에 무관해야 한다. 하지만 component-wise $(f\theta_x,f\theta_y)$를 기준 좌표로 쓰면 $\sqrt{\theta_x^2+\theta_y^2}$가 $\theta$와 같지 않아 azimuth 의존값이 생긴다. 실제 ray tracing 자체의 비대칭을 뜻하지 않는다.

## 확인한 공식 자료

1. Ansys OpticStudio 2025 R2, **FFT MTF**. Field별 FFT 분석, T/S target orientation, image-space 주파수 및 sampling 한계.
   https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html
2. Ansys OpticStudio 2025 R2, **Sagittal and Tangential**. True meridional plane과 제품 fixed-axis 규약, TAN으로 MTF 규약이 재정의되지 않는 점.
   https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Sagittal_and_Tangential.html
3. Ansys, Mark Nicholson, **Why are FFT and Huygens MTF results different on tilted image surfaces?** FFT의 chief-normal 계산면과 Huygens 계산면 차이.
   https://optics.ansys.com/hc/en-us/articles/42661833682707-Why-are-FFT-and-Huygens-MTF-results-different-on-tilted-image-surfaces

자료 확인일: 2026-10-18. 문서의 설명을 요약했고 긴 원문 인용은 하지 않았다.
