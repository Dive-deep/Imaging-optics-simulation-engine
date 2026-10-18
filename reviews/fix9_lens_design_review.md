# R1–R9 수정 후 독립 렌즈 설계 관점 리뷰

검토일: 2026-10-18. 역할: 카메라·투사 결상 광학 설계 엔지니어 관점의 독립 AI 리뷰. 코드 작성자의 설명만을 근거로 승인하지 않고, 수정 전 `audit_round2/input/optical_engine`와 작업 중 루트 소스를 비교하고 독립 수치 대조를 실행했다. 본 검토자는 엔진 코드를 수정하지 않았다.

## 판정과 범위

**최종 확인한 R3–R6 수정은 원래 결함의 물리·계산 의도를 만족한다. 검토 범위에서 미해결 차단 결함은 발견하지 못했다.** 초기 수정본에서 추가로 발견한 곡면 STOP rim 허용오차 불일치는 작성자에게 전달하여 수정되었고, 같은 독립 대조를 다시 실행해 해소를 확인했다.

핵심 범위는 STOP aiming/실제 aperture 분리, 곡면 교점, LSA 기준점, spot 표본/적분 가중치 구분이다. EFL/BFL, uniform thermal scaling, 공통 geometric focus 및 새 원리 문서도 대조했다. R1/R2와 광량·회절 전체의 승인 판단은 해당 전문 역할의 별도 리뷰를 함께 읽어야 한다. 테스트 통과는 Zemax 상관 검증 또는 실제 렌즈 성능 보증이 아니다.

## 수정 내용의 타당성

| 항목 | 수정 전 문제 | 최종 코드 판단 |
|---|---|---|
| R3 | STOP 앞의 유한 aperture 때문에 aiming 잔차 자체가 정의되지 않아 정상 clipping을 수치 실패로 분류 | 가상 연장한 analytic sag의 real-ray 잔차로 발사 광선을 찾고, 별도 `trace_ray`가 실제 aperture를 적용한다. 곡면 STOP의 실제 sag 높이를 목표로 사용하며, 마지막에 STOP 첫 교점을 다시 검증한다. 실제 surface를 근축으로 바꾸지 않는다. |
| R4 | 보조 cylinder 판별식의 미소 음수가 유효 rim을 누락 | 항의 크기에 따른 roundoff bound와 안정적인 이차근, 평면 직접 해를 사용한다. 추가 확인 과정에서 곡면 cylinder와 aperture의 길이 허용오차도 일치시켰다. 실제 optical tangency와 진짜 miss는 계속 거부한다. |
| R5 | 사용자가 바꾼 LSA 보고 기준을 q=0의 물리적 축 교점으로 사용 | `paraxial_focus_z`와 `reference_focus_z`를 분리했다. 기준면 변경은 LSA에 상수 offset만 준다. 유한 물체의 paraxial conjugate도 올바른 기준이다. |
| R6 | 0 적분 가중치의 chief/rim이 plot와 최대 반경에서 빠짐 | 성공한 모든 기하 표본을 plot/최대 반경에 포함한다. centroid/RMS/common focus에는 양의 적분 가중치만 사용하여 물리적 적분 측도를 바꾸지 않는다. Angular spot에도 같은 구분이 적용된다. |

이 방식에서 `escaped`는 지정된 유한 sequential surface를 통과하지 못한 광선이다. 형상 정의 자체가 없거나 local aiming branch를 찾지 못한 상태를 임의로 광량 0으로 바꾸지는 않는다. 연속 pupil 전체의 극대 spot 반경을 보증하는 기능으로 확대 해석해서도 안 된다.

## 이 리뷰에서 새로 발견하여 해결한 곡면 rim 문제

독립 처방은 첫 면 z=0의 평면, clear radius 2 mm, air→n=1.5, 두 번째 면인 STOP은 z=10 mm, curvature 0.02/mm, radius 2 mm, n=1.5→air다. 입력 field는 X=0.2 rad다. pupil=(1,0)의 광선은 평면 Snell의 독립 해에서 첫 면 x=0.658418753 mm로 aperture 안에 있다.

초기 수정본은 가상 STOP 교점을 r=2+5.4×10⁻¹³ mm에서 얻어 aiming tolerance를 통과시켰으나, physical cylinder를 정확히 r≤2에서 잘랐다. 실제 근 직전 sag residual이 −1.15×10⁻¹² mm로 남으면서 `escaped`가 되었다. 이 상태는 실제 비네팅이 아니므로 수정자와 주 검토자에게 즉시 전달했다.

최종 수정은 bounding interval의 radial 길이 허용오차를 실제 aperture predicate와 맞췄다. 동일 독립 사례 및 추가된 `test_refracting_curved_stop_rim_uses_same_physical_radial_tolerance`가 통과한다. 실제 optical tangent를 단순히 허용하는 수정은 하지 않았다.

## 직접 실행한 검증

모든 실행은 저장소 루트에서 `PYTHONPATH=src`로 수행했다. 실제 import 경로가 `/workspace/scratch/e8215b7f31f6/optics_core_v0_2/src/sequential_optics/__init__.py`임을 출력 확인했다. 고정 과거 소스가 아닌 수정 소스를 실행한 결과다.

- `test_regression_r3_r4.py`: **13개 통과**. 평면 겹침, 굴절 sphere의 독립 clipping threshold, 곡면 STOP rim, 실제 tangent, later-root 방지, 정의되지 않은 continuation의 실패 구분을 포함한다.
- `test_regression_r5_r6.py`: **9개 통과**. LSA 기준 이동, 유한 물체의 평면 Snell 가상 공액, spot/angular plot와 quadrature 분리, 공통 초점의 가중치 유지 등을 포함한다.
- 위 두 모듈의 개수는 전체 suite 개수와 별도다. 본 리뷰에서 전체 suite를 다시 실행했다고 표시하지 않는다.

테스트 파일과 별개로 수행한 독립 대조:

| 대조 | 독립 근거 및 결과 |
|---|---|
| 평면 굴절→곡면 STOP | Snell에서 sinθ′=sin(0.2)/1.5를 구하고 STOP의 실제 sag 높이에서 역산한 첫 면 위치와 비교. 12 azimuth×3 radius=36광선, 최종 22통과/14물리 miss 모두 독립 aperture 판정과 일치. 발사 XY 최대 오차 1.54×10⁻¹³ mm. |
| 곡면의 보조 cylinder 접선 | sphere rim을 정확히 지나는 직선을 독립 구성. 길이 배율 0.001, 1, 1000에서 hit 오차 각각 6.62×10⁻¹⁷, 1.78×10⁻¹⁵, 0 mm. 보조 cylinder 접선이 실제 optical surface 접선과 다름을 확인. |
| 두꺼운 렌즈 | R₁=50, R₂=−50, thickness=5 mm, n=1.5의 lensmaker/ABCD 독립 해와 일치. EFL=50.847457627118644 mm, BFL=49.152542372881356 mm. |
| 굴절률 고정 균일 팽창 | CTE=10⁻⁴/K, ΔT=100 K에서 길이 배율 1.01. EFL 비율도 정확히 1.01. 이는 n(T)와 형상 변화를 분리한 similarity 대조다. |
| 실제 triplet의 R6 | 8 rings/24 arms, 성공 표본 217개 중 양의 가중치 192개. 독립 image-plane intercept로 centroid/RMS/최대를 다시 계산. RMS=0.004556286083508109 mm, sampled max=0.009032283705502662 mm와 일치. |
| 실제 triplet의 R5 | 보고 기준을 +7.3 mm 이동해도 축 교점 배열은 그대로이고 LSA만 −7.3 mm 이동. q=0은 −7.299999999999997 mm. |
| 공통 재초점 | 축상/0.03 rad 두 field의 가중치 1:2. 공통 z=91.73896769696545 mm에서 직접 재계산한 RMS² merit 1.1746183121×10⁻⁵ mm². ±0.1 mm에서는 각각 2.4532279094×10⁻⁵ mm²로 증가. 한 센서면의 공통 최소점이라는 정의와 일치. |

## 핵심 원리 문서의 정의 검토

`docs/core_principles_ko.md`의 다음 구분은 코드와 일치하고 구현 참고 자료로 적절하다.

1. 실제 sag 교점·vector Snell이 spot과 mapping의 기본이며, 근축 행렬은 EFL/BFL·초기 aiming에 별도로 쓰인다. A₂가 vertex curvature에 2A₂로 들어가는 점도 맞다.
2. BFL 원점은 마지막 물리면 vertex다. 유한 공액의 `paraxial_image_z`를 무한 공액 BFL과 구별한다.
3. Spot centroid/RMS의 pupil 측도와 geometric maximum의 표본 집합을 구별하고, LSA 물리 교점과 보고 기준면을 분리한다.
4. Distortion은 실제 chief mapping의 작은 field 미분을 기준으로 하며, field curvature는 chief 주변 parabasal 미분 초점이다. Full-pupil 최적 초점과 같은 지표로 부르지 않는다.
5. FFT의 chief-normal 계산면과 실제 global-z 센서면을 구별하며, 이 구현의 fixed X/Y FFT 축과 field-meridional 기하 T/S를 혼동하지 않는다.
6. Uniform thermal state는 기준 처방에서 새로 만들고, 고정 IMA와 모든 field의 공통 재초점을 구분한다. 원리 문서는 예외·버그 기록으로 흐르지 않고 기본 계산 흐름을 다룬다.

관련 Ansys 공식 문서를 직접 열어 확인했다. `Field Curvature and Distortion`의 parabasal Z-coordinate 기준 및 실제 small-field chief 기준, `FFT PSF`의 chief-normal 계산면, `FFT MTF`의 image-space frequency와 별도 FFT T/S 규약이 문서의 설명을 뒷받침한다. 다만 자체 pupil mapping과 수송 모델까지 vendor 동일성으로 확대하지 않은 현재 표현을 유지해야 한다.

- https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Field_Curvature_and_Distortion.html
- https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_PSF.html
- https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html

## 검토한 소스 식별자

다음은 위 최종 대조 시점의 SHA-256이다. 이후 문서 편집이 있으면 문서 hash는 변경될 수 있다.

```text
geometry.py  93fef04a5f27d69f70bcd41b99fc477d561b697942d146c5c05f605e9df296fe
metrics.py   1d9dde79413f65a6c10700e9773318a37bcd2ce7f9a6ef3b7a4be21f22f33659
paraxial.py  3e2f885a04bb0b1cb41e8e184533a464c2cf593a585cd312d9a89068c5247248
thermal.py   a173e504be65ebee4c1f3f90b49c0660853943a0683a87a7e21ff7c291f637b0
core_principles_ko.md  491bf5bf0e3187f4d0351032608d417bb6a0cc802f7a0b4ccbb923de804f608a
```
