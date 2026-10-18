# C1–C4 수정 및 D1·field별 MTF 확인 — v0.2.1

2026-10-18. v0.2.0 독립 리뷰의 C1–C4를 수정했다. Cloud 마이그레이션은 보류했다. D1은 정의와 원인을 그림으로 설명하며 계산 코드는 이번에 변경하지 않았다.

## 최종 결과

기존 235개에 이번 회귀 42개(C1 13, C2 7, C3/C4 16, MTF 6)를 더한 **277개 테스트가 모두 통과**했다. 실패·오류·skip 0, 전체 기록 시간 99.641초다. 광학 물리·렌즈 설계·소프트웨어 역할의 별도 AI 에이전트 3명이 구현 후 코드를 독립 검토했다. 후속 리뷰에서 발견한 STOP chain aiming과 complex 파장 입력 문제까지 수정·재검증했다. 각 리뷰의 부분 테스트는 전체 277개와 중복되므로 합산하지 않는다.

## 수정 사항

| ID | 원인과 수정 | 검증 |
|---|---|---|
| C1 | 성공한 교점에서 동일 sag의 다음 면을 다시 풀면 미세 잔차 때문에 정상 광선이 사라졌다. 정확히 동일한 연속 면에서 인증된 교점을 거리 0으로 재사용한다. 동일 STOP chain은 첫 면을 aiming 위치 제약으로 사용한다. 실제 모든 면의 aperture·Snell·TIR는 유지한다. | 기존 217-ray 3 fields, 유한 물체·A2…A12, dummy 전후 STOP, 더 작은 aperture, 실제 추가 굴절, 미소 양의 간격, tangent/first-root. 독립 30-ray chain 최대 위치 차이 6.11e-16 mm. |
| C2 | 작은 polygon 면적을 큰 절대 좌표의 곱의 차로 구해 유효 숫자가 소실됐다. 검출기 원점과 셀 중심을 뺀 좌표에서 clipping하고 국소 triangle cross product로 면적을 구한다. Crop된 전력을 강제 정규화하지 않는다. | 정상 convex cell의 0/5/1000 mm 이동, 부분 crop, 비등방/기울어진 cell, 실제 source-wave 경로. 독립 12개 tilt/azimuth의 구간 적분 대조에서 최대 pixel 질량 오차 1.04e-17. |
| C3 | 파장별 영상 배열의 잘못된 shape가 NumPy broadcasting으로 허용됐다. 격자와 정확한 shape, 실수·유한·비음수 전력/response, 양의 실수 파장 및 합산 overflow를 검사한다. | 비균일 pixel 면적에서 전력/조도 보존, 잘못된 shape·NaN·음수·overflow·complex 파장 거부. |
| C4 | 입력으로 허용한 NumPy 정수가 JSON 저장 시 실패했다. 검증된 count/STOP을 Python int로 정규화하고 공통 JSON 기본형 변환을 사용한다. | 설정 9개 count 및 STOP, 저장/복원 후 실제 evaluate·hash 일치. bool/float는 계속 거부. 독립 정수 8종 및 import 6순서 확인. |

### C2는 geometry가 꼬여서 생긴 문제인가?

**아니다.** 렌즈와 광선 추적이 전혀 없는 대조에서도 정상 convex 직사각형 cell과 정상 detector grid를 함께 평행 이동하면 재현된다. 상대 위치와 면적은 그대로다. 수정 전 5 mm 위치에서는 총 확률 초과 오류, 1000 mm 위치의 특정 작은 셀에서는 전체를 포착해도 0.8614733815가 나왔다. 수정 후 완전 포착은 1이며 부분 crop은 실제 잃은 전력을 유지한다.

잘못된 detector edges·특이 투영을 수용하도록 완화한 수정이 아니다. 사용자가 큰 전역 좌표를 생성할 때 이미 잃은 float 비트까지 복구할 수는 없다. 비축 detector로의 국소 평행 투영이라는 물리 근사도 그대로다.

## D1: 광축 둘레로 입사 방향을 회전한 비교

![D1 방향 좌표 설명](docs/d1_azimuth_explanation_ko.png)

입사 방향과 광축의 실제 각도를 10°로 고정하고 방위각을 90°에서 45°로 바꾼다. X/Y field angle 입력은 A=(0°,10°), B=(7.107076°,7.107076°)다. B의 입력 성분 norm은 10.050923°로 실제 입사각과 다르다. 성분별 F-θ 기준은 이 norm을 사용하므로 방사 방향 F-θ와 일반 2D field에서 같지 않다.

현재 D1에서 실제 chief 영상 반경은 두 방향에서 같았다. 차이는 왜곡을 나누는 기준 좌표에서 생긴다. **한 field의 tangential/sagittal 성능 비교가 아니다.** 상세 식과 원래 0.15 rad 재현 수치는 `docs/D1_explanation_ko.md`, 수치 그림은 PNG/SVG로 제공한다. `metrics.py`는 이전 버전과 동일하며 D1 코드 변경은 하지 않았다.

## Field별 T/S MTF

기존 엔진도 각 field에 대해 T/S를 별도로 계산했다. 이번에는 축 metadata, 명시적인 field 정보, 사용 예제 및 6개 회귀를 추가했다. 렌즈 단독 FFT MTF에 source pixel/fill factor를 곱하지 않는다.

- `mtf_t`, `mtf_s`: 두 MTF 곡선.
- `otf_t`, `otf_s`: 위상까지 보존한 complex OTF.
- `diffraction_mtf_t`, `diffraction_mtf_s`: 같은 통과 pupil·진폭에서 수차 위상을 제거한 기준.
- `frequencies_lp_mm`: focal image-space 주파수. Afocal은 별도 `frequencies_cycles_rad`.
- `axis_convention_id='fixed_chief_normal_xy'`, `mtf_frequency_axes_global`, `mtf_target_line_axes_global`, `field_radial_axes=False`: 계산축의 의미.

Ansys FFT MTF는 object-X로 뻗은 줄무늬의 응답을 T, object-Y 줄무늬의 응답을 S로 정의한다. 주파수 방향은 줄무늬와 수직이다. Field 방위각을 따라 자동 회전하는 radial T/S와 구별한다. 본 엔진은 global X를 chief-normal 면에 투영하여 만든 직교 frame을 사용한다. **임의 2D field에서 Zemax의 pupil mapping 및 수치와 동일하다는 검증은 없다.** `vendor_equivalence='unverified'`를 유지한다.

`examples/run_field_mtf.py`는 축상·Y·X·대각 field의 두 곡선과 회절 기준, 방향 벡터, 품질 상태를 `results/field_mtf_demo.json`에 저장한다. 사용 방법은 `docs/field_mtf_usage_ko.md`다. 이번 예제의 네 field는 모두 `conditional`이며 수렴/적용성 검토를 대체하지 않는다. 90° 회전한 X/Y field에서 고정 축의 T/S가 교환되는 회귀도 포함했다.

## 독립 리뷰와 재현 기록

- 물리: `reviews/postfix_physics_c.md`, 독립 면적/전력 대조 script와 JSON.
- 렌즈: `reviews/postfix_lens_c.md`, 독립 STOP-chain script와 JSON.
- 소프트웨어: `reviews/postfix_software_c.md`.
- 구현 기록: `reviews/fix_c1_implementation.md`, `fix_c2_implementation.md`, `fix_c3_c4_implementation.md`.
- 전체 suite: `results/test_results.json`, `results/test_run.log`.
- 세 예제와 G-only CLI 모두 종료 코드 0. 기록: `results/fix_c_integration.json`.
- 원리 문서 MD/TeX의 수식 231개는 앞뒤 공백을 제외한 AST의 종류·순서·원문이 일치하고 기존 수식도 보존됐다. 기록: `results/core_math_preservation.json`.
- 이전 9개 수정 보고서는 `FIX9_REPORT.md`, 이번 수정 전 리뷰는 `history/CODE_REVIEW_v0.2.0.md`다. 기존 `audit_round2` 역사적 파일은 변경하지 않았다.

## 남는 범위

검토한 수정 범위에서 미해소 차단 결함을 확인하지 않았다는 결과다. 사람 전문가 인증, Maxwell 해 또는 상용툴/측정 상관 검증이 아니다. D1의 방사 방향 정의 변경, 기존 별도 개선사항인 underresolved source kernel의 영상 합성 상태 전파, 무한 물체의 흡수성 object medium 정책 등은 이번 C1–C4 수정에 포함하지 않았다. Scalar·단일 경로·sampled pupil·chief-normal FFT/국소 detector 투영의 한계는 유지된다. `underresolved`를 정답 데이터로 채택하지 않으며 `conditional`도 처방별 수렴 확인이 필요하다.

## 공식 정의 출처

1. [Ansys OpticStudio — FFT MTF](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html)
2. [Ansys OpticStudio — Sagittal and Tangential](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Sagittal_and_Tangential.html)
3. [Ansys OpticStudio — Field Curvature and Distortion](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Field_Curvature_and_Distortion.html)

공식 정의를 확인한 날짜는 2026-10-18이다. 출처를 실행 결과 일치의 근거로 사용하지 않는다.
