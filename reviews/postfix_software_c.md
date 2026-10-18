# C1–C4 수정 후 독립 소프트웨어 리뷰

2026-10-18. 구현자가 아닌 별도 소프트웨어 리뷰 역할의 에이전트가 검토했다. 검토 대상은 루트 `src/sequential_optics/`이며 비교 기준은 v0.2.0 고정 사본이다. 이 리뷰에서는 소스 코드를 수정하지 않았다.

## 판정

**검토한 변경 범위에서 미해결 차단 결함을 발견하지 않았다.** C3의 파장 입력 검증에서 NumPy complex scalar가 실수로 강제 변환되는 잔여 구멍을 찾아 구현자에게 전달했고, 최종 수정 후 거부 동작을 독립 재검증했다. 상용툴 수치 일치나 전체 광학 모델의 정확성 인증을 의미하지 않는다.

## C3 — 파장별 영상 합성

- 공개 입력 경계에서 `cell_power.shape == PixelGrid.shape`를 강제하므로 원래 1×1 전력이 2×2 격자로 broadcasting되던 문제가 차단된다.
- 유한·실수·비음수 전력과 response를 검사한다. 조도는 검증된 전력/물리적 셀 면적으로 다시 계산하여 잘못된 입력 조도 캐시를 사용하지 않는다.
- 입력이 각각 유한해도 합산·response 곱·전체 채널 전력·조도가 overflow할 수 있으므로 결과의 비유한 값을 검사하는 절차가 필요하며 구현되어 있다.
- 양쪽 축 모두 비균일한 2×2 격자에서 두 파장 및 두 response 채널을 독립 구성했다. `irradiance × cell_area == cell_power`와 원본 전력 배열 비변경을 확인했다. 실제 물리적 파장/채널 등록 의미를 보존한다.
- 리뷰 중 `np.complex128(550+20j)`가 `550.0`으로 수용되는 것을 재현했다. 구현자가 `np.iscomplexobj(wavelength_nm)` 거부를 추가한 뒤 complex64, complex128, 허수부가 0인 complex128 모두 `InputError`로 거부하고 경고도 발생하지 않는 것을 재확인했다.

## C4 — NumPy 정수와 JSON

- 허용되는 정수 count와 STOP index를 Python `int`로 정규화하는 수정은 계산 계약을 보존한다. bool/float를 정수로 조용히 받는 변경이 아니다.
- 설정의 `dataclasses.replace`는 전달받은 frozen 하위 설정을 변경하지 않고 정규화된 새 객체를 만든다.
- `system_to_dict`가 JSON 기본형으로 변환된 처방을 반환하므로 `save_system`뿐 아니라 사용자가 직접 `json.dumps(system_to_dict(...))`를 호출하는 경로도 정상화된다.
- 독립 점검에서 NumPy signed/unsigned 8·16·32·64-bit 정수 총 8종을 설정 및 STOP에 적용하고 JSON 왕복을 확인했다.
- `settings`, `models`, `io`의 import 순서 6가지 각각을 별도 Python 프로세스에서 실행하고 `EngineSettings.to_json`을 호출했다. 순환 import/부분 초기화 오류가 없었다. `to_json`의 지역 import는 호출 시점에 수행되므로 현재 의존 관계에서 안전하다.

## C2 및 field별 MTF 변경의 소프트웨어 검토

C2는 PSF 원점과 셀 중심을 뺀 국소 좌표에서 clipping하고, 면적을 국소 삼각형 cross product로 계산하도록 바뀌었다. 작은 셀에 전역 좌표 곱의 차를 적용하던 수치 상쇄를 제거하는 적절한 수정이다. Half-plane 교점은 해당 경계 좌표로 고정하고 출력 합을 강제로 정규화하지 않는다. 물리적 타당성·geometry 여부의 독립 검증은 별도 광학 리뷰가 담당한다.

기존 MTF 배열 계산은 유지하면서 `analyze_wave` 결과에 field와 축 metadata를 추가했다. focal/afocal/source 경로의 기존 T/S 배열과 새 metadata가 일치하며 JSON 변환 가능한 배열·기본형이다. 축 정보는 공간 주파수 방향과 수직인 표적 선 방향을 명시한다. 현재 정의는 **fixed chief-normal X/Y**이며 일반 field-radial/meridional T/S로 이름을 바꾼 수정이 아니다. 해당 제한은 결과의 `field_radial_axes=False`와 `axis_convention_id`로 드러난다.

## 실제 실행한 검증

1. `test_fix_c3_c4` 16개와 `test_settings_io` 5개: **21개 통과**, 0.099 s.
2. `test_mtf_field_axes`: **6개 통과**, 13.756 s. field별 두 곡선·회절 기준, 직교 축, 직사각 pupil 독립 해, X/Y field 회전, afocal, source 경로를 포함한다.
3. 위 별도 import 6순서, NumPy 정수 8종 JSON 왕복, 비균일 2축 격자 전력 보존 점검: 모두 통과.
4. complex 파장 보완 후 해당 회귀 1개 및 별도 complex64/complex128 점검: 모두 통과. 이 회귀는 1번 테스트와 중복되므로 총 테스트 수에 더하지 않는다.

명령은 `PYTHONDONTWRITEBYTECODE=1`, `PYTHONPATH=src:examples:tests`로 실행하여 작업 소스를 사용했다. 전체 suite는 상위 작업자가 별도로 실행한다. 이번 리뷰가 이전 보고서의 별도 범위 항목인 `underresolved` 영상 상태 전파 또는 기타 광학 모델 한계를 해결했다는 뜻은 아니다.
