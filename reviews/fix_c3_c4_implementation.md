# C3/C4 구현·검증 기록

2026-10-18. v0.2.0 독립 소프트웨어 리뷰의 C3/C4를 수정한 구현자 기록이다. 원래 리뷰 증거 및 역사적 테스트 결과를 변경하지 않았다.

## C3 — 파장별 영상 합성 경계의 검증

`registered_spectral_sum`은 외부에서 전달되는 각 영상의 `grid`, `cell_power`, `wavelength_nm`를 검사한다. `cell_power.shape`는 `PixelGrid.shape`와 **정확히** 같아야 한다. 유한한 실수·비음수 전력 및 양의 유한 파장이 필요하다. 이를 통해 2×2 물리 격자에 1×1 전력이 들어와 NumPy broadcasting으로 네 셀의 조도로 복제되던 문제를 막는다.

파장별 전력이 이미 광원의 스펙트럼을 포함한다는 기존 계약은 유지한다. 선택적 response의 크기는 `(channels, wavelengths)`이며 적어도 한 채널이 있어야 한다. 이 response도 실수·유한·비음수여야 한다. 합산·response 곱·채널별 전체 전력 및 `cell_power / cell_area`의 비유한 결과는 `InputError`로 거부한다. 원본 입력 전력과 response가 각각 유한하더라도 합산 또는 곱이 overflow할 수 있다는 리뷰 반례를 포함한다.

계산 기준량은 `cell_power`다. 입력 dictionary의 캐시된 irradiance 값을 신뢰하지 않고 검증된 면적으로 다시 계산한다. 다른 파장의 lateral color를 재중심화하거나 PSF tail을 다시 정규화하지 않는다. 비균일 픽셀 면적의 양성 대조에서 두 입력의 검출 전력은 `[4.5, 1.5]`, `[1, 4]`, 합은 `[5.5, 5.5]`다. 픽셀 면적 `[2, 4] mm²`에서 조도 `[2.75, 1.375]`를 얻고 조도×면적의 합은 검출 전력 합인 11과 일치한다. 별도 3채널 response도 예상된 선형 합과 일치한다.

## C4 — 허용된 NumPy 수치 입력의 JSON 왕복

엔진이 이미 허용하던 NumPy 정수 scalar count 및 STOP index를 검증한 후 Python `int`로 정규화한다. 설정 내부의 3개 ray count와 6개 sampling count에 공통 적용한다. 원래 사용자가 전달한 frozen 설정 객체는 변경하지 않고 `dataclasses.replace`로 엔진의 새 설정을 생성한다. bool·NumPy bool·float·NumPy float index/count는 기존과 같이 거부한다.

`EngineSettings.to_json`과 `system_to_dict`는 기존의 공통 `json_value` 변환기를 사용해 그 밖의 허용된 NumPy scalar도 JSON 기본형으로 변환한다. 후자는 공개 `system_to_dict`의 반환값을 직접 `json.dumps`에 넘기는 경로도 지원한다. `settings.to_json`의 로컬 import는 현재 `io -> settings` 의존에서 순환 초기화를 피하기 위한 것이다.

NumPy 정수 STOP 및 NumPy scalar 처방을 저장하고 복원한 뒤 실제 `evaluate(..., metrics=('first_order',))`의 성공을 확인했다. 저장 전후의 입력 SHA256도 일치한다. 설정의 9개 count 모두와 NumPy 정수 G 파장도 JSON 왕복 후 같은 설정을 유지한다.

## 실제 실행 결과

- 신규 `tests/test_fix_c3_c4.py`: **16개 통과**, 실패·오류·skip 0. 0.052 s.
- 기존 관련 회귀 `test_settings_io`, `test_review_software`, `test_integration`, `test_radiometry`: **43개 통과**, 실패·오류·skip 0. 11.912 s.
- 실행 환경에서는 `PYTHONPATH=src:examples[:tests]`, `PYTHONDONTWRITEBYTECODE=1`을 명시해 이번 작업 소스를 실행했다.
- 전체 통합 테스트 및 구현자와 다른 역할의 독립 재리뷰는 상위 작업자가 수행한다. 이 기록의 59개는 전체 테스트 수에 추가해 이중 계산하면 안 된다.

독립 소프트웨어 재리뷰에서 NumPy complex 파장을 `float()`가 허수부를 버리고 받는 작은 경계 누락을 추가 확인했다. 파장도 변환 전에 복소수 여부를 검사하도록 수정했고 `np.complex64`와 `np.complex128` 반례를 기존 파장 검증 테스트에 추가했다. 이 최종 경계 수정 후 같은 신규 회귀 **16개를 0.047 s에 재실행해 모두 통과**했다.

## 변경 범위와 제한

수정 파일은 `imaging.py`의 spectral composition 경계, `settings.py`, `models.py`의 STOP 정규화, `io.py`의 처방 변환과 신규 테스트다. 광선 교점·회절·MTF·열해석 계산식은 이 작업에서 변경하지 않았다. C3는 물리 모델의 새로운 정확도 보증이 아니라 잘못된 입력과 수치 overflow를 정상 물리 출력으로 반환하지 않도록 하는 수정이다. `underresolved` 상태 전파, 극단적 complex OTF 상대 가중치 등 원래 리뷰에서 별도 개선으로 분류한 항목은 이 수정에 포함하지 않았다.
