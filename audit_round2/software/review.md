# 2차 독립 리뷰 — 과학 계산 소프트웨어 개발자

검토일: 2026-10-18. 이 검토자는 배포 코드 작성에 참여하지 않았다. 먼저 구현과 테스트를 읽고 재현을 수행한 뒤, 마지막에 기존 `reviews/software_review.md`의 범위와 비교했다. 배포 코드와 기존 결과 파일은 수정하지 않았다.

## 판정

모듈 분리, 기본 G 파장, 필요한 자료만 요구하는 물성 호출, JSON 처방 재현성과 명시적 오류 분류는 좋은 기반이다. 그러나 **설정이 고수준 API에서 일부 사라지는 실제 결함**, **영상 데이터 객체의 생성 후 불변식이 보장되지 않는 결함**이 남아 있다. 잘못된 thermal 입력을 조기에 차단하지 못하는 낮은 우선순위 결함도 재현했다.

전체 자동 테스트가 통과했다는 사실과 이 결함들은 양립한다. 이번 반례는 기존 테스트의 범위를 벗어난다. 이 보고서는 수식 전체나 Zemax 일치를 인증하지 않으며, 다른 역할의 광학 물리/렌즈 설계 검토와 합쳐 판단해야 한다.

## 검토 대상과 재현성

- 배포 ZIP SHA-256: `d854cbaba525581088fe7b890ee54e2f1fb663c8ef8c6fb562b06a66a65284ef`
- 코드: `../input/optical_engine/src/sequential_optics/`.
- 새 결과와 재현 스크립트는 이 `software/` 폴더에만 저장했다.
- 기존 editable 설치를 사용하지 않도록 frozen source 경로를 `PYTHONPATH` 및 스크립트 `sys.path[0]`로 명시했다. 실제 import 경로도 결과 JSON에 기록했다.
- Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0.
- 배포 테스트 재실행: **183개, 실패/오류 0개, 39.735초**. 원본 `run_tests.py`는 기존 결과를 덮어쓰므로 사용하지 않았다. 로그: `baseline_tests.log`.
- 추가 재현: `python reproduce_software_findings.py --output reproduction_results.json`. 현재 배포의 세 결함과 선택적 데이터 정상 동작을 모두 assertion으로 확인하며 약 1초 걸린다.

## SW-01 — 고수준 API가 spot/refocus의 사용자 설정을 누락

**우선순위 P2 / 확신 높음 / 실제 재현 완료.** 일반적인 광축 부근 기본 설정 결과가 모두 잘못되었다는 뜻은 아니다. 사용자 설정을 바꿨을 때 저수준 API와 고수준 API가 서로 다른 판정을 내리는 결함이다.

위치:

- `src/sequential_optics/api.py:117`: `gm.spot_diagram(current,bundle)`에 `settings=settings` 누락.
- `src/sequential_optics/api.py:139`: `gm.common_geometric_refocus(current,bundles)`도 동일.
- `src/sequential_optics/metrics.py:119`, `679`: 설정이 없으면 새 `EngineSettings()` 생성.
- `src/sequential_optics/metrics.py:86–95`: 이 설정의 `direction_tol`이 영상면 평행 판정에 쓰인다.

최소 재현 조건:

1. 굴절 없는 평면 aperture 한 장, 반경 1 mm, IMA z=0.001 mm.
2. field X=`π/2−1e−13` rad. 극단적인 grazing 수치 경계 사례이지만 선언된 forward field 범위 안에 있다.
3. `direction_tol=1e−15`, boundary clearance=`1e−12 mm`, rings=2, arms=8, boundary rays 없음.
4. 같은 처방·광선·설정으로 직접 `spot_diagram(...,settings=cfg)`와 `evaluate(...,settings=cfg)`를 비교.

| 항목 | 관찰값 |
|---|---|
| 최종 chief z 방향 성분 | `9.998130455622147e−14` |
| 실광선 성공 | 16/16 |
| 직접 spot | `ok` |
| 고수준 spot | `not_applicable`, `Ray is parallel to the image plane` |
| 직접 common focus | `no_unique_finite_focus`: 평행 광선의 올바른 적용성 진단 |
| 고수준 common focus | `not_applicable`, `Ray is parallel to the image plane` |

원인은 고수준 호출에서 기본 `direction_tol=1e−12`로 되돌아가기 때문이다. 결과 파일에는 사용자의 `1e−15` 설정이 기록되므로 기록된 계산 조건과 실제 판정 조건까지 불일치한다.

수정 제안: 두 호출에 설정을 전달하고, 모든 orchestration 경로의 설정 전달을 계약으로 검사한다. 물리량에 영향을 주는 사용자 설정을 바꾼 뒤 직접 모듈 호출과 batch 호출의 상태·값을 비교하는 회귀 테스트를 추가한다. 이 사례를 수정하려고 광선 조건이나 허용 범위를 임의로 바꿀 필요는 없다.

## SW-02 — frozen 영상 객체의 NumPy 배열을 바꾸면 보존 조건 우회

**우선순위 P2 / 확신 높음 / 실제 재현 완료.** 이 문제는 **호출자가 생성된 객체의 배열을 in-place 수정한 경우** 발생한다. 정상 실행 중 저절로 광량이 증가한다는 주장은 아니다.

위치:

- `src/sequential_optics/imaging.py:17–27`: `PixelGrid`는 frozen dataclass지만 내부 edge 배열은 writable copy.
- `src/sequential_optics/imaging.py:39–58`: `SourceKernel`의 확률 질량 배열도 writable copy.
- `src/sequential_optics/imaging.py:82–89`: 합성 때 grid 동일성과 파장만 확인하며 kernel mass/면적의 불변식을 재검사하지 않는다.

최소 재현:

```python
grid = PixelGrid([0., 1.], [0., 1.])
kernel = SourceKernel(grid, [[1.]], 550.)
kernel.probability_mass[0, 0] = 2.       # 예외 없이 성공
r = synthesize_image([kernel], [1.], [1.])
# collected_power=1, detected_roi_power=2, outside_roi_power=0

grid.x_edges_mm[1] = -1.              # 단조 증가 조건 붕괴; 예외 없음
r = synthesize_image([kernel], [1.], [1.])
# irradiance = [[-2.]]
```

기대 동작: immutable 값 객체로 제공한다면 수정이 거부되어야 한다. 가변 객체를 허용한다면 계산 경계에서 재검증하여 음수 면적 또는 합계 1 초과 확률을 거부해야 한다.

영향: 영상 커널을 재사용하면서 스케일/정규화 또는 좌표 변환을 in-place로 하는 통상적인 NumPy 사용자 코드가, 생성 시 검증을 우회할 수 있다. 원래 전력 초과를 드러내야 할 `outside_roi_power`가 `max(0,...)`로 0이 되어 전력 ledger만 확인해서 발견하기도 어렵다.

수정 제안: `Ray` 객체에서 이미 사용한 방식처럼 방어적 copy에 `setflags(write=False)`를 적용한다. 확장 API에서 외부 커널을 받는 경우에는 합성 경계에서도 비음수/유한/총확률/양의 면적을 확인한다. 새로운 커널은 생성자를 거쳐 만드는 계약을 문서화한다. 정상 커널의 전력 합성식 자체는 이 재현에서 문제가 없었다.

## SW-03 — 잘못된 thermal surface index가 raw TypeError로 batch 중단

**우선순위 P3 / 확신 높음 / 실제 재현 완료.** 물성 데이터가 생략되었을 때의 결함이 아니라 **잘못된 입력 자료형 검증의 결함**이다.

위치:

- `src/sequential_optics/thermal.py:30–36`: `LensGroup`가 앞/뒤 index의 정수형 여부를 확인하지 않는다.
- `src/sequential_optics/thermal.py:71`: 검증되지 않은 index로 tuple 접근.
- `src/sequential_optics/api.py:45–53`: `_attempt`가 `OpticsError`만 구조화하므로 이 TypeError는 밖으로 전파된다.

최소 입력:

```python
RunRequest(metrics=('first_order',), temperatures_c=(20.,), fail_fast=False,
    thermal_options={'lens_groups': [{'front_index': 0.0, 'back_index': 1.0}]})
```

요청 생성은 성공한다. 이후 `evaluate`는 `TypeError: tuple indices must be integers or slices, not float`로 종료한다. `LensGroup(False, True)`도 수치 비교를 통과하여 index로 사용될 수 있다.

기대 동작: 정수형 index만 허용한다면 생성 시 명확한 `InputError`로 거부한다. 정수값인 실수를 허용할 설계라면 명시적으로 정규화해야 한다. `OpticalSystem.stop_index`에 이미 적용한 `Integral` 검사와 bool 배제 방식을 일관되게 적용할 수 있다.

수정 제안: 입력 스키마를 보완한다. 모든 Python 예외를 무차별적으로 catch해서 물리 실패로 바꾸는 대응은 권장하지 않는다. 내부 프로그래밍 오류와 자료 부족을 구별하는 현재 오류 계층은 유지하는 것이 맞다.

## 잘 구현되어 있음을 직접 확인한 항목

새로운 positive control에서 단일 G 굴절률만 있는 singlet에 대해 다음을 확인했다.

| 자료/요청 | 결과 |
|---|---|
| CTE·dn/dT·흡수 없음, nominal G first-order/spot | 모두 `ok` |
| 같은 요청에 B 추가, B 굴절률 없음 | G 보존, B만 `missing_data` |
| 공통 20°C 기준 온도 명시, CTE 없음 | G 계산 가능; 불필요한 팽창 자료 요구 안 함 |
| 40°C 요청, CTE 없음 | 해당 온도만 `missing_data`; 기준 결과 보존 |

코드에서 확인한 좋은 설계:

- 단일 파장 측정 표와 명시적인 무분산 `n_constant`를 구별한다. Abbe 수나 재료명만으로 임의 분산을 만들지 않는다.
- CTE의 instantaneous/reference-derivative 정의와 기준 형상 온도를 나눠 처리한다. 기준 형상으로부터 온도별 상태를 생성하고 hot geometry의 누적 재팽창을 거부한다.
- 열 형상 생성과 n(T)/흡수 호출이 분리되어 있어 필요한 분석에서만 누락 자료를 검사한다.
- `FrozenDict`는 caller가 원래 설정/material registry dictionary를 바꾸는 단순 alias 문제를 막는다. Ray 배열도 읽기 전용이다. SW-02는 다른 데이터 객체에서 아직 같은 방어가 빠진 사례다.
- 복소 OTF를 실수로 버리지 않고 real/imag를 보존하며 nonfinite 결과는 JSON에 명시적으로 표시한다.
- 주요 batch API는 수치/적용성/자료 부족 오류를 구분한다. `status=ok`는 실행 성공이라는 문서 규약이며, 파면 `conditional`을 검증 성공으로 몰래 바꾸지 않는다.

## 선언된 한계 및 개선 권고 — 확인된 구현 결함과 구별

- `input_sha256`는 처방·설정·요청 identity이며 코드 revision identity가 아니다. 코드 버전/배포 ZIP 해시를 함께 저장하는 것이 데이터셋 재현성에 필요하다. 현재 README와 이전 리뷰가 이 의미를 설명하므로 새 버그로 세지 않았다.
- 원시 `state.status='ok'` 또는 `fr.rays.status='ok'`만으로 학습 데이터에 편입하면 안 된다. nested analysis status와 수렴 지표를 해석해야 한다는 현재 계약을 유지하되, 응용 계층에 명시적인 데이터 채택 정책을 두는 편이 좋다.
- RI/source-image의 일부 sampling 인자는 함수 keyword 기본값으로 남아 있다. 사용자에게 중요한 재현 설정을 하나의 RunSettings/요청 스키마에 모으는 것은 확장성 개선 사항이다. 이를 수식 오류라고 판단하지 않았다.
- 단일 sample 수렴/테스트 통과는 상용툴 일치, 실제 제조 물성, off-axis scalar 근사의 정확도를 보장하지 않는다. 기존 명세에 선언된 무효부 형상·scalar/forward 모델 범위를 확대해 해석하면 안 된다.
- 다음 테스트 우선순위는 기본값 반복 검증보다 (1) custom-settings API 동등성, (2) 생성 후 배열 불변식, (3) 실제 입력 스키마별 경계 타입, (4) traced bundle과 처방/온도 provenance 일치다. 마지막 항목은 별도 광학 물리 리뷰와 협의했으며 여기서 중복 결함으로 등록하지 않았다.

## 종합 권고

현 배포는 해석 연구용 프로토타입으로 실행·조합할 수 있지만, 검증된 범위를 벗어난 자동 데이터 생산에 바로 신뢰 표시를 붙일 단계는 아니다. 우선 다른 도메인 리뷰의 물리 결함과 함께 SW-01을 수정하고, SW-02처럼 재사용 객체의 조건을 무너뜨리는 경로를 막는 것이 합리적이다. 이 검토는 코드 수정 없이 문제와 재현만 남긴다.
