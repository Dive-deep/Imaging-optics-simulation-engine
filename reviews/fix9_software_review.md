# R1–R9 수정 후 독립 소프트웨어 리뷰

검토일: 2026-10-18. 대상은 `optics_core_v0_2/src/sequential_optics/` 작업 소스다. 고정 과거 사본 `audit_round2/input/optical_engine/src/`와 변경 내용을 비교했다. 구현 담당자와 별도인 소프트웨어 검토 에이전트가 코드·설정 전달·상태 무결성·전력 불변식·오류 격리·직렬화를 검토했다. 이 검토자는 엔진 소스를 직접 변경하지 않았다.

## 판정

**R7–R9와 R2의 소프트웨어 계약 수정은 타당하다. 최종 확인한 범위에서 미해소 blocker를 발견하지 않았다.** 리뷰 중 BoundaryState 자체의 배열 변경을 탐지하지 못하는 잔여 문제 1건을 확인해 구현자에게 전달했고, 그 보완도 재현으로 확인했다. 물리 모델의 정확성, 특정 실제 처방의 수렴 또는 Zemax와의 수치 동등성을 승인하는 판정은 아니다.

## 집중 검토 결과

| 항목 | 확인한 구현과 판단 |
|---|---|
| R7 설정 전달 | `api.py:121,143`에서 spot 및 common refocus에 동일 settings 전달. 고수준 결과에 기록된 tolerance와 실제 저수준 계산의 tolerance가 일치한다. `direction_tol=1e-15` 사례의 직접 호출/고수준 호출 결과 동등성 통과. |
| R8 배열·광량 | `imaging.py:20–47,70–87`의 bytes-backed 복사로 원본 alias 및 일반적인 WRITEABLE 재활성화를 막는다. `synthesize_image`는 사용 시 다시 grid·mass를 검증하고 비유한 산술/수집 전력 초과를 거부한다. 허용한 부동소수점 수준 residual은 결과에 공개한다. |
| R9 index·격리 | `thermal.py:33–48,80–86`에서 Python/NumPy integer scalar만 수용하고 bool·float·비scalar를 거부한다. 요청 생성 및 실행 경계에서 모두 검사해 고정된 malformed index 문제를 잡는다. batch 실행 중 외부에서 훼손한 group은 해당 열 상태만 실패하고 nominal 결과는 유지한다. |
| R2 추적 provenance | `provenance.py:29–101`은 추적 시점 처방/STOP/공액/실제 n/파장/온도/설정과 ray 내용에 각각 identity를 둔다. 소비 API에서 새 identity를 만들며 오래된 광선을 승인하는 문제가 제거됐다. |
| R2 선택적 데이터 | identity는 현재 단색 추적에 실제 필요한 굴절률을 읽는다. 사용하지 않는 RGB·CTE·흡수를 추가하지 않아도 nominal G 계산이 된다. 흡수 데이터만 바뀌면 기하 광선 재사용을 허용하고 physical source 경로는 광량을 다시 계산한다. |
| R2 IMA 재사용 | image_z는 trace identity에서 제외하고 physical boundary를 고정해 사용한다. IMA 이동만 한 비교 결과가 일치하며 STOP/면/물성/field 변경과 구별된다. |
| R1 및 R3–R6 API | 요청 각주파수의 underresolved 표시, 실제 clipping/수치 실패 분리, LSA 기준면, zero-weight 표본과 quadrature의 분리에 대한 회귀를 실행했다. 새로운 API 결과 구조와 오류 종류의 의미를 숨기는 회귀를 발견하지 않았다. |

## 리뷰에서 발견해 수정한 잔여 문제

**P2, 해소됨 — BoundaryState 배열·scalar의 캐시 내용 변경.** 최초 검토한 `wave.py`는 bundle identity와 source amplitude만 비교하고, cached state의 `opl_mm`, `uv_mm`, reference distance 등의 내용을 비교하지 않았다. 호출자가 cache 배열을 in-place 변경하면 다른 파면을 같은 처방의 결과로 승인할 수 있었다. 이는 정상 호출만으로 저절로 발생하는 결함이 아니라 외부 배열 변경/잘못된 cache 조립이 필요한 조건이다.

독립 재현은 synthetic triplet, 4 rings/16 arms, wave grid 32, padding 4, 60 lp/mm를 사용했다. 정상 MTF는 `0.5060573390158289`였다. `opl_mm += 2e-4 * r²/rmax²`로 cache OPL만 변경하자 `0.5195487885622005`가 되었으며 두 결과 모두 `conditional`, `verified_engine_trace`로 승인됐다. 이때 입력 bundle과 처방은 바꾸지 않았다.

구현자는 `wave.py:388–429,595`에서 전체 snapshot의 방어적 불변 복사와 content seal 검증을 추가했다. 생산 단계만 새 seal을 만들고 `dataclasses.replace`는 기존 seal을 유지하므로 변경된 cache를 재승인하지 않는다. 이후 독립 재검증에서 다음을 확인했다.

- `uv_mm`, `opl_mm`, basis, mask, triangle 등의 WRITEABLE 재활성화는 거부된다.
- `dataclasses.replace`로 OPL, output index, reference distance를 바꾸면 `InputError`다.
- 수정하지 않은 deepcopy/pickle 왕복은 정상 동작한다. 복사본 배열이 쓰기 가능해지는 경우라도 그 내용을 바꾸면 소비 시 content seal 검사가 거부한다.
- 정상 baseline MTF60은 `0.5060573390158289`로 유지됐다.

해시는 accidental mutation/stale cache를 찾기 위한 장치다. 호출자가 내부 private 함수를 호출해 seal까지 다시 만들 수 있는 Python 프로세스에서 악의적인 조작의 보안 경계를 제공한다는 주장은 하지 않는다.

## 실제 독립 실행

모든 실행에서 `PYTHONPATH=src`를 지정했고 실제 import 경로가 작업 루트의 `src/sequential_optics/__init__.py`임을 확인했다. 과거 audit 재현 스크립트를 수정 정답 테스트로 세지 않았다.

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_regression_*.py' -q
```

최종 소프트웨어 재검토 시점 결과: **51 tests, 11.955 s, OK**. 먼저 실행한 R7–R9 12개, R2 최종 12개의 부분 실행도 통과했으며 위 51개와 중복되므로 추가 합산하지 않는다. 전체 원래 테스트까지 포함한 통합 실행은 주 작업자의 최종 검증 기록을 따른다.

회귀 외 독립 probe에서는 새 BoundaryState의 위 변경/복사/왕복 경우를 검사했다. 정상 wave 결과 전체를 `json_value` 후 `json.dumps(..., allow_nan=False)`로 직렬화하는 것도 성공했다(검토 예제 2,678,517 characters). 내부 복소수·비유한 값은 명시적 encoding으로 기록되어 비표준 NaN JSON에 의존하지 않는다. 이 JSON은 결과 보관 형식이며 임의 JSON을 검증된 engine bundle로 자동 복원하는 loader가 있다는 의미는 아니다.

## 적용 범위

R1의 phase-span gate는 부족한 표본에서 그럴듯한 MTF를 내보내지 않도록 하는 필요조건 검사다. 통과만으로 적분 오차가 작다고 보증하지 않는다. R2의 `verified_engine_trace`도 해당 trace의 출처/입력 일치 여부이지 물리적 정확성 인증이 아니다. 공급한 analytic boundary의 명시적 unverified 경로와 일반 engine trace 경로는 계속 구별한다. 기존 scalar 회절·local detector·sampled topology의 한계는 별도 물리/렌즈 리뷰를 따른다.

## 검토 시점 소스 SHA-256

아래 hash는 보고서 작성 순간의 소스다. 이후 수정 시 최종 배포 manifest를 기준으로 다시 연결한다.

| 파일 (`src/sequential_optics/`) | SHA-256 |
|---|---|
| `api.py` | `91efa03580a6b5415e1fa5b476cb51fc004c92860eca77c561285302bd3758ab` |
| `imaging.py` | `f533c66059f066655d456e1fe8089b6fc440fb87fc8c44cee9455fbf494d5c9f` |
| `thermal.py` | `a173e504be65ebee4c1f3f90b49c0660853943a0683a87a7e21ff7c291f637b0` |
| `provenance.py` | `1326c89602233c0ad68d06a962663d16573f545f9a391a0007986302f381b342` |
| `wave.py` | `79727b3d941b4a7ba771d2baa16bbe649c061589c5a58d938657220681b53961` |
| `models.py` | `644917452d9b6002fab1281e1225c2344676d2a736ada9a2a64ffbf24662fd8b` |
| `metrics.py` | `1d9dde79413f65a6c10700e9773318a37bcd2ce7f9a6ef3b7a4be21f22f33659` |
| `source_wave.py` | `5b691fc38c8b2602f51fa21659049b70f543892b1e31f96480cfdb6fad9433b4` |
