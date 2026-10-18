# Codex Cloud 인수인계 — Sequential Optics Core

기준일: 2026-10-18. 이 파일은 대화 전체를 복원하는 대신 구현 목적·현재 상태·다음 작업의 근거를 전달한다.

## 현재 요청과 상태

- 사용자의 현재 요청은 **지금까지의 작업을 Codex Cloud로 옮기는 것**이다. 이번 이관 과정에서 결함을 수정하지 않았다.
- 루트의 엔진은 배포된 **v0.1.0**이다. 기존 **183개 테스트는 전부 통과**했지만, 이후 독립 리뷰에서 **미수정 결함 R1–R9**가 재현됐다.
- `audit_round2/`는 2차 리뷰 자료와 당시 고정 입력 소스를 보존한 영역이다. 여기의 재현 코드는 현재 버그가 발생함을 확인한다. 올바른 결과를 검증하는 회귀 테스트와 구별해야 한다.
- **Zemax/Code V와의 수치 상관 검증은 수행하지 않았다.** 테스트 통과나 `conditional` 상태를 실제 처방의 정확성 승인으로 해석하지 않는다.
- 대상 GitHub 저장소는 아직 정해지지 않았으며, GitHub push와 Codex Cloud 환경 생성은 수행하지 않았다.

## 사용자의 목적과 물리 모델

- Zemax를 사용할 수 없는 환경에서, Sequential 모드를 참고하는 재사용 가능한 광학 계산 엔진을 만들어 데이터 생성에 사용한다.
- 기본 대상은 공통 광축의 렌즈 3매. 각 면은 Standard 또는 Even Asphere이며 curvature/radius, conic, A₂·A₄·A₆·A₈·A₁₀·A₁₂를 사용한다.
- 앞면 thickness는 해당 렌즈 중심 두께, 뒷면 thickness는 다음 렌즈까지 거리라는 입력 관례를 지원한다. 내부 모델의 좌표와 adapter 규칙은 README를 따른다.
- active/mechanical diameter와 무효부 높이 입력을 구별한다. 현재 무효부는 analytic sag를 mechanical radius까지 연장한 불투명 영역 모델이며 임의 단차 형상 모델이 아니다.
- 사용자가 object·image·STOP 면, field, EPD 또는 F-number를 정할 수 있어야 한다. 교점·굴절은 **real ray**로 계산하며 EFL/BFL 등 일차 특성은 근축 정의를 사용한다.
- 기본 파장은 설정 파일에 둔 **G=550 nm**이고 변경 가능하다. RGB·열·투과 데이터가 없더라도 해당 데이터가 불필요한 분석은 정상 실행해야 한다. 필요한 물성이 없으면 그 분석만 명확히 실패시킨다.
- 재료는 사용자 제공 모델을 사용한다. 현재 예제의 crown/flint 및 열 물성은 **합성 검증용**이며 실제 상용 재료 데이터라고 주장하지 않는다.
- 광원은 광축 중심의 12.8×3.2 mm, 320×80 정사각 픽셀, 40 μm pitch, fill factor 64%, 각 픽셀 Lambertian이다. RI에는 연속 균일 광원 근사를 사용하고 광량·fill factor를 중복 적용하지 않는다.
- 일반 결상 RI와 무한대 출력의 광도/상대 분포를 구별한다. 보통 FFT MTF는 위 픽셀 광원과 분리한 광학계 성능 지표이며, 기본 G에서 image 공간 lp/mm(요청 범위 0–120)와 diffraction reference를 계산한다. Afocal angular MTF의 cycles/rad를 혼동하지 않는다.
- 균일 주변 온도 변화로 인한 굴절률·자유 열팽창과 focus/MTF/spot 변화를 본다. 발열, 비균일 온도, 공간적 GRIN, 응력·복굴절, mount/주변 구조물 변형, 제조 공차는 현재 범위에서 제외한다.
- 모듈 재사용, 설정 분리, 필요 데이터만 요구하는 API가 핵심 요구다. 수치 실패를 물리적 차단으로 숨기거나 검증되지 않은 근사를 확정값으로 내보내면 안 된다.

## 먼저 읽을 파일

1. `README.md`, `VALIDATION.md`: 현재 API, 실행법, 단위, 모델 제한.
2. `docs/calculation_spec_ko.md`: 수식·참고 문헌·미해결 논점이 포함된 계산 명세. 실제 구현과 차이는 코드 및 후속 리뷰로 확인한다.
3. `audit_round2/CODE_REVIEW.md`: 9개 결함의 근거·코드 위치·우선순위·재검증 기준.
4. `audit_round2/physics/review.md`, `audit_round2/lens_design/review.md`, `audit_round2/software/review.md`: 역할별 독립 검토.
5. `reviews/`: 첫 구현 단계의 샘플 선정 토론·독립 리뷰와 검증 결과. 2차 결함 발견 이전 기록이므로 최신 판정은 `audit_round2/`와 함께 읽는다.

## 미수정 결함 목록

| ID | 우선순위 | 확인된 문제 | 수정 방향 |
|---|---|---|---|
| R1 | P1 | Angular MTF의 요청 각주파수에 대한 출력 표본 부족/aliasing 미검출 | 실제 각도 좌표에서 적분 위상·수렴 검사, 부족하면 refinement 또는 underresolved |
| R2 | P1 | 변경 전 ray bundle/경계장을 변경된 처방에 재사용 가능 | trace 시점 identity 보존·비교, 허용되는 IMA 이동과 재추적 필요 변경 구분 |
| R3 | P1 | STOP 앞 정상 비네팅이 aim_unresolved로 분류되어 분석 전체 실패 | aiming 잔차와 실제 차단 판정 분리; 수치 실패를 무조건 차단으로 바꾸지 않기 |
| R4 | P1 | 보조 bounding cylinder 판별식의 반올림으로 정상 STOP rim 교점 누락 | 크기별 반올림 한계·안정적인 교점 처리; 모든 음수 판별식의 강제 0 처리 금지 |
| R5 | P2 | LSA 기준면 변경이 q=0 실제 축 교점을 덮어씀 | 물리적 paraxial limit와 보고 기준면 분리 |
| R6 | P2 | 최대 spot 반경/plot에서 0 적분 가중치의 rim·chief 표본 제외 | RMS용 적분 mask와 기하 표본 mask 분리; rim에 임의 양의 가중치 부여 금지 |
| R7 | P2 | evaluate의 spot/common refocus 호출에 사용자 settings 전달 누락 | 설정 전파와 고수준/저수준 동등성 검증 |
| R8 | P2 | PixelGrid/SourceKernel 내부 배열 변경으로 전력·면적 불변식 우회 | 방어적 copy와 배열 잠금, 계산 경계 불변식 검증 |
| R9 | P3 | thermal group의 float/bool index 검증 누락으로 batch 중단 | 명시적인 정수형 정책과 입력 단계 검증 |

## 실행 및 재현

Python ≥3.10, NumPy ≥1.24, SciPy ≥1.10. 다음 명령은 저장소 루트에서 실행한다.

```bash
python -m pip install -e .
python run_tests.py
python examples/run_demo.py
```

- R1/R2: `audit_round2/physics/repro_angular_otf_alias.py`, `repro_stale_bundle.py`.
- R3–R6: `audit_round2/lens_design/reproduce_lens_findings.py`.
- R7–R9: `audit_round2/software/reproduce_software_findings.py`.
- 재현 스크립트의 import 경로는 **리뷰 당시 고정 사본**을 가리킬 수 있다. 수정 검증에서는 루트 `src/sequential_optics`를 실제로 읽는지 먼저 확인한다. 고정 사본 결과와 수정본 결과를 섞지 않는다.
- 과거 결과는 `results/`, `audit_round2/reproductions/` 및 역할별 JSON에 있다. 실행 환경·처방·설정·주파수·표본 수가 다른 숫자를 직접 비교하지 않는다.

## 수정 작업을 요청받았을 때의 권장 순서

1. P1 중 R2 provenance와 R1 angular MTF 검증을 먼저 처리하고, 이어 R3/R4의 비네팅·교점 분류를 수정한다.
2. R5–R9를 수정하고 각 재현 사례를 **기대되는 올바른 동작**을 검사하는 회귀 테스트로 전환한다.
3. 기존 183개와 새 회귀 테스트를 함께 실행하고, 대표 처방의 ray density·FFT grid·padding·요청 주파수를 각각 바꾸어 수렴을 확인한다.
4. 광학 물리학자·렌즈 설계 엔지니어·소프트웨어 개발자 역할의 독립 에이전트로 수정 코드와 계산 절차를 다시 리뷰한다. 샘플의 타당성도 검토한다.
5. 수정·검증한 범위와 남은 논점을 문서화한다. Scalar 회절, 비축 detector 투영, 실제 재료 데이터 및 상용툴 상관 검증의 제한을 별도로 유지한다.

## Codex 시작 프롬프트

아래 프롬프트는 **9개 결함 수정까지 요청하고 싶을 때** 복사해 사용한다. 이관만 한다면 여기까지 실행할 필요는 없다.

> HANDOFF.md, README.md, VALIDATION.md와 audit_round2/CODE_REVIEW.md를 읽고 현재 상태를 확인해줘. v0.1.0의 기존 183개 테스트가 통과해도 R1–R9는 아직 미수정이다. 이번에는 이 9개를 P1부터 수정해줘. 먼저 각 재현의 물리·수치적 근거를 검토하고, 고정 사본이 아닌 작업 중 소스에 대해 기대 동작의 회귀 테스트를 만들어라. 선택적 열/분광 데이터, 설정 분리, real-ray 계산과 사용자 범위를 유지해라. 수정 후 전체 테스트와 필요한 수렴 검증을 실행하고, 광학 물리학자·렌즈 설계 엔지니어·소프트웨어 개발자 역할의 독립 에이전트 리뷰를 거쳐 발견 문제를 해결해라. 수정되지 않은 사항과 모델 근사를 분명히 기록하고 Zemax와 수치 일치했다고 주장하지 마라.
