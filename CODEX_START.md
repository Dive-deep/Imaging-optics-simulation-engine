# Codex에서 시작하기

이 저장소는 독립 sequential 광학 엔진 v0.2.1의 후속 개발 기준이다. GitHub 이전은 2026-10-18 사용자 요청으로 승인되었다. 이번 이전은 광학 엔진 소스를 변경하지 않는다.

## 환경 준비

Python 3.12를 권장한다. 저장소 루트에서 다음을 실행한다.

```bash
python -m pip install -e . -c requirements-validation.txt
python run_tests.py
```

실행 전 기록은 277개 통과다. 새 환경의 실제 결과를 별도로 확인한다. 핵심 의존성은 NumPy/SciPy이며 pytest는 필요 없다. 설정 중 의존성 설치를 위한 패키지 저장소 접근이 필요하다. 계산 자체는 외부 서비스, API 키, Zemax 또는 GPU를 요구하지 않는다. 서버를 계속 실행할 필요가 없다.

`requirements-validation.txt`는 이전 검증 환경의 선택적 버전 고정이다. 일반 호환 설치는 `python -m pip install -e .`로 가능하다. 그림을 다시 생성할 때만 matplotlib과 한국어 폰트가 추가로 필요하며 엔진 사용에는 필요 없다.

## 먼저 읽을 파일

1. `AGENTS.md`, `HANDOFF.md`: 구현 규칙과 현재 상태.
2. `README.md`, `FIX_C_REPORT.md`, `VALIDATION.md`: API, 변경과 검증 범위.
3. `docs/core_principles_ko.md` 또는 `.tex`: 기본 원리/수식.
4. `docs/calculation_spec_ko.md`: 상세 절차와 출처/미결 쟁점.
5. `docs/field_mtf_usage_ko.md`, `docs/D1_explanation_ko.md`: 축 규약과 미수정 D1.

## 시작 프롬프트 예시

> AGENTS.md, HANDOFF.md, CODEX_START.md, FIX_C_REPORT.md를 읽고 이 광학 엔진의 개발을 이어갈 준비를 해줘. 우선 설치 후 python run_tests.py로 기준 상태를 확인하고, 계산 모델의 범위와 남은 개선사항을 요약해줘. 이번 단계에서는 광학 계산식이나 D1 정의를 변경하지 마.

GitHub 파일은 대화 전체를 자동 전달하는 기능이 아니다. 지속할 요구조건·결정·제한은 위 문서에 정리했고, 새 작업의 구체적인 목표는 프롬프트로 지정한다.

## 미해결 사항과 계약

- D1은 성분별 f_theta chart의 정의 문제로 설명만 추가했고 코드는 변경하지 않았다. 방사 방향 F-θ 정의 변경은 별도 작업이다.
- Field마다 T/S MTF, complex OTF, diffraction reference를 반환한다. 고정 chief-normal X/Y 규약이며 field-radial T/S로 자동 회전하지 않는다. 임의 2D field의 Zemax 수치 동치는 미검증이다.
- Real-ray 교점·Snell·OPL과 근축 정의 EFL/BFL, scalar FFT를 구별한다.
- 기본 G는 설정으로 관리한다. 사용하지 않는 RGB/열/흡수 자료는 요구하지 않는다.
- Lens-only MTF와 실제 Lambertian source RI/영상은 별도 모델이다.
- 균일 열 변화만 범위다. 응력·비균일 온도·GRIN·기구 FEA·제조 공차는 포함하지 않는다.
- 수치 실패를 물리적 clipping으로 숨기지 않는다. underresolved/conditional은 물리 정확도 인증이 아니다.
- 실질적인 엔진 변경 후 물리학자·렌즈 설계자·소프트웨어 역할의 별도 에이전트 리뷰를 수행한다.

## 예제와 재생성 결과

```bash
python examples/run_demo.py
python examples/run_source_wave.py
python examples/run_field_mtf.py
python -m sequential_optics examples/triplet_G_only.json --settings examples/settings.json --output results/fix_c_G_only_cli.json
```

현재 소스/테스트/문서/작은 검증 기록은 Git에 포함했다. 큰 실행 결과 23개(약49 MB)는 테스트 입력이 아니므로 Git에서 생략했고, 원본 `sequential_optics_core_v0.2.1.zip`에 그대로 보존되어 있다. 목록과 원본 해시는 `MIGRATION_MANIFEST.json`의 `omitted_generated_files`에 있다. 과거 audit의 큰 결과는 해당 ZIP에서 복구해야 당시의 바이트가 유지된다. 새로 실행한 결과를 옛 증거로 덮어쓰지 않는다.

`history/SOURCE_MANIFEST_v0.2.1_full_zip.json`은 전체 ZIP의 기록이며 Git 체크아웃의 전체 목록이 아니다. 루트 `SOURCE_MANIFEST.json`은 이번 이전 파일의 기록이다. `audit_round2/input/optical_engine`는 과거 고정 소스이고 현재 작업 대상은 루트 `src/`, `tests/`다.

## Cloud 연결

웹 또는 데스크톱의 Codex에서 Cloud 환경을 만들고 이 GitHub 저장소를 선택한다. 위 설치/테스트를 준비 단계에서 실행하게 한 뒤 환경을 Publish하고 새 작업을 시작한다. 모바일에서는 준비된 환경을 선택하여 작업을 이어갈 수 있다.

공식 절차: https://learn.chatgpt.com/docs/environments/cloud-environments (2026-10-18 확인). 이 문서는 GitHub 자료 준비를 완료하는 것이며 사용자의 Codex 환경을 대신 생성하거나 Publish했다고 주장하지 않는다.
