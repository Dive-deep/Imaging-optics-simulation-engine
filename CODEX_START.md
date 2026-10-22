# 결상 계산 Engine&DB — Codex에서 시작하기

새 대화에는 [붙여넣기용 프롬프트](docs/new_chat_prompt_ko.txt)를 사용하고 [요구조건·결정·결과·남은 작업 인계](docs/conversation_handoff_20261022_ko.md)를 먼저 읽는다. 이전 대화 원문과 미업로드 파일이 자동 이관되는 것은 아니다.

> 프로젝트 「결상 계산 Engine&DB」를 이어가자. 이 저장소의 AGENTS.md, CODEX_START.md, HANDOFF.md, docs/conversation_handoff_20261022_ko.md와 CT20 METHOD_KO.md를 읽어줘. 최신 23개 후보 전부와 기존 결과·hash를 보존하고, 검증된 기하 결과와 미수렴 파동 결과를 구분해 현재 상태를 복원한 뒤 내 다음 지시에 따라 진행해줘.

## 최신 작업: 2026-10-22 각 렌즈 CT20 제한 후보군

최신 입력은 `experiments/candidates_ct20_20261022/`에 있다. 먼저 [생성·재계산 방법](experiments/candidates_ct20_20261022/METHOD_KO.md)과 [후보 JSON](experiments/candidates_ct20_20261022/candidates.json)을 읽는다. 이전 기하 합격 부모 38개를 대상으로 각 렌즈의 축상 중심 두께 **2 < CT ≤ 20 mm** 조건을 적용해 재계산했고, 합격한 서로 다른 부모의 **23개 후보를 전부 유지**했다. 전체 순위나 단일 winner를 정하지 않았다. 실제 CT는 2.0317579927–19.9999970367 mm이며, 첫 렌즈 전면 clear diameter 50 mm 대안 22개와 EPD 50 mm 대안 1개, 재료 순서 16종이다.

후속 연산은 [동결 처방과 hash](experiments/candidates_ct20_20261022/results/frozen_selection.json)를 기준으로 한다. [최종 집계](experiments/candidates_ct20_20261022/results/final_summary.json), [물리 검증](experiments/candidates_ct20_20261022/results/final_selection_physics_audit.json), [원본 엔진 대조](experiments/candidates_ct20_20261022/results/native_final_candidates/summary.json)를 함께 확인한다. 독립 물리 검사는 23/23개, 새 광선 732,849개 실패 0이며 원본 엔진/native 커널 대조 6,900개 광선도 상태·수치 실패 0이다. 원본 `src/`와 `tests/`는 변경하지 않았다.

주 파동 MTF의 field/파장별 T/S 276쌍 중 200쌍은 `underresolved_quadrature`, 3쌍은 `conditional`, 73쌍은 `not_applicable`이다. 주 경로는 scalar pupil autocorrelation 적분이며 별도 actual FFT 진단과 구별한다. 기하 MTF 대체 표시에는 회절이 없고 파동 실패 상태를 유지한다. 이 후보군은 고NA 파동 정확도가 인증된 설계나 확정 surrogate 학습 label이 아니다. 후보 주변 DB 생성이나 MySQL 저장도 수행하지 않았다.

전체 23개를 [HTML](experiments/candidates_ct20_20261022/reports/candidate_report.html), [PDF](experiments/candidates_ct20_20261022/reports/candidate_layouts_and_mtf.pdf), [전체 단면 모음](experiments/candidates_ct20_20261022/reports/candidate_contact_sheet.png)에서 확인한다. 최종 전달 묶음 경로는 [CT20 후보 23개 ZIP](https://github.com/Dive-deep/Imaging-optics-simulation-engine/raw/refs/heads/main/deliverables/20261022/optics_initial_candidates_ct20_23_20261022.zip)이다. 아래 이전 24개 캠페인과 [기존 ZIP, 약 27 MB](https://github.com/Dive-deep/Imaging-optics-simulation-engine/raw/refs/heads/main/deliverables/20261022/optics_initial_candidates_24_20261022.zip)은 그대로 보존한다.

## 이전 작업: 2026-10-22 고NA 3매 초기 후보 24개

새 클라우드 환경에서 GitHub의 v0.2.1을 복원했다. 이전 대화의 미업로드 속도 최적화 파일은 이 체크아웃에 없었으므로 복원했다고 간주하지 않는다. 원본 코어와 과거 결과는 보존하고, 이번 실험은 `experiments/candidates_20261022/`에 분리했다. 복원 환경에서 원본277개 test가 통과했고 새 기록은 해당 실험의 `results/restored_core_tests.json`에 있다.

후속 작업은 먼저 [후보 생성 과정](experiments/candidates_20261022/METHOD_KO.md)과 [후보 JSON](experiments/candidates_20261022/candidates.json)을 읽는다. 24개 모두 첫 렌즈 opening50mm 대안이며 EPD50mm 설계가 아니다. 최종 처방은`results/frozen_selection.json`으로 고정했고 독립 물리·원본 광선 추적·처방 export 대조를 수행했다. 일부 파동 모델 적용 불가 및 미수렴 곡선은 기하 MTF 진단과 분리되어 있으므로 확정 학습 label로 승격하지 않는다. 그림은`reports/candidate_report.html`과PDF에 있다. 후보 주변 DB를 실제로 생성하거나 MySQL에 저장하지 않았다.

이 저장소는 독립 sequential 광학 엔진 v0.2.1의 후속 개발 기준이다. GitHub 이전은 2026-10-18 사용자 요청으로 승인되었다. 이번 이전은 광학 엔진 소스를 변경하지 않는다.

## 환경 준비

Python 3.12를 권장한다. 저장소 루트에서 다음을 실행한다.

```bash
python -m pip install -e . -c requirements-validation.txt
PYTHONPATH=src:examples python -m unittest discover -s tests -v
```

실행 전 기록은 277개 통과다. 새 환경의 실제 결과를 별도로 확인한다. 핵심 의존성은 NumPy/SciPy이며 pytest는 필요 없다. 설정 중 의존성 설치를 위한 패키지 저장소 접근이 필요하다. 계산 자체는 외부 서비스, API 키, Zemax 또는 GPU를 요구하지 않는다. 서버를 계속 실행할 필요가 없다.

새 로그는 별도 경로에 저장한다. `run_tests.py`는 기존 `results/test_results.json`과 `results/test_run.log`를 덮어쓰므로 역사적 증거를 보존할 때는 위 unittest 명령이나 임시 checkout을 사용한다. 컴파일된 `.so`와 캐시는 Git 배포에서 제외하므로 최신 CT20 커널은 저장소 루트에서 다음 명령으로 빌드한다.

```bash
python experiments/candidates_ct20_20261022/native/build.py
```

CT20 `native/README.md`에는 이전 실험 경로가 남아 있으므로 최신 빌드 경로는 위 명령과 [CT20 방법 문서](experiments/candidates_ct20_20261022/METHOD_KO.md)를 따른다. JSON·manifest의 과거 `/workspace/...` 경로는 당시 환경의 감사 기록이다. 저장된 처방 읽기·재추적은 현재 체크아웃 경로에서 진행할 수 있지만, 부모 원자료를 찾아 전체 이력을 재탐색하려면 과거 경로 prefix를 현재 체크아웃으로 대응시키는 별도 처리가 필요할 수 있다. 이 때문에 기존 증거의 경로 문자열이나 hash를 바꾸지 않는다.

`requirements-validation.txt`는 이전 검증 환경의 선택적 버전 고정이다. 일반 호환 설치는 `python -m pip install -e .`로 가능하다. 그림을 다시 생성할 때만 matplotlib과 한국어 폰트가 추가로 필요하며 엔진 사용에는 필요 없다.

## 먼저 읽을 파일

1. `AGENTS.md`, `HANDOFF.md`: 구현 규칙과 현재 상태.
2. `README.md`, `FIX_C_REPORT.md`, `VALIDATION.md`: API, 변경과 검증 범위.
3. `docs/core_principles_ko.md` 또는 `.tex`: 기본 원리/수식.
4. `docs/calculation_spec_ko.md`: 상세 절차와 출처/미결 쟁점.
5. `docs/field_mtf_usage_ko.md`, `docs/D1_explanation_ko.md`: 축 규약과 미수정 D1.

## 시작 프롬프트 예시

> AGENTS.md, HANDOFF.md, CODEX_START.md, FIX_C_REPORT.md와 docs/conversation_handoff_20261022_ko.md를 읽고 이어갈 준비를 해줘. 새 환경 검증은 과거 결과를 덮어쓰지 않는 별도 기록으로 남기고, 현재 CT20 후보군과 모델의 남은 한계를 요약해줘. 구체적인 후속 작업은 내가 이어서 지정할게.

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
