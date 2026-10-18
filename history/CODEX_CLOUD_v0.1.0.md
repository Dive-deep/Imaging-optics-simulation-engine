# Codex Cloud로 이어가기

2026-10-18 공식 문서를 확인하여 작성했다. 이 묶음은 업로드할 수 있는 저장소 내용이며, GitHub 업로드나 Codex Cloud 환경 생성은 아직 수행하지 않았다.

## 절차

1. 이 폴더의 **압축을 푼 내용**을 사용할 GitHub 저장소에 넣고 commit/push한다. `pyproject.toml`, `src/`, `HANDOFF.md`, `AGENTS.md`가 저장소 최상위에 있도록 한다. ZIP 파일만 올리면 소스 탐색과 설치가 번거로워진다.
2. 웹 또는 데스크톱의 `Work in > Cloud > Select environment > Create environment`에서 그 저장소를 선택한다. 또는 `Settings > Codex Cloud > Environments`에서 생성한다.
3. 환경 준비 대화에 아래 내용을 전달한다.

   ```text
   이 저장소의 CODEX_CLOUD.md와 HANDOFF.md를 읽어 주세요.
   Python 3.12 환경을 준비하고 bash scripts/cloud_setup.sh를 실행해 주세요.
   python run_tests.py로 기본 테스트를 실행해 설치를 확인해 주세요.
   현재는 환경 준비 단계입니다. audit_round2의 9개 결함은 미수정 상태이므로
   테스트 통과를 물리 정확성 승인으로 표시하지 마세요.
   ```

4. 환경 준비 결과를 검토하고 `Publish`한다. 준비된 환경에서 새 작업을 시작한다. 모바일에서는 준비된 환경을 선택해 작업을 이어갈 수 있으며, 환경 생성은 웹/데스크톱에서 한다.
5. 다음 작업에는 `HANDOFF.md`의 시작 프롬프트를 사용할 수 있다. 결함 수정을 원한다면 해당 수정 프롬프트를 명시적으로 사용한다.

## 파일과 대화 맥락

- 코드·기술 명세·설정·테스트·샘플·이전 리뷰·최신 리뷰를 함께 포함했다.
- `HANDOFF.md`는 현재 상태와 사용자 결정의 인수인계 요약이고, `AGENTS.md`는 저장소 작업 지침이다. 전체 대화의 축어록은 아니다.
- 공식 문서는 New chat에서 기존 ChatGPT 채팅을 Codex 채팅에 추가하는 경로도 안내한다. 해당 기능이 표시되면 원 대화를 추가할 수 있다. 이것만으로 이 작업의 임시 파일, 실행 중 프로세스, 서브 에이전트 상태까지 복제된다고 가정하지 않는다.
- 새 Cloud 작업에는 별도 작업 공간이 생긴다. 같은 Cloud 작업을 다시 열면 그 작업의 저장 상태를 이어간다. 지속적으로 보존할 코드 변경은 Git으로 관리한다.

## 설치와 검증 환경

이 엔진은 Python/NumPy/SciPy로 실행되며 Zemax나 GPU가 필수인 코드는 없다. Python 최소 버전은 기존 pyproject 기준 3.10이며, 실제 검증 환경은 Python 3.12.14 / NumPy 2.3.5 / SciPy 1.17.0이었다.

`requirements-validation.txt`는 기존 검증의 NumPy/SciPy 버전을 재현하는 선택적 constraint다. 이 constraint를 사용하는 Cloud 환경에는 Python 3.12를 권장한다. 일반 설치는 기존 `pip install -e .`로도 가능하다. 스크립트는 서비스나 서버를 띄우지 않는다.

`audit_round2/`는 변경하지 않은 과거 검토 기록이다. 해당 폴더의 재현 스크립트 일부는 의도적으로 내부의 frozen v0.1.0을 import한다. 루트 엔진을 수정한 뒤 이 스크립트만 실행하고 수정 여부를 판정하지 않는다. 루트 `tests/`에 기대 동작을 검증하는 테스트를 추가한다.

## 공식 참고

- [Cloud environments — 환경 생성, GitHub 연결, 설치 및 게시, 저장 상태](https://learn.chatgpt.com/docs/environments/cloud-environments)
- [Projects and chats — 기존 ChatGPT 채팅을 Codex 채팅에 추가](https://learn.chatgpt.com/docs/projects)
- [Iterating Development Workflows with Codex — AGENTS.md와 맥락 파일](https://developers.openai.com/cookbook/examples/codex/iterating-development-workflows-with-codex)
