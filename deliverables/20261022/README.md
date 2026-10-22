# 2026-10-22 후보 생성 결과 전달 파일

최신 작업은 **CT20 후보 23개 전체**다. 이전 후보 24개는 두께 상한을 적용하기 전의 보존 자료다.

- [CT20 후보 23개 ZIP](https://github.com/Dive-deep/Imaging-optics-simulation-engine/raw/refs/heads/main/deliverables/20261022/optics_initial_candidates_ct20_23_20261022.zip)
- [이전 후보 24개 ZIP](https://github.com/Dive-deep/Imaging-optics-simulation-engine/raw/refs/heads/main/deliverables/20261022/optics_initial_candidates_24_20261022.zip)
- [파일 크기·SHA256·다운로드 URL](archive_assets.json)

두 ZIP은 생성 당시의 스냅샷을 바이트 변경 없이 보존한다. ZIP 내부 문서의 당시 GitHub 업로드 상태나 `/workspace/outputs` 경로는 역사적 기록이며, 이후 GitHub 이전 상황은 최신 main 브랜치의 [인계 문서](../../docs/conversation_handoff_20261022_ko.md)와 [새 대화 프롬프트](../../docs/new_chat_prompt_ko.txt)를 따른다.

저장소에는 두 캠페인의 코드·입력·결과·검토·최종 HTML/PDF/전체 단면 그림을 둔다. 재생성 가능한 native `.so`, Python 캐시와 HTML/PDF에 이미 포함된 개별 `reports/plots` PNG는 Git에서 제외한다. 개별 그림은 각 실험의 `render_report.py`로 다시 만들 수 있다. 실행 환경을 바꿔도 이전 감사의 원자료나 hash를 새 결과로 덮어쓰지 않는다.

GitHub 파일 화면에서는 자체 포함 HTML을 실행하지 않으므로 다운로드한 `candidate_report.html`을 브라우저에서 연다. 코드와 후보 데이터를 사용하려면 저장소를 clone하고 `CODEX_START.md`의 설치 절차를 따른다. ZIP을 받기만 하면 ChatGPT 프로젝트의 새 대화가 자동으로 해당 내용을 읽는 것은 아니다. 새 대화에서 저장소 접근을 연결하거나 인계 문서와 필요한 파일을 제공한다.
