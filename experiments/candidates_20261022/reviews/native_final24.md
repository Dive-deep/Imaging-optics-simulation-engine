# 최종24개 처방의 원본 엔진 독립 대조

최종 `results/frozen_selection.json`의 C01–C24를 **재선별하거나 제외하지 않고**
모두 원본 v0.2.1 Python 엔진으로 다시 계산했다. 일부 처방이 앞선38개 감사
이후 비구면·초점·clear aperture 조정을 받았으므로 이전 결과를 그대로
적용하지 않았다. 이번24개 검증은 `results/native_final24/`에 별도로 보존한다.

실행 전에 전체 선정 파일의 SHA256과24개 각각의 `design_sha256`을 확인했다.
검사 시작·종료 사이 선정 파일, 모델, 원본 geometry/models/materials,
네이티브 C/Python/binary의 변경은 없었다. 입력/소스 hash와 감사 script hash는
`results/native_final24/plan.json`에 있다.

각 처방의 RGB486/588/656nm, field0/0.5/0.8/1F에서 chief1개·rim8개·내부16개를
새로 계산했다. 원래16rings×48arms 표본에서 고정한25개 pupil 점이며
처방당300개, 전체 **7,200개 광선**이다. 각 광선의6개 실제 굴절 경계면
위치·굴절 방향·누적OPL 및 고정 상면 좌표를 비교했다.

| 항목 | 최종24 결과 |
|---|---:|
| 새로 계산한 원본/네이티브 광선 | 7,200쌍 |
| 원본 방사 구간 geometry 검사 | 24/24 통과 |
| 상태·수치 대조 실패 | 0 |
| 최대 경계면 위치 차이 | 2.92779e−11mm |
| 최대 방향 성분 차이 | 1.19693e−12 |
| 최대 OPL 차이 | 7.01306e−11mm |
| 최대 상면 좌표 차이 | 3.65610e−11mm |
| 선정 입력·엔진/모델 소스 변경 | 0 |
| 단일 CPU 감사 wall time | 56.9597초 |

위치·OPL의 비교 한계는 절대1e−8mm, 방향은 절대2e−10이다. 경계면 좌표
최대는 C04/486nm/1F의 sagittal rim에서, 방향·상면 좌표 최대는
C04/656nm/0.8F의 tangential rim에서 관측됐다. OPL 최대는
C15/486nm/1F의 tangential rim에서 관측됐다. 모두 각 case JSON에 기록했다.

고차 비구면의 첫 전방 교점 선택과 실제 Snell 추적에서 원본과 불일치한
표본은 없었다. 원본 `validate_geometry`의 방사 구간 검사도24개 모두
성공했으므로 유한 점 표본의 면 간격 검사만으로 통과시킨 결과와 구분된다.
검사 영역은 인접면의 공통 clear radius이며, 외곽 기계적 렌즈 형상·chamfer·
가공 공차·마운트 간섭은 여기서 인증하지 않는다. 이 수치 구간 알고리즘을
임의 계수에 대한 형식적 증명이라고 표현하지 않는다.

검사는 원본과 네이티브의 선택된 광선 결과를 비교한 것이다. 연속 pupil의
전역 극값, 모든817개 광선에 대한 원본 재계산, FFT MTF 수렴, 고 NA 벡터 회절,
전역 최적 설계 또는 Zemax 수치 동치를 증명하지 않는다. 파동 출력의
`conditional`/`underresolved` 판정은 별도 보고를 그대로 따라야 한다.

근거 파일:

- `native/audit_final24.py`: 이전 감사에서 입력만 최종 freeze로 바꾼 별도 script.
- `results/native_final24/plan.json`:24개 선정·처방 hash·표본·소스 hash.
- `results/native_final24/summary.json`: 전체 판정과 최대 오차 위치.
- `results/native_final24/C01.json` … `C24.json`: 처방별 독립 결과.

이번 감사 중 원본 코어·모델·네이티브 커널을 수정하지 않았다. 수치 작업은
최종24개 대조를 마친 뒤 종료했다.
