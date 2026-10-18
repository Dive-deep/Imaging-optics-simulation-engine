# 후속 리뷰 9개 결함 수정 보고서 — v0.2.0

2026-10-18. 대상은 원래 v0.1.0에 대한 `audit_round2/CODE_REVIEW.md`의 R1–R9다. 마이그레이션은 보류했고 엔진과 문서만 수정했다.

## 결과

**9개 결함을 모두 수정했으며 기존 183개 + 신규 회귀 52개 = 235개 테스트가 통과했다.** 실패·오류·skip은 0이며 전체 실행 시간은 약 58.14초였다. 두 예제 스크립트(`run_demo.py`, `run_source_wave.py`)와 G-only CLI도 종료 코드 0으로 실행 완료했다. 실행한 소스는 루트 `src/sequential_optics/`다. 원래 버그의 존재를 확인하는 audit 스크립트를 정상 동작 테스트로 합산하지 않았다.

수정 구현은 역할별로 분담했고, 이후 구현자와 다른 광학 물리·렌즈 설계·소프트웨어 에이전트 3명이 코드와 실제 재현을 다시 검토했다. 이는 역할을 부여한 AI 리뷰이며 사람 전문가 인증이나 상용툴 상관 검증은 아니다. 최종 검토 범위에서 미해소 차단 결함은 발견되지 않았다.

## 결함별 수정과 대조

| ID | 변경 | 검증 근거 |
|---|---|---|
| R1 | 각 angular cell의 실제 각도 범위와 요청 주파수의 적분 위상을 검사. 저수준은 부족한 해상도를 거부하며 고수준은 해당 주파수 NaN와 `underresolved`를 반환. 실제 pupil과 reference, physical source 경로에 공통 적용. | Gaussian alias 사례에서 기존 거짓 고MTF를 거부. Padding 16/24에서 같은 주파수가 수렴. High-angle cell 대조로 direction-cosine cutoff와 구별. DC 보존. |
| R2 | 추적 시점의 처방·STOP·OBJ·field·파장·온도·n·ray 설정과 전체 ray 내용 기록. 소비 API에서 비교. BoundaryState는 불변 snapshot 및 내용 seal로 보호. | 마지막 곡률 0.5%, STOP 반경/index, n(T), field 등 변경 거부. IMA-only/FFT sampling 재사용 허용. 수동 analytic bundle은 명시 opt-in과 unverified 표시. 흡수 변경 시 전력 재계산. |
| R3 | STOP 목표의 실제 sag로 정한 3D 점을 향해 aperture와 분리한 가상 연장 잔차로 aiming. 수렴 후 실제 STOP 첫 교점을 확인하고 별도 물리 추적에서 clipping 적용. | 두 평면 aperture의 독립 overlap mask 134/192, RMS와 quadrature fraction 일치. 구면 굴절의 독립 Snell rim mapping. 가상 conic-domain 실패 및 뒤쪽 root를 정상 처리로 숨기지 않음. |
| R4 | Cylinder 판별식의 스케일 기반 roundoff 한계와 안정적 이차식 해. 평면은 직접 교점. Bounding과 aperture에 같은 길이 tolerance 적용. | 음의 미소 판별식 및 curved STOP rim 유지. 3개 길이 스케일의 실제 miss와 실제 optical tangency는 계속 거부. |
| R5 | 실제 근축 축 교점과 보고 기준면 분리. q=0은 실제 물리적 극한 유지. | 실제 focus 91.92430937348901 mm, 기준면 +10 mm에서 q=0 LSA=−10 mm. 인접 작은 pupil ray의 극한과 연속. |
| R6 | 적분 표본과 기하 표본 mask 분리. RMS/centroid/covariance에는 기존 양의 가중치, 좌표와 표본 최대 반경에는 성공한 rim/chief도 포함. | 기본 triplet 표본 최대 반경 0.00903228370550245 mm, RMS 0.004556286083574049 mm. Rim 포함 여부와 weighted 통계 동일. |
| R7 | 상위 spot/common focus에 사용자 settings 전달. | direction tolerance 경계의 실제 처방에서 고수준/저수준 spot·공통 focus 일치. |
| R8 | PixelGrid/SourceKernel 배열을 bytes-backed 불변 복사로 생성하고 계산 경계에서 불변식 재검사. 전력 위반 거부, 미세 반올림 residual 공개. | 원본 alias, in-place 변경, write flag 재활성화 차단. 음수/비유한 면적·전력 거부. 전력 ledger 대조. |
| R9 | LensGroup index는 Python/NumPy 정수 scalar만 허용하고 bool/float 거부. 생성 및 실행 경계 검증. | 잘못된 입력은 InputError. 런타임 손상 그룹은 열 상태에서 실패하고 nominal 상태 보존. 정상 NumPy 정수 정규화. |

R1의 guard는 적분 수렴의 증명이 아니다. 필요 해상도를 확보하지 못한 데이터를 채택하지 않기 위한 판정이며 설정에서 phase-span 한계와 허용 미해상 전력 비율을 조절한다. R6의 최대 반경 역시 추적 표본의 최대값이지 연속 pupil의 전역 극값 보증은 아니다.

## 독립 리뷰에서 추가 확인하고 반영한 사항

1. **BoundaryState의 내용 변경:** bundle의 출처만 확인하면 cache 내부 OPL 변경을 놓친다는 소프트웨어 리뷰 지적을 반영했다. 배열의 불변성, 전체 snapshot seal, dataclasses.replace 및 deepcopy/pickle 후 변경 탐지를 확인했다.
2. **곡면 STOP rim의 tolerance 불일치:** aiming이 허용한 미소 위치 오차를 physical cylinder가 잘라 유효 교점을 놓치는 렌즈 리뷰 반례를 반영했다. 독립 plane Snell로 목표 광선을 구해 수정 전 실패/수정 후 통과를 확인했다.

세 독립 보고서는 `reviews/fix9_physics_review.md`, `reviews/fix9_lens_design_review.md`, `reviews/fix9_software_review.md`다. 각 보고서의 부분 테스트 수는 해당 검토 시점의 수이며 서로 중복된다. 최종 전체 수는 235개다.

## 새 원리 문서

`docs/core_principles_ko.md`는 예외 목록이나 버그 사례 대신 입력→real ray→OPL/OPD→scalar FFT→광량→영상→균일 열 상태의 기본 원리를 설명한다. 핵심 번호식 44개와 Sellmeier 식, 참고문헌 15개를 담았다. 같은 내용의 `.tex`를 제공한다.

Markdown와 LaTeX를 각각 Pandoc AST로 읽어 **총 231개 수식 표현이 순서·내용을 유지**하는지 검사했다. 별도 물리·렌즈 리뷰에서 핵심 식과 실제 구현의 정의를 대조했다. LaTeX 파일은 XeLaTeX 및 koTeX가 있는 환경에서 컴파일하도록 구성했다. 이 실행 환경에는 koTeX가 없어 완전한 문서 조판 검증은 수행하지 않았으며, PDF는 제공하지 않는다. Markdown/LaTeX의 수식 원문과 의미 보존 검사는 통과했다.

## 유지되는 모델 범위

기본 G에서 선택적 물성 없이 실행하는 계약과 사용자 설정 분리는 유지했다. Real-ray 추적, 1차 EFL/BFL, scalar 회절 모델의 역할은 서로 다르다. 원래의 scalar·단일 경로·sampled pupil·국소 chief-normal FFT 및 선택적 비축 detector 투영 제한은 바꾸지 않았다. 실제 상용 재료 물성·Zemax 수치·측정 데이터와의 상관 검증은 아직 없다.

검증 원본은 `results/test_results.json`, 사람이 읽는 로그는 `results/test_run.log`다. `audit_round2/`의 역사적 파일은 원래 manifest와 비교하여 그대로 보존한다.
