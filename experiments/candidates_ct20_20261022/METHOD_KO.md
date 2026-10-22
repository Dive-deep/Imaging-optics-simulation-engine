# 3매 광학계 초기 후보군 — 각 렌즈 중심 두께 20 mm 제한 반영

작성일: 2026-10-22 UTC. 사용자의 추가 조건을 반영하여 **각 렌즈의 광축 방향 중심 두께를 2 mm 초과, 20 mm 이하**로 제한했다. 조건을 통과한 서로 다른 출발점의 후보들을 모두 제공한다. 단일 우승 설계를 선정하는 작업이 아니다. 기존 `../candidates_20261022/` 처방과 보고서는 변경하지 않았다.

`candidates.json`은 전체 후보의 입력·출력·생성 이력을 담은 파일이다. `reports/candidate_report.html`은 전체 후보 단면을 먼저 보여주고 각 후보의 상세 그림을 제공한다. PDF에는 후보마다 한 페이지가 있다. 이 문서는 재생성 과정과 결과의 해석 범위를 설명한다. DB 서버나 MySQL 테이블은 생성하지 않았다.

## 1. 설계 조건과 판정 기준

| 항목 | 이번 처방의 정의 |
|---|---|
| 물체·좌표 | 무한 물체, 중심 회전대칭, +Z 진행, 첫 면 꼭짓점 Z=0 |
| 화각 | ±12°, 전체 24°; 성능 field는 0/6/9.6/12° = 0/0.5/0.8/1F |
| 상면·상고 | 모든 RGB와 field가 공유하는 평면. 588 nm의 12° chief ray가 +6.4 mm, 허용 잔차 0.001 mm |
| 파장 | 486/588/656 nm, 입사 RGB 파워 동등 가중 |
| 매수·재료 | 공기 간격을 둔 3매. PMMA/PC 중 두 장 + FDS18W/TAFD45 중 한 장 |
| 렌즈 두께 | 각 매의 축상 중심 두께 **2 < CT ≤ 20 mm**. 총 렌즈군 길이나 주변 두께의 상한을 뜻하지 않음 |
| 간격 | 두 축상 렌즈 간 air gap ≥1 mm |
| BFL | 마지막 면 꼭짓점에서 공통 상면까지 ≥3 mm; 별도로 588 nm 근축 BFL ≥3 mm |
| NA | 공기 중 실제 출사 방향으로 계산. **축상 전체 rim의 최소 NA가 세 RGB 각각 ≥0.7** |
| 구경 A | 앞 STOP 지름 50 mm = EPD 50 mm |
| 구경 B | 첫 렌즈 전면 clear diameter 50 mm. 실제 EPD를 별도로 보고 |
| STOP | 첫 렌즈 앞, 위치를 변수로 탐색. 실제 면과 간섭하지 않고 EFL보다 가까운 조건 확인 |
| 플라스틱 | standard 또는 even_asphere, conic과 A4/A6/A8/A10/A12 |
| 유리 | standard 구면만 허용. conic=0, 모든 비구면 계수=0 |
| 제외 | 열해석, 접합면, 다른 면 형식, A14 이상 항 |

NA는 `1/(2 F-number)`로 대체하지 않는다. 주변 field의 cone NA도 기록하지만 전 시야 모든 방향에 NA≥0.7을 강제한 조건은 아니다. 6.4 mm는 기준 파장 chief의 위치이며 RGB 모두의 chief나 전체 spot이 12.8 mm 센서 안에 들어온다는 보장이 아니다. detector crop은 적용하지 않았다.

최종 일반 광선 평가는 STOP 16 rings ×48 arms의 면적 가중 annuli에 chief와 rim을 더한 817개 광선/field/파장을 사용한다. 구경·시야 판정은 0–12°의 13개 field ×3파장 ×817=31,863개 광선/후보다. 검사한 표본 안에 실패가 있으면 평균값이 좋아도 통과시키지 않는다.

## 2. 재료와 수식

| 재료 | nd, d=587.5618 nm | Vd |
|---|---:|---:|
| PMMA | 1.492 | 57.47 |
| PC | 1.586 | 32.00 |
| FDS18W | 1.94594 | 17.98 |
| TAFD45 | 1.95375 | 32.32 |

사용자가 채택한 normal-line 모델을 유지했다. `PgF=0.6438−0.001682 Vd`, `nF−nC=(nd−1)/Vd`를 사용하고 `n(λ)=a+b/λ²+c/λ⁴`의 세 계수를 nd·F−C·g−F 조건으로 결정했다. 이 식의 λ 단위는 μm, 기준선은 d587.5618/F486.1327/C656.2725/g435.8343 nm다. 요청한 486/588/656 nm의 실제 모델값과 계수는 JSON에 있다. 재료 이름을 이유로 사용자 수치를 catalog 값으로 바꾸지 않았다. 측정 분산이나 anomalous partial dispersion까지 복원한 모델은 아니다.

표면 sag는 `z(r)=c r²/[1+sqrt(1−(1+k)c²r²)]+Σ A_m r^m`, m=4/6/8/10/12다. 길이는 mm, c는 mm⁻¹, A_m은 mm^(1−m). 최적화에서는 수치 스케일을 맞추기 위해 `b_m(r/25)^m`으로 쓰며 export 시 `A_m=b_m/25^m`로 변환한다. curvature=0인 평면의 radius는 JSON null이다. 각 면의 thickness는 다음 면까지의 축상 거리이며 마지막 면에서는 공통 상면까지의 거리다.

실제 교점·벡터 Snell·OPL은 별도 실험의 FP64 C 커널로 계산한다. clipping/TIR/교점 실패/뒤쪽 진행 상태를 구분한다. 이번 처방 검증을 FP16에서 수행했다고 표현하지 않는다. 원본 `src/`와 `tests/`는 변경하지 않았다.

## 3. 서로 다른 출발 구조를 만든 과정

이전 탐색의 112개 출발점(수식 기반 96개, 공개 특허 변환 16개)을 그대로 추적 가능한 부모 이력으로 사용했다. 그중 이전 구경·광선 조건을 통과했던 서로 다른 38개 부모의 가장 좋은 저장 checkpoint를 이번 두께 제한으로 다시 탐색했다. `repair_inputs.json`에 원래 경로·SHA256·처방·이전 최대 CT가 있다. 기존 최종 24개 중 CT20 조건을 이미 만족한 것은 두 개였으며, 전체 38개 부모 중에는 세 개였다. 나머지를 두께 숫자만 잘라서 후보로 내보내지 않았다.

수식 기반 출발점은 재료 순서, 유리 앞/중간/뒤 위치, 양/음 파워 배열과 bending을 바꿨다. Rectilinear 기준 목표 EFL은 `6.4/tan(12°)=30.1096327 mm`다. `Σφ≈1/f`, `Σφ/V≈0`을 출발 조건으로 사용하고 reduced-ray `(y,n u)` ABCD 행렬로 실제 두께를 포함한 파워를 보정했다. 근축 일치는 실제 광선 합격을 대신하지 않는다.

특허 기반 출발점은 [US3944337A](https://patents.google.com/patent/US3944337A/en)의 수치 예 네 가지다. 앞 STOP·플라스틱 두 장의 conic 네 면·구면 유리 한 장이라는 구조를 참조했다. 원래 f/1.4–1.57 처방을 길이 조정하고 허용 재료로 대체한 뒤 구경과 실제 광선을 다시 최적화했다. 원문 성능이 그대로 유지된다고 가정하지 않았다. 원문에는 12차 초과 항이 없으며 이번 export도 A12까지만 존재한다.

원래의 수치 전사·재료 대체·구경 확대 과정은 [이전 방법 기록](../candidates_20261022/METHOD_KO.md), [출처 조사](../candidates_20261022/references/RESEARCH_KO.md)에 보존되어 있다. 각 최종 후보의 `provenance.generation_history`에는 이전 과정과 이번 수정이 연결된다.

## 4. CT20 제한으로 다시 계산한 과정

`repair_ct20.py`는 38개 부모를 전부 처리했다. 먼저 현재 최대 중심 두께에서 상한 20 mm까지 한 번에 최대 2.5 mm씩 상한을 내리는 continuation을 수행했다. 각 단계에서 여섯 면 곡률, 플라스틱 conic/A4–A12, 중심 두께, 공기 간격, 상면, STOP을 함께 재최적화했다. 큰 CT를 단순히 20으로 자르면 주변에서 렌즈면이 겹칠 수 있다는 렌즈 설계 검토를 반영한 것이다.

각 단계는 3×8 pupil과 4개 성능 field/RGB, 최대100 nfev, 종료 tolerance 2e−6을 사용했다. 뒤이어 6×16 pupil에서 공통 RGB centroid와 복소 geometric OTF를 포함한 merit로 최대150 nfev, tolerance1e−6의 polish를 수행했다. nfev는 유한차분 Jacobian을 포함한 전체 residual 호출 수와 같지 않다. 로그에 실제 단계별 종료 상태와 시간을 보존했다.

기본 탐색 범위는 플라스틱 curvature ±0.1/mm, 유리 curvature ±0.98/clear-radius, conic −40…2, 정규화 b4…b12 각각 ±30 mm, CT 2.02…20 mm, gap1…20 mm, 상면 BFL3…50 mm, 앞 STOP 거리0.1…28 mm다. 중간 continuation에서만 이전 CT를 연결하기 위해 상한20 초과가 허용되며, 최종 합격 검사에는 절대로 그 중간 상한을 사용하지 않는다. 합격 판정은 실제 저장된 CT≤20을 정확히 검사한다.

상고의 마지막 미세 보정은 모든 field/RGB가 공유하는 상면을 한 번 이동하는 방식이다. CT를 늘릴 수 있는 전체 similarity scaling은 새 CT20 재계산에서 사용하지 않았다. 이전 처방 이력에 있는 scaling은 이전 단계의 기록으로만 남는다.

빠른 merit의 geometric OTF는 광선 상면 좌표의 복소 합이며 회절이 없다. 20 lp/mm 한 점만 높이는 대신 여러 주파수의 응답, T/S, RMS와 RI 편차를 함께 사용한다. 이 단계가 파동 MTF 최적화 또는 전 대역 단조 MTF를 보장하지는 않는다.

## 5. 실패 사례, 보완 실험과 STOP 스캔

초기 재계산 뒤 seed0059/0048/0012에는 축상 NA가 충분하지만 물리 clear 가장자리에서 차단되는 광선이 남았다. `rescue_ct20.py`로 광선 연결을 돕기 위한 가상 trace aperture 확대와 **실제** clear의 footprint penalty를 사용하여 6×16/최대250 nfev 보완을 시도했다. 가상 aperture는 최적화 중에만 사용하며 후보 수락에는 원래 물리 clear와 16×48 재추적을 사용했다. 세 결과에는 각각13/1/40개의 실패 광선이 남아 그대로 통과시키지 않았다.

소프트웨어 검토는 무조건적인 가상 aperture 확대가 일부 conic의 실수 domain을 넘어설 수 있음을 확인했다. 별도 `rescue_ct20_guarded.py`에 domain 제한과 최대 footprint penalty를 구현하고 seed0075의300-ray control로 확인했다. 이 guarded 버전은 이번 최종 후보의 추가 최적화에는 실행하지 않았다. 구현·control 증거와 실제 사용 경로를 구분한다.

대신 첫 렌즈 opening50 대안에서 STOP 반경을 원래 값 주변 ±0.60 mm, 0.02 mm 간격으로 검사했다. 여섯 면 곡률/conic/A, CT/gap과 상면은 고정하고, 첫 전면 clear radius25 mm도 유지했다. 후속 qualification은 사용하지 않는 다른 면 clear만 다시 축소할 수 있다. EPD50 모드에는 이 스캔을 적용하지 않았다. STOP 범위의 기존 하한21.6 mm는 사용자 제약이 아니라 이전 탐색의 box bound였으므로 그보다 작은 값도 실제 NA≥0.7과 광선 조건을 다시 만족하면 허용했다.

각 trial은 먼저 16×48,4field/RGB를 검사하고 합격 가능성이 있는 경우 13field 전체와 두께·NA·구경을 다시 판정했다. seed0059는3개, seed0048은9개 STOP 표본이 최종 기준을 통과했고 seed0012는0개였다. 같은 부모 안에서만 한 checkpoint를 남겼다. 두 부모가 추가되어 최종 서로 다른 후보는23개다. `results/stop_scan/scan_log.json`, 후보 이력의 `ct20_stop_scan`에 원래 반경과 변경량이 있다. 실패한0012를 후보 수에 포함하지 않았다.

## 6. 기하 및 독립 검증

`qualify.py`는 13field/RGB 전체를 실제 clear로 새로 추적한다. 모든 광선이 통과하는 처방만 수락하며, 각 면의 최대 footprint에 max(0.30 mm,2%) 여유를 둬 사용하지 않는 clear를 축소한다. 조정 전후 광선 수·상태·경계 좌표·방향·OPL을 대조하고 통과하던 광선을 제거하는 조정은 거부한다.

각 면의 실수 sag domain, 렌즈면/공기 간격의 radial 표본 간섭, STOP/상면 간섭, 재료와 면 유형, CT의 열린 하한/닫힌 상한, gap, 두 종류의 BFL, 실제 RGB rim NA, 상고를 별도 gate로 검사한다. `validate_ct20_boundary.py`의15개 control은 20 통과, nextafter(20,+∞) 실패, 2 실패, nextafter(2,+∞) 통과 및 중간 최적화 상한으로 최종 검사를 완화할 수 없음을 확인했다.

최종 동결 처방의 독립 물리 검사는 `validate_selection_physics.py`, 원본 Python 엔진과 별도 C 커널의 대조는 `native/audit_final_candidates.py`로 실행한다. 결과는 `results/final_selection_physics_audit.json`, `results/native_final_candidates/`, 역할별 `reviews/`에 기록한다. 검증은 원래 캠페인의24개 결과를 재사용하지 않고 CT20 최종23개 입력으로 다시 수행한다. 4097 radial samples와 여러 pupil 표본의 통과는 연속 영역 전체의 수학적 증명이나 제조 공차 인증을 뜻하지 않는다.

실제 최종 독립 검사는 **23/23 통과, 732,849개 광선 중 실패0개**였다. 원본 엔진의 연속 radial-interval geometry 검사도23/23 통과했고, 원본/native의6,900개 광선 대조에서 상태·수치 실패0개였다. 최대 차이는 면 경계 좌표2.80e−11 mm, 방향1.64e−12, OPL6.99e−11 mm, 상면 좌표3.98e−11 mm였다. 이는 두 구현의 일치 검사이며 물리 모델의 외부 인증은 아니다. 전체 검증 전후 입력·관련 소스 hash는 변하지 않았다.

## 7. RI·왜곡·spot의 의미

RI는 균일한 원거리 radiance의 국소 상면 irradiance에 대한 chief-map 근사다. `E/L=Aeff(θ) cosθ sinθ/[r(θ) dr/dθ]`, 축상 극한은 `Aeff(0)/f0²`를 사용한다. ray 면적 가중치와 무코팅 Fresnel 전력을 반영하며 축상 값으로 나눈다. chief-map 미분 간격을 바꾼 대조도 기록한다. 전 PSF가 검출기에서 재분포하는 과정이나 완전한 exit-pupil solid-angle 적분을 대체한 것은 아니다. RGB Fresnel RI와 파장별 RI를 모두 저장한다.

왜곡은 실제 chief 상고를 588 nm 근축 EFL의 `f tanθ` 기준과 비교한 백분율이다. 별도 파장별 값도 저장한다. Spot은 면적 가중 실제 상면 intercept의 centroid 기준 RMS radius이며, RGB 값은 공통 centroid를 사용하여 lateral color를 보존한다. radius와 diameter를 혼용하지 않는다. 보고서 요약 그림에 사용한 파장/합성 종류를 명시한다.

clear semi-diameter는 광학 면의 사용 범위다. 외경, bevel, 고정구, 사출 구배, 공차나 taper CAD는 정의하지 않았다. 단면 그림에서 면 가장자리를 잇는 선은 개략적인 연결이다.

## 8. MTF의 계산 범위와 제한

각 후보에 0/0.5/0.8/1F, 486/588/656 nm의 T/S 곡선을 0–20 lp/mm로 요청한다. 주 파동 경로는 physical-UV 역매핑 OPD와 pupil autocorrelation을 사용하는 **연속 주파수 scalar 적분**이며 실제 FFT라는 이름으로 부르지 않는다. 128/256 적분 격자 결과를 비교한다. 별도로588 nm의0/0.8F에1024 FFT 진단을 시도한다. 기존 채팅의128×128 FFT 속도 benchmark와 같은 작업량이라는 의미는 아니다.

고NA에서 이 모델은 균일 angular pupil amplitude, local chief-normal 투영을 사용하는 조건부 근사다. 단일 boundary ray-map에 fold/caustic이 있으면 `not_applicable`로 기록하고 주 파동 곡선을 만들지 않는다. 128→256 비교에서 요청 T/S 주파수 중 하나라도 상대 오차0.5%를 초과하면 `underresolved_quadrature`로 기록한다. 평균 오차가 작다는 이유로 통과시키지 않는다. `conditional`도 vector diffraction이나 Zemax와의 일치 인증은 아니다.

모든 후보에는 별도로 geometric ray-intercept OTF/MTF를 제공한다. 단색 RGB와 공통 상면 RGB 합성 곡선이 있으며, 회절을 포함하지 않는다는 상태를 유지한다. 하나의 field 패널에 RGB 중 파동 적용 불가가 있으면 패널 전체를 명확히 표시된 기하 진단으로 보여준다. 파동 실패 값을 기하 값으로 덮어쓰지 않는다. 따라서 이 후보군을 곧바로 수렴한 회절 MTF의 surrogate 정답 데이터로 사용하면 안 된다.

좋은 축상 단색 곡선과 좋은 전 시야 RGB 성능은 다르다. 이번 데이터는 후속 설계 공간 탐색을 위한 초기 후보군이며 전 후보가 우수한 완성 설계라는 보장이 없다. 서로 다른 재료·유리 위치·EPD 대안을 유지하여 DB 탐색의 출발 범위를 제공한다.

## 9. 파일과 재현

- `candidates.json`: 후보 전체의 처방, 재료 분산, 입력, 기하 출력, 파동/기하 MTF, status와 이력.
- `requirements.json`, `repair_inputs.json`: 요구조건과 이전38개 부모의 입력.
- `results/frozen_selection.json`: 파동 계산 전에 동결한23개 처방과 hash.
- `results/repair/`, `results/rescue/`, `results/stop_scan/`: 성공·실패를 포함한 변경 과정.
- `results/final_wave/`, `results/final_summary.json`: 입력hash와 연결된 파동 결과 및 집계.
- `reports/candidate_report.html`, `reports/candidate_layouts_and_mtf.pdf`, `reports/candidate_contact_sheet.png`: 전체 후보 탐색, 모든 개별 layout/MTF와 전체 단면.
- `reviews/`: 렌즈 설계, 물리, 소프트웨어 검토.

재현은 저장소 루트에서 `python experiments/candidates_ct20_20261022/native/build.py`로 native 커널을 준비한다. `prepare_candidates.py`에 qualified_existing/qualified_repair/qualified_stop_scan 폴더를 전달하고 `--count 20 --output .../results/frozen_selection.json`을 지정한다. **count는 최소 개수이지 잘라낼 개수가 아니다.** 같은 부모 중복 checkpoint만 하나로 묶고 합격한 서로 다른 부모를 모두 보존한다. `evaluate_candidates.py`에 frozen_selection을 전달하여 파동 결과를 생성한 뒤 `finalize_candidates.py`, `render_report.py` 순으로 실행한다. BLAS/OMP 스레드는1로 제한한다. 이전 원자료를 참조하는 전체 재탐색에는 이전 캠페인 폴더도 필요하다.

처방·주파수 배열·검사 상태는 JSON을 직접 읽어 후속 연산이나 DB 입력에 사용할 수 있다. 단, 입력hash가 달라지면 기존 성능 label을 재사용하지 않는다. `CT20_XXXX`는 원래 seed에서 유도한 고정 ID이고 품질 순위가 아니다.

## 10. 최종 선정·검증 요약

새 CT 제한으로 다시 계산·검증한 합격 checkpoint 26개를 독립 parent별로 묶었다. **서로 다른 seed의 합격 후보 23개를 전부 보존**했으며 상위 일부나 단일 winner를 고르지 않았다. 최종 재료 순서는 16종, 유리 위치별 수는 `{'1': 10, '2': 4, '3': 9}`다. 각 렌즈는 **2 mm < CT ≤ 20 mm**를 실제 처방과 새 수치 검사에서 모두 만족한다.

기하 screening 진단 점수는 `0.6×평균 기하 MTF 정규화 AUC + 0.4×최악 field/axis AUC − 0.25×|RGB 무코팅 RI(1F)−1|`이다. 같은 parent의 여러 checkpoint 중 하나를 고를 때만 사용하며, 서로 다른 후보들의 전체 순위나 파동 MTF 우열을 뜻하지 않는다. 표와 그림은 유리 위치·재료 순서·원본 seed ID 순서로 배치한다. `CT20_XXXX` ID는 원본 seed에서 유도하여 후보 추가나 진단 점수 변화로 바뀌지 않는다. `provenance.generation_history`는 보존된 이전 출발점·최적화 이력과 새 `ct20_continuation`/`ct20_polish`/`ct20_rescue`/`ct20_stop_scan`, `parent_input.source_result` 연결을 담는다. STOP 스캔 이력에는 원래 반경과 변경량을 보존하며, 최종 변경 후에도 전체 광선·NA·CT 조건을 다시 검사한다.

구경 모드 집계는 **front_diameter50: 22개、epd50: 1개**다. `front_diameter50`은 첫 렌즈 전면 clear diameter 50 mm 대안이며 EPD 50 mm 설계라고 부르지 않는다. 축상 NA≥0.7은 세 파장 각각 검사했지만 주변 시야의 모든 방향 cone NA≥0.7을 강제한 것은 아니다. 상고6.4 mm는588 nm chief 기준이다. 다른 파장의 lateral color와 spot이 센서 가장자리를 넘어가는 양은 데이터에 남기며, detector crop·pixel response를 적용하지 않았다.

| 집계 항목 | 최소 | 최대 |
|---|---:|---:|
| 실제 EPD [mm] | 43.00047 | 50 |
| RGB 축상 complete-rim NA | 0.7000361 | 0.7715565 |
| 실제 상면 BFL [mm] | 3.000343 | 7.705986 |
| 근축 BFL,588nm [mm] | 3.085424 | 8.686622 |
| EFL,588nm [mm] | 29.64858 | 31.32232 |
| 중심 두께 [mm] | 2.031758 | 20 |
| 축상 공기 간격 [mm] | 1.004063 | 10.14153 |
| 앞 STOP 거리 [mm] | 0.1057946 | 9.040022 |
| 첫 꼭짓점–상면 TTL [mm] | 45.83016 | 70.16402 |
| 1F 무코팅 RGB RI | 0.8614039 | 1.013422 |
| 1F 왜곡,588nm [%] | -3.871629 | 1.555047 |
| 축상 RGB spot RMS radius [mm] | 0.05409114 | 1.847196 |
| 1F RGB spot RMS radius [mm] | 0.1223166 | 1.562118 |
| 1F chief 상고,486nm [mm] | 6.30461577 | 6.35982932 |
| 1F chief 상고,588nm [mm] | 6.4 | 6.4 |
| 1F chief 상고,656nm [mm] | 6.41701513 | 6.43852857 |

파동 평가는 전체 23개 후보에서 276개의 field/파장별 T/S 쌍을 요청했고, 주 파동 곡선을 계산할 수 있었던 것은 203쌍, RGB 합성은 64field다. 적용 불가 결과는 단일 boundary ray-map의 fold/caustic 등 상세 사유와 함께 보존했다. 128→256 적분 격자 점별0.5% 기준 통과는 3/276이다. 상태 집계: `{'underresolved_quadrature': 200, 'conditional': 3, 'not_applicable': 73}`. 1024 FFT 진단은588nm의0/.8F 두field씩 46개를 계획했으나 primary 적용 조건 실패 시 실행을 생략했고, 실제 저장된 진단은 39개, 상태는 `{'not_applicable': 22, 'underresolved': 17}`다.

모든 23개 후보의 field/RGB에 대해 기하 ray-intercept OTF/MTF 진단 곡선도 따로 제공한다. 기하 단색 T/S 쌍은 276개, RGB 합성 T/S 쌍은 92개다. 파동 모델 적용 불가 패널은 이를 명시하고 기하 곡선을 표시한다. 기하 MTF에는 diffraction이 없으므로 FFT나 Huygens MTF의 계산 성공으로 바꾸어 표시하지 않는다. 실패한 파동 값을 기하 값으로 덮어쓰지 않는다.

수렴 통과도 scalar 모델·균일 angular pupil amplitude 가정 안의 적분 비교다. 최종 처방 전체의 ray-mesh·vector 고NA 물리 정확도 검증으로 확대하지 않는다. 이전 실험의 격자 probe는 이력 자료이며 CT20 최종처방의 입력hash가 다르면 새 후보 인증으로 재사용하지 않는다.

| 후보 | 시작 seed / 방법 | 재료 순서 | CT 3매 mm | EPD mm | NA588 | 실제 BFL mm | RGB RI1 | 1F RGB RMS radius μm |
|---|---|---|---|---:|---:|---:|---:|---:|
| CT20_0039 | seed_0039 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | FDS18W / PC / PC | 2.0420 / 12.6457 / 19.9965 | 43.424 | 0.71510 | 7.3378 | 0.8780 | 212.85 |
| CT20_0040 | seed_0040 / 수식·ABCD + CT20 연속 재최적화 + CT20 polish | FDS18W / PC / PC | 16.1219 / 19.9953 / 19.9995 | 43.200 | 0.71684 | 6.1202 | 0.8787 | 201.34 |
| CT20_0025 | seed_0025 / 수식·ABCD + OTF 재최적화 | FDS18W / PC / PMMA | 10.3713 / 9.4358 / 16.4169 | 50.000 | 0.76368 | 6.9015 | 0.9849 | 615.77 |
| CT20_0015 | seed_0015 / 수식·ABCD + CT20 연속 재최적화 + CT20 polish | FDS18W / PMMA / PC | 5.4216 / 17.6661 / 17.6653 | 43.200 | 0.71456 | 7.5004 | 0.8975 | 174.89 |
| CT20_0016 | seed_0016 / 수식·ABCD + OTF 재최적화 | FDS18W / PMMA / PC | 4.2845 / 18.9609 / 18.0309 | 43.239 | 0.72133 | 6.7506 | 0.8614 | 180.21 |
| CT20_0088 | seed_0088 / 수식·ABCD + CT20 연속 재최적화 + CT20 polish | TAFD45 / PC / PC | 14.7022 / 19.9986 / 19.3818 | 43.442 | 0.70630 | 3.0894 | 0.9964 | 175.86 |
| CT20_0076 | seed_0076 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | TAFD45 / PC / PMMA | 19.9981 / 19.9980 / 13.7975 | 43.253 | 0.71261 | 5.2718 | 0.9824 | 1562.12 |
| CT20_0064 | seed_0064 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | TAFD45 / PMMA / PC | 12.1977 / 19.9063 / 19.0727 | 43.200 | 0.72266 | 5.4468 | 0.8742 | 171.46 |
| CT20_0051 | seed_0051 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | TAFD45 / PMMA / PMMA | 2.0318 / 16.5667 / 19.9998 | 43.200 | 0.71317 | 7.7060 | 0.8802 | 151.38 |
| CT20_0052 | seed_0052 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | TAFD45 / PMMA / PMMA | 19.9975 / 19.9884 / 10.2203 | 43.272 | 0.70877 | 6.9427 | 0.9487 | 257.27 |
| CT20_0091 | seed_0091 / 수식·ABCD + CT20 연속 재최적화 + CT20 polish | PC / TAFD45 / PC | 19.9959 / 14.3439 / 19.9942 | 43.702 | 0.71169 | 3.3109 | 1.0128 | 207.56 |
| CT20_0092 | seed_0092 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | PC / TAFD45 / PC | 20.0000 / 20.0000 / 6.9025 | 43.538 | 0.71010 | 3.3676 | 0.9179 | 184.65 |
| CT20_0079 | seed_0079 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | PC / TAFD45 / PMMA | 19.7956 / 15.3415 / 19.9867 | 43.301 | 0.70728 | 6.4543 | 1.0134 | 188.00 |
| CT20_0080 | seed_0080 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | PC / TAFD45 / PMMA | 19.9999 / 19.9989 / 8.0495 | 43.200 | 0.71286 | 4.1231 | 0.8951 | 187.60 |
| CT20_0047 | seed_0047 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | PC / PC / FDS18W | 19.9998 / 19.8927 / 16.9093 | 43.273 | 0.70817 | 3.0232 | 0.9575 | 180.02 |
| CT20_0048 | seed_0048 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish + CT20 rescue + STOP 반경 스캔 | PC / PC / FDS18W | 19.9969 / 12.3166 / 19.9668 | 43.000 | 0.70697 | 3.0453 | 0.9602 | 210.30 |
| CT20_0036 | seed_0036 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | PC / PMMA / FDS18W | 19.9844 / 17.5071 / 19.9991 | 43.659 | 0.70994 | 3.0134 | 0.9937 | 141.25 |
| CT20_0084 | seed_0084 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish | PC / PMMA / TAFD45 | 19.9995 / 16.6620 / 19.9986 | 43.796 | 0.70883 | 3.0003 | 0.9973 | 126.70 |
| CT20_0098 | seed_0098 / 특허 변환 + CT20 연속 재최적화 + CT20 polish | PMMA / PC / FDS18W | 19.9978 / 4.1388 / 5.2119 | 43.200 | 0.70914 | 7.2075 | 0.9458 | 334.26 |
| CT20_0102 | seed_0102 / 특허 변환 + CT20 연속 재최적화 + CT20 polish | PMMA / PC / FDS18W | 19.9976 / 4.6453 / 6.6236 | 43.200 | 0.70734 | 5.9216 | 0.8980 | 163.06 |
| CT20_0072 | seed_0072 / 수식·ABCD + CT20 연속 재최적화 + CT20 polish | PMMA / PC / TAFD45 | 19.8628 / 9.0058 / 12.0662 | 43.282 | 0.71193 | 6.4485 | 0.9567 | 183.03 |
| CT20_0011 | seed_0011 / 수식·ABCD + 공통초점/scale + CT20 연속 재최적화 + CT20 polish | PMMA / PMMA / FDS18W | 19.9989 / 18.5210 / 16.3788 | 43.383 | 0.70852 | 3.0211 | 0.9599 | 122.32 |
| CT20_0059 | seed_0059 / 수식·ABCD + OTF 재최적화 + CT20 연속 재최적화 + CT20 polish + CT20 rescue + STOP 반경 스캔 | PMMA / PMMA / TAFD45 | 19.9849 / 19.8527 / 18.2069 | 43.426 | 0.70613 | 3.4514 | 0.9805 | 179.95 |

이 표와 전체 단면 그림은 합격한 독립 후보군 전체를 같은 수준으로 제시한다. 유리 위치·재료·seed ID 배치는 품질 순위가 아니다. 좋은 축상 단색 성능을 가진 후보도 주변 시야·RGB 합성은 크게 낮아질 수 있으므로 전 시야 고MTF 완성 설계로 해석하지 않는다.
