# 기하·근축·균일 열 모델 독립 코드 리뷰

검토일: 2026-10-18. 검토자는 metrics 구현을 담당했으며, 여기서 검토한 geometry/paraxial/materials/thermal 코드의 작성자는 아니다. 아래 평가는 실제 코드 읽기와 독립 회귀 실행에 근거한다. Zemax 대조나 실측 검증을 수행했다는 뜻은 아니다.

## 검토 범위와 결론

`geometry.py`, `paraxial.py`, `thermal.py`와 그 의존성인 `Material.expansion_scale`를 계산 명세 G1–G7, 열 형상 장과 비교했다. 실제 sag 교점, 벡터 Snell, STOP 재조준, 파면용 launch phase, EPD/F-number, 균일 팽창 계수와 위치 규칙을 확인했다.

**명시한 순차·단일 경로·균질 재료 모델 안에서 주 계산식과 연결 순서는 타당하다.** 새로 발견한 물성 유효 온도 범위 누락은 작성자에게 재현 방법을 전달해 수정했고, 별도 테스트로 수정 동작을 확인했다. 아래 두 적용 범위와 온도 baseline 규약은 반드시 결과 해석과 문서에 남겨야 한다.

## 발견 사항과 처리

| ID | 중요도 | 발견 내용 | 처리/상태 |
|---|---|---|---|
| GT-01 | 중간 | `temperature_range_c=(0,85)`가 선언되어도 `expansion_scale(200,20)`은 계산을 계속했다. index/absorption만 범위를 검사해 geometry-only 요청은 유효 범위 밖으로 외삽했다. | 물성 작성자가 비동일 온도 변화의 목표·baseline 온도를 모두 검사하도록 수정. 독립 테스트 `test_expansion_obeys_generic_model_validity_range` 통과. 정확한 T=T0 형상 항등 연산은 물성 예측을 하지 않는 별도 no-op 규약이다. |
| GT-02 | 규약 주의 | `temperature_c=None`은 각 재료 reference 온도의 n을 쓴다. geometry `reference_temperature_c`만 30으로 바꾸고 재료 reference가 20이면, `thermal_system(T0=30)`의 형상은 동일해도 index와 EFL이 변한다. | root에 재현 및 thermal sweep baseline을 동일 실제 T0로 평가할 필요 전달. 이는 현재 None 규약의 결과이며 열 팽창 식 오류는 아니다. 이후 통합 API 상태는 아래 추가 기록 참조. |
| GT-03 | 적용 범위 | `mechanical_radius`까지 동일 sag를 연장하여 최초 opaque hit 후보를 찾는다. 사용자가 처음 언급한 별도의 무효부 높이·단차·edgewall 형상은 모델 필드에 없다. | geometry 작성자 확인. 같은 sag를 연장한 불투명 rim이라는 제한을 명시해야 한다. 임의의 실제 rim 형상에 대한 ray occlusion 정확성을 주장할 수 없다. |
| GT-04 | 적용 범위 | STOP aiming은 paraxial/direct seed와 감쇠 Newton으로 하나의 국소 regular 해를 찾는다. 성공했다고 전체 정의역에 다른 aiming branch가 없음을 증명하지는 않는다. | single-branch 모델 제한. 특이 Jacobian·실패는 explicit status이며 전역 다중 branch 검출/합산은 지원 범위 밖이다. |

## 확인한 물리·수치 연결

1. **Sag와 미분**: base conic에 A2…A12를 더한다. `c=0`이어도 A 계수가 있으면 평면으로 처리하지 않는다. vertex curvature에는 `c+2*A2`가 들어가며, 실제 finite ray는 이 근축 curvature로 대체하지 않는다.
2. **최초 교점**: 유한 cylinder와 sag bound에서 전진 t 구간을 만들고, 작은 t 쪽 interval을 먼저 탐색한다. 앞선 tangent/미해결 구간을 건너뛰어 뒤쪽 root를 채택하지 않는다. 원래 sag branch의 residual을 다시 검사한다. 이 interval 구현은 formal proof 라이브러리가 아니며 arbitrary-coefficient 전역 보증으로 부를 수 없다.
3. **Active aperture**: geometric first hit를 먼저 선택한 뒤 STOP/clear radius를 검사한다. 첫 opaque hit를 무시하고 뒤의 active hit를 골라서는 안 된다는 기존 반론이 코드에 반영되어 있다. 단, GT-03의 rim 형상 제한이 있다.
4. **Snell과 OPL**: 절대 n의 접선 운동량 보존과 투과 반공간을 사용하고 TIR/critical-angle 수치 미결정을 구분한다. OPL 증가는 입사 매질 n×실제 거리다. 동일 매질 dummy surface는 불필요한 굴절을 하지 않지만 교점·개구·경로 의미는 보존한다.
5. **Launch phase**: 유한 물점은 공통 zero phase, 무한 물체는 `n*d·launch_origin`이다. 따라서 서로 다른 launch 위치의 raw geometric path 차이 때문에 가짜 경사 OPD가 생기지 않는다. 추가 테스트는 단일 공기 평면보다 강한, 비단위 ambient n을 갖는 경사 평행판에서 전체 output eikonal 식과 직접 대조했다.
6. **STOP aiming**: STOP에서 굴절 후 방향이 아니라 실제 STOP 입사 교점을 목표로 삼는다. STOP의 TIR을 조준 실패로 잘못 분류하지 않는 처리도 확인했다. 내부 STOP 앞쪽의 inverse map이 실패하면 numerical unresolved를 남긴다.
7. **근축량과 aperture**: reduced state `(y,n_abs*theta)`의 굴절·이동 행렬 determinant는 1이다. image-space EFL과 보고용 index-normalized EFL이 구별된다. F-number는 명시한 nominal infinite-conjugate `abs(report_EFL)/EPD`이며 working F/#라는 이름으로 사용해서는 안 된다. EPD를 실제 STOP 반경으로 바꾸는 데 STOP 앞 prefix matrix를 사용한다.
8. **열 형상**: `z_T(r)=s*z_0(r/s)`, `c_T=c/s`, `A_(2m,T)=s^(1-2m)*A_(2m)`을 적용한다. radius, 두께, clear/mechanical 반경은 s배, conic은 보존된다. 미분 불변량 `grad z_T(s*r)=grad z_0(r)` 및 `H_T(s*r)=H_0(r)/s`도 별도 확인했다.
9. **열 위치 규칙**: lens anchor와 고정 gap은 다른 정책이다. fixed gap과 독립 downstream anchor 운동을 동시에 강제하지 않는다. STOP radius의 follow/fixed, OBJ/IMA offset의 유지 여부도 분리되어 있다. 형상 baseline을 매 온도에 새로 사용하며 이미 가열된 geometry를 다시 팽창시키지 않는다.
10. **필요한 데이터만 요청**: 기하 팽창은 CTE만 요청하고 dn/dT·transmission은 호출하지 않는다. 순수 nominal G ray trace는 관련 없는 CTE·RGB·흡수 정보가 없어도 실행된다. 비기준 온도의 실제 굴절 계산은 n(T) 또는 명시적 invariance를 요구한다.

## 추가한 독립 테스트와 근거

`tests/test_review_geometry.py`에 아래 6개를 추가했고 전부 통과했다. 기존 작성자 테스트와 같은 output snapshot을 복제하지 않고 다른 물리 항등식을 선택했다.

| 테스트 | 독립 기대값/불변량 | 검출 목적 |
|---|---|---|
| 경사 평행판 output eikonal | `S_exit = n0*d_xy·x + sqrt(ng²-n0²*|d_xy|²)*thickness` | launch phase와 OPL의 연결 오류, absolute/relative n 혼용 |
| 6면 비구면 처방의 정확한 역경로 | 마지막 실제 ray를 뒤집으면 원래 모든 교점과 입사 방향으로 복귀 | 단순 반대쪽 plane-wave bundle과 정확한 ray reversal 혼동, Snell branch/normal 오류 |
| 다른 입출력 n의 단일 굴절면 | R=50, n0=1.2, n1=1.6 → image EFL=200 mm, ambient-normalized EFL=150 mm; 실제 작은 ray limit=200 mm | index gauge와 물리 초점거리의 혼동, nominal F/# aperture 환산 오류 |
| CTE 온도 유효 범위 | 0–85°C 자료로 200°C 비항등 팽창을 요구하면 실패 | 발견 GT-01 회귀 |
| 열 비구면 미분 similarity | gradient 불변, Hessian 1/s | A2…A12의 지수 오류가 surface normal/굴절까지 전파되는 문제 |
| 전체 3매 assembly similarity | 명시 anchor 운동으로 모두 s배 → 위치/OPL s배, 방향 동일, EFL s배 | 개별 sag 단위 테스트를 통과해도 geometry 배치/재료/trace 연결이 틀리는 문제 |

실행: `PYTHONPATH=src python -m unittest discover -s tests -p test_review_geometry.py -v`.

## 남는 한계

- 정확한 Maxwell 전파, 비균일 온도, 열응력, mount 변형을 이 검증으로 확인한 것이 아니다.
- standalone thermal radial-separation 검사는 sample+local refinement이며 proof certificate가 아니다. 실제 tracing에서는 geometry의 clear-aperture interval 검사가 추가된다.
- 미세한 접촉/근접 접선·ill-conditioned pupil에서는 값을 억지로 만드는 대신 미결정 상태를 유지하는 것이 맞다. 그런 상태를 physical vignetting 또는 성능 0으로 합치면 안 된다.
- 실제 제품과의 수치 일치는 vendor 설정, pupil measure, 공기/파장 규약을 맞춘 독립 대조가 필요하다.

참조: `docs/calculation_spec_ko.md`의 G1–G7, thermal geometry/CTE 장과 해당 장의 primary-source URL 목록. 위 새 테스트의 기대식은 Snell, eikonal continuity, 균일한 공간 similarity에서 직접 유도했으며 합성 재료 수치를 catalog 값으로 사용하지 않았다.
