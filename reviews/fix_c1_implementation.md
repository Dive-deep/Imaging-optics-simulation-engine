# C1 수정 — 동일 곡면의 순차 교점 재사용

대상 소스: `src/sequential_optics/geometry.py`. 신규 회귀: `tests/test_fix_c1.py`.

## 원인

일반 교점 솔버가 허용오차 내에서 구한 점에는 0이 아닌 sag residual이 남을 수 있다. 그 점에서 동일한 곡면을 다시 추적하면, 기존 first-root 솔버는 `t=0`에서 정확한 부동소수점 0 또는 부호 변화를 찾는다. residual 부호에 따라 실제로는 같은 경계에 있는 광선을 이미 경계를 지난 광선으로 판단했다. 따라서 검증에 통과한 합법적 동일 곡면의 무굴절 dummy가 광선을 소실시켰다.

## 수정

- 직전 면의 교점이 정상 인증되었고 다음 면의 vertex z, curvature, conic, A2…A12가 **정확히 동일한 경우**에 한해 그 교점을 거리 0으로 재사용한다. 생략된 다항식 항과 명시적 0은 같은 sag로 비교한다.
- 새 면의 mechanical domain, residual, 비접선 조건을 확인한다. Active aperture 및 현재 매질에 따른 Snell 검사는 기존 추적 절차대로 수행한다.
- 같은 규칙을 STOP aiming의 가상 prefix에도 적용한다. 이 가상 경로에서는 원래대로 physical aperture가 aiming residual을 차단하지 않으며, 최종 실제 추적은 physical boundary를 적용한다.
- STOP이 동일 sag의 연속 chain 안에 있으면 aiming 위치 제약은 그 chain의 첫 면을 기준으로 한다. 후속 독립 리뷰가 찾은 rear STOP 사례에서, 아직 수렴 전인 trial ray가 첫 rear 면에서 목표 z보다 뒤에 닿는 것을 별도 prefix 실패로 처리하던 문제도 수정했다. 같은 위치에서의 연속 굴절은 그 위치를 바꾸지 않으므로 이 기하 사건의 첫 교점을 인증하면 된다. 실제 추적은 chain의 모든 면을 그대로 처리한다.
- 일반 `intersect`의 판정은 바꾸지 않았다. 약간 다른 면, 임의 외부 광선, 작은 residual만 있는 점을 동일 면이라고 간주하지 않는다. 전역 허용오차도 넓히지 않았다.

교점의 z를 sag에 맞춰 다시 투영하지 않아 광선 경로를 바꾸거나 인공 OPL을 더하지 않는다. 무굴절 dummy에서는 좌표·방향·OPL이 유지되며, 굴절률이 실제로 다른 동일 면은 거리 0인 정상 Snell 경계로 처리한다.

## 검증

실행:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:examples python3 -m unittest tests.test_fix_c1 tests.test_geometry tests.test_regression_r3_r4 tests.test_review_geometry -v
```

결과: **54개 테스트 통과**, 22.227초. 이 중 신규 C1 테스트는 **13개**다. 후속 독립 리뷰의 rear STOP 반례를 포함한 최종 재실행이다. 전체 테스트는 통합 작업자가 별도로 실행한다.

- 원래 sphere singlet의 `Field(0,0)`, `Field(.01,.02)`, `Field(0,.1)` 각각에서 원본·dummy 추가본 모두 **217/217** 정상 추적. Pupil mask/weight, surface positions, directions, OPL, normals, spot centroid/RMS/표본 최대 반경 및 transmitted fraction 비교.
- 원래 실패 사례의 양의 작은 residual을 실제로 확인한 후, 동일 교점을 좌표 투영 없이 그대로 유지함을 검사.
- 내부 STOP 이전 prefix, 동일 면 자체가 STOP인 경우, 유한 물체와 A2…A12 even-asphere를 검사. 서로 다른 Newton aiming 경로의 허용오차는 테스트 비교에 반영했다.
- Front 3중·rear 2중 동일 sag chain에서 rear STOP을 마지막 duplicate로 옮긴 뒤 비축 두 field와 4개 pupil을 회귀 검사. 별도 독립 리뷰 제어에서는 front/rear STOP, 3개 field, 5개 pupil의 **30개** 비교를 통과했고 위치 차이 최대 6.11e-16 mm, 방향 차이 2.78e-17, OPL 차이 3.12e-17 mm였다.
- Dummy의 더 작은 active aperture 및 mechanical domain을 각각 여전히 적용하는지 검사.
- 동일 경계의 두 번 굴절 `1 → 1.2 → 1.5`가 단일 `1 → 1.5`의 Snell 방향과 일치하며 추가 OPL이 0인지 검사.
- 1e-6 mm 떨어진 다른 면과 같은 z이나 A2가 다른 면은 양의 전파거리와 OPL을 유지하는지 검사.
- 작은 residual의 독립 ray를 임의로 경계에 붙이지 않는지, 실제 tangency 및 첫 번째 교점을 유지하는지 검사.
- 기존 R3/R4의 물리 clipping/가상 aiming/curved STOP rim, 원래 first-root/tangency 및 독립 역추적·열 similarity 대조 모두 통과.

## 범위

이 수정은 선언된 동일 sequential sag에 대한 수치 일관성 수정이다. 서로 다른 처방 표현이 우연히 같은 곡면을 만드는 모든 대수적 동치까지 탐색하지 않는다. 예를 들어 곡률 0일 때 무관한 conic 값이 다르더라도 별도 면으로 다룬다. 큰 허용오차를 사용해 임의의 외부 ray를 surface에 자동 스냅하는 기능은 추가하지 않았다.
