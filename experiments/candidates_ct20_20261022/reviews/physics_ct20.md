# CT ≤20 mm 최종 후보군 독립 광학 물리 검토

검토일: 2026-10-22 UTC. 검토자: `candidate_physics`. 대상은 이 폴더의 동결된 23개 후보 전부이며, 기존 `../candidates_20261022/`의 처방·코드·24개 감사 결과는 수정하거나 새 후보의 검증으로 재사용하지 않았다.

최종 입력은 `results/frozen_selection.json`, SHA-256은 `19a6f2e21fff64b2a36947426ffb42a75683617fb3fd41971f7fcf396464497d`다. 아래 검증은 계산상의 요구조건을 만족하는 초기 후보군임을 확인한다. 완성 설계의 고MTF, 제조 공차, 전 시야 NA 또는 Zemax와의 일치를 인증하는 결과는 아니다.

## 1. 새로 수행한 독립 검증

`validate_selection_physics.py`를 CPU 프로세스 하나, BLAS/OMP/MKL 스레드 하나로 실행했다. 저장된 합격 표시만 읽지 않고, 후보마다 0–12°의 13개 field와 486/588/656 nm에서 STOP 16 rings ×48 arms에 chief·rim을 더한 817개 광선을 다시 추적했다.

- 후보 23개, 서로 다른 부모 23개 모두 통과했다.
- 새 광선 수는 **732,849개**이며 실패 광선은 **0개**다. 유효 광선만 남겨서 성공 비율을 다시 정의하지 않았다.
- RGB 각각의 축상 complete rim 전체가 통과했고, 그중 최소 실제 출사 NA가 각각 0.7 이상이다.
- 표에 내보낸 radius, A4–A12, conic, material-after, vertex, thickness, clear semi-diameter, STOP 및 상면으로 native 면 데이터를 재구성했다. 최대 면 데이터 차이는 **6.9389×10⁻¹⁸**이다.
- 저장된 chief 상고, RI, 왜곡, 투과 면적 비율, RMS radius를 현재 처방에서 다시 계산한 값과 대조했다. 최대 차이는 **0**이었다.
- 각 면 sag domain과 인접 면 사이 간격을 공통 clear aperture의 반경 4,097점에서 독립 sag 식으로 확인했다. 최소 인접 면 간격은 **1.0040630077 mm**였다. 이 유한 표본 검사는 연속 영역의 수학적 증명을 대신하지 않는다. 원본 Python 엔진의 별도 연속 radial-interval geometry 검사는 `reviews/software_ct20.md`에 있다.
- 무한 물체, degree 단위 0/6/9.6/12°, 전체 FOV 24°, RGB 파장, 이미지 지름 12.8 mm의 export 계약을 확인했다.
- 같은 부모의 checkpoint 진단 점수는 수식과 대조하되, 이전의 전체 score 내림차순 요구는 제거했다. 표시 순서는 유리 위치·재료 순서·원본 seed이며 전체 rank나 단일 winner는 없다.

근거는 `results/final_selection_physics_audit.json`과 `results/final_selection_physics_statistics.json`이다.

## 2. 중심 두께의 엄격한 상·하한

각 매의 축상 중심 두께는 **2 < CT ≤20 mm**를 정확한 부동소수점 비교로 검사한다. 상한에 +ε를 더하거나 평균 두께로 판정하지 않는다. 최적화 continuation의 임시 상한이 20을 넘더라도 최종 gate는 바뀌지 않는다.

`validate_ct20_boundary.py`의 15개 public-evaluator control은 세 렌즈 각각에 대해 다음을 확인했다.

- 정확히 20 mm는 상한 통과, 바로 다음 표현 가능한 값 `20.000000000000004` mm는 실패한다.
- 20보다 바로 아래 값은 통과한다.
- 정확히 2 mm는 열린 하한에 실패하고, 바로 위 값 `2.0000000000000004` mm는 통과한다.
- `optimization_ct_upper_mm=35`를 넣어도 최종 판정은 완화되지 않는다.

최종 CT 범위는 **2.0317579927–19.9999970367 mm**다. 최대 CT와 20의 차이는 약 2.96 nm에 불과하다. 이는 저장된 명목 처방의 수치 통과이며, 그 작은 여유가 제조 공차를 허용한다는 의미는 아니다. Boundary control JSON에 저장된 selection-auditor hash는 표시 순서 검사를 바꾸기 전 버전이며, CT evaluator는 바꾸지 않았다. 최종 처방 전체는 변경된 감사 스크립트로 별도 재검사했다.

## 3. 구경, NA, 상고와 BFL의 해석

22개는 **첫 렌즈 전면 clear diameter 50 mm** 대안이고, 1개는 **EPD 50 mm** 대안이다. EPD 범위는 43.000474–50 mm다. EPD50 후보의 첫 전면 clear diameter는 약 53.199181 mm이며, 50 mm opening 대안과 혼동하지 않는다. 전면 STOP 앞에는 광학 파워가 없으므로 nominal STOP 지름이 입사동 지름이며, 이번 표본에서는 축상 전체 원형 rim과 pupil이 실제로 통과한다.

NA는 공기 중 출사 방향으로부터 `sqrt(dx²+dy²)=sin(u)`를 구하며, 가장 큰 한 광선이나 근축 `1/(2 F-number)`로 합격시키지 않는다. 축상 complete-rim 최소 NA의 588 nm 범위는 **0.706134–0.763681**, RGB 전체 범위는 **0.700036–0.771556**다. 최소 NA 여유가 작으므로 재료 분산 모델·제조 변화까지 견딘다고 주장할 수 없다.

이 기준은 **축상 NA 조건**이다. 1F의 588 nm tangential opposite-rim cone NA는 **0.663625–0.778230**이며, 일부 후보는 주변 시야에서 0.7보다 작다. 사용자가 향후 전 시야 T/S cone NA≥0.7을 요구하면 현재 통과 판정과는 다른 제약으로 다시 탐색해야 한다.

588 nm의 1F chief 상고는 모든 후보에서 6.4 mm이고 최대 절대 잔차는 8.88×10⁻¹⁶ mm다. 0.001 mm는 계산 gate의 수치 근접 기준이며 사용자 승인 설계 오차가 아니다. 다른 파장에서는 lateral color가 남는다.

| 항목 | 최소 | 최대 |
|---|---:|---:|
| 1F chief 상고, 486 nm [mm] | 6.30461577 | 6.35982932 |
| 1F chief 상고, 656 nm [mm] | 6.41701513 | 6.43852857 |
| 실제 고정 상면까지 vertex BFL [mm] | 3.00034343 | 7.70598620 |
| 588 nm 근축 BFL [mm] | 3.08542414 | 8.68662198 |
| 588 nm EFL [mm] | 29.64858322 | 31.32231665 |
| 축상 공기 간격 [mm] | 1.00406301 | 10.14153434 |

적색 chief는 모든 후보에서 +6.4 mm를 넘는다. spot 전체가 센서 안에 들어간다는 뜻이 아니며 detector crop을 적용하지 않았다. BFL 두 값은 서로 다른 정의라 별도로 검사했다. 모든 렌즈는 공기로 분리되어 있고, PMMA/PC 두 장과 구면 standard 유리 한 장이며 유리에 conic이나 A4–A12를 넣지 않았다. 유리 위치는 앞 10개·중간 4개·뒤 9개, 재료 순서는 16종이다.

## 4. RI, spot과 왜곡의 물리적 범위

RI 식은 균일 원거리 radiance의 chief-map 국소 근사다.

`E(θ)/L = Aeff(θ) cosθ sinθ / [r(θ) dr/dθ]`, 축상 극한은 `Aeff(0)/f0²`이다.

`Aeff`는 STOP 면적 가중치와 투과 광선의 Fresnel 전력을 반영한다. 막힌 광선이 있다면 그 면적을 그대로 잃으며 surviving rays만으로 면적을 재정규화하지 않는다. 반면 spot의 centroid와 RMS는 정의상 투과광 분포의 조건부 모멘트다. RGB spot은 하나의 공통 centroid를 사용하므로 lateral color를 지우지 않는다.

1F RGB RI는 기하 모델에서 **0.860721–1.037504**, 무코팅 scalar Fresnel을 넣으면 **0.861404–1.013422**다. RI>1은 기준 축상보다 국소 image mapping이 집중되는 경우 가능하며 값을 1로 잘라내지 않았다. 다만 이 결과는 각 field의 PSF가 주변 detector 위치로 재분포하는 전 과정을 적분한 값이 아니다. 특히 큰 aberration이 있는 후보에서 실제 센서 irradiance의 정답으로 쓰기 어렵다. 고 NA의 exit-pupil projected-solid-angle 적분이나 완전한 radiometric image simulation과의 검증이 추가로 필요하다.

단색 왜곡은 해당 파장의 근축 EFL과 `f tanθ`를 기준으로 계산한다. 보고서의 대표 왜곡은 588 nm다. 그 1F 범위는 **−3.871629%…+1.555047%**이다.

RGB RMS spot radius는 축상 **54.09–1847.20 μm**, 1F **122.32–1562.12 μm**다. 1F 반경이 1.56 mm인 후보까지 포함되므로 모든 후보를 고해상도 완성 렌즈라고 표현하면 부정확하다. 조건을 만족하는 서로 다른 설계 공간의 출발점을 모두 보존한 데이터로 해석해야 한다.

## 5. MTF와 학습 label에 대한 비판적 판단

기하 OTF는 실제 상면 intercept의 복소 Fourier 합으로서 탐색에 유용하지만 diffraction MTF가 아니다. `abs(OTF)`의 AUC만 높다고 전 주파수에서 좋은 대비를 보장하지 않는다. 일부 곡선은 영점 이후 다시 올라가며 실수 OTF가 음수가 되는 구간도 있다. 예를 들어 CT20_0025와 CT20_0076에는 각각 133개와 153개의 음수 실수 OTF 표본이 기록되어 있다. 따라서 AUC 진단 점수는 같은 부모의 checkpoint 선택에 제한하고, 후속 최적화에서는 curve 전체와 복소 OTF, 저주파 대비·T/S·RGB·전 시야를 함께 봐야 한다.

이번 고 NA에서는 균일 angular pupil amplitude와 scalar 회절, chief-normal 국소 투영을 사용한 파동 경로도 조건부 모델이다. 물리적 boundary ray map이 fold를 만들면 단일값 pupil 보간을 그대로 사용할 수 없다. 해당 파동 결과의 `not_applicable`을 기하 MTF로 성공 처리해서는 안 된다. 격자 비교의 점별 0.5% 통과 역시 그 scalar 모델 안의 quadrature 비교이며, 광선 mesh 수렴·벡터 회절·실측·Zemax와 0.5% 이내라는 뜻은 아니다.

최종 보고서의 conditional/underresolved/not_applicable 상태 및 별도 기하 fallback 표시는 이 구분을 유지한다. 현재 후보의 입력과 기하 출력은 DB 탐색 기록으로 저장할 수 있지만, 전체를 신뢰할 만한 고 NA diffraction-MTF surrogate의 정답 label로 바로 사용하기에는 검증이 부족하다.

## 6. 재현 명령

저장소 루트에서 다음을 실행한다. 결과 파일은 기존 증거를 덮어쓰지 않도록 배타적으로 생성하므로 재실행 시 새로운 output 이름을 사용한다.

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python experiments/candidates_ct20_20261022/validate_selection_physics.py experiments/candidates_ct20_20261022/results/frozen_selection.json --output experiments/candidates_ct20_20261022/results/physics_audit_repeat.json
```

동일 source와 frozen 처방에서 다시 계산하는 감사이며, 생성·최적화·파동 계산 전체 시간의 benchmark는 아니다.
