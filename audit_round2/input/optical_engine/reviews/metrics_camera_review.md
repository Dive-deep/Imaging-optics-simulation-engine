# Independent camera-design review of geometric metrics

검토 역할: 결상광학·카메라 광학 설계 관점. 검토 대상은 `metrics.py`, 해당 API 연결, `tests/test_metrics.py`이다. 이 검토자는 해당 metrics 모듈을 작성하지 않았다. 추가 검증은 `tests/test_review_metrics.py`에 남겼다.

## 결론

현재 구현에서 spot/실광선 fan/국소 T·S 초점/왜곡/공통 기하 초점/Zernike 진단의 계산 정의와 수식은 명시된 독립 광학 모델에 일관된다. 실광선 결과를 EFL로 치환하지 않고, centroid와 chief, local focus와 동시에 달성 가능한 common focus, 단순 even fan과 Seidel coma 계수를 구분한다. 이 평가는 알고리즘·해석 검토와 아래 물리 fixture 결과이며 Zemax 일치나 실제 렌즈 측정 검증이 아니다.

## 발견하여 반영한 사항

1. **API–fan 인자 불일치.** `evaluate`가 `angular=`를 넘기는데 기존 `ray_fans`가 받지 않아 분석 요청이 `TypeError`로 끝나는 상태를 확인했다. 모듈 작성자가 afocal 실광선 fan 경로를 추가했다. 추가 회귀시험은 focal/afocal 두 모드의 실제 API 요청을 통과한다.
2. **유한 각도의 T/S 좌표변환.** chief-local 두 atan2 각도를 global XY field T/S 행렬로 선형 회전하면 일반적인 비축 2축 field에서 정확한 meridional angular fan이 아니다. 검토 의견에 따라 작성자가 실제 3D tangent frame T3/S3를 만들고 각 출사 방향에서 `atan2(d·T3,d·dc)`, `atan2(d·S3,d·dc)`를 직접 평가하도록 수정했다. 원래 XY angular chart도 별도로 보존한다. 독립 단일 평면 굴절 fixture로 유한각 T/S를 확인했다.

## 독립 정답이 있는 fixture

실제 렌즈의 정답 데이터 없이도 검증할 수 있도록, 평면 굴절면의 벡터 Snell 법칙을 독립 미분했다. 축 방향 물체 거리 a, 물체 높이 h, 굴절률비 η=n1/n2, `ci=a/sqrt(a²+h²)`, `ct=sqrt(1-η²(1-ci²))`라 하면 surface plane 기준의 parabasal 초점은 다음이다.

- Tangential: `delta_T = -a ct³/(η ci³)`.
- Sagittal: `delta_S = -a ct/(η ci)`.

이 경우 영상은 virtual이며 유한 입사각에서는 평면 굴절면에도 T/S virtual focus 차이가 생긴다. a=10 mm, h=5 mm, η=1/1.5에서 실제 aiming·Snell·finite-difference field curvature가 위 값과 허용 상대오차 2e-8/절대오차 1e-7 mm 이내로 일치한다. field를 0.71 rad 회전한 경우 동일한 T/S 초점을 얻는다. 이 fixture는 구면 렌즈 초점 예제를 되풀이하는 대신 부호·실광선 미분·기준면·T/S 분리를 동시에 확인한다.

## 기존 시험에서 확인한 방어 장치

- 기하 spot의 가중치는 pupil quadrature이며 LED radiance나 Fresnel을 몰래 넣지 않는다. 전달된 가중치와 실제 통과 비율이 결과에 남는다.
- chief와 centroid 기준 RMS가 다른 경우를 구분한다. 물리 clipping과 수치 실패를 혼합하지 않는다.
- common focus는 field마다 pupil weight를 먼저 정규화한 뒤 field weight를 적용한다. 한 field의 ray 개수/weight 전체 스케일이 목적함수를 임의로 바꾸지 않는다. 모든 field를 한 IMA에서 재평가한다.
- distortion은 선언한 input chart와 작은 field의 실제 chief Jacobian을 사용한다. nominal map을 따로 보존하므로 온도별 배율을 다시 맞추어 왜곡값이 0이 되어도 원래 배율 drift가 숨지 않는다.
- Zernike는 full-disk RMS 정규화된 실제 polynomial의 weighted SVD fit이다. 불규칙/차단 pupil에서는 rank와 condition을 남기며 제거한 defocus 등의 값은 원래 OPD를 변경하지 않는 진단값이다.
- LSA는 실제 축 교차이고 parabasal focus와 별개다. Coma even/odd descriptor는 Seidel 계수라는 이름으로 보고되지 않는다.

## 남는 해석 범위

- 공통 기하 RMS focus는 FFT MTF 최적 focus와 다르다. 사용자가 MTF 기준 최적화를 원하면 그 목적함수를 별도로 평가해야 한다.
- 가상 IMA 연장은 spot/field-curvature 진단으로 유효하지만 실제 detector irradiance가 아니다. `virtual_intercepts` 및 `allow_virtual` 구분을 유지해야 한다.
- 각도 spot의 RMS는 명시한 chief-local atan2 chart의 RMS이며 구면 geodesic RMS 또는 제품의 모든 afocal 표기와 동일하다고 볼 수 없다.
- 한 field의 T/S 국소 초점은 작은 pupil 주변의 differential result이며 전체 비구면 pupil이 동시에 만드는 최적 spot을 대신하지 않는다.
- 분석 request의 `status=ok`는 계산 함수 반환 상태이다. rank deficiency·coupling·virtual geometry 등 내부 상태도 사용자에게 보여야 한다.

## 실행 기록

`PYTHONPATH=src python -m unittest tests.test_review_metrics tests.test_metrics -q`

2026-10-18: **28 tests passed**, 1.473 s. 새 독립 검토 시험은 4개이다. 기존 시험을 포함한 전체 엔진 합격 여부는 root의 최종 integration gate에서 별도로 기록한다.
