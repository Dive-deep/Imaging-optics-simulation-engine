# C2 수정 — PSF/검출기 겹침 면적의 좌표 이동 안정성

2026-10-18. 수정 대상: `src/sequential_optics/source_wave.py`의 polygon clipping, 면적, `wave_to_source_kernel`. 회귀 테스트: `tests/test_fix_c2.py`.

## Geometry가 꼬인 현상인가?

**아니다.** 렌즈나 광선 추적 없이, 한 개의 정상적인 볼록 직사각형 FFT cell과 증가 순서가 엄격하게 유지된 검출기 격자만으로 재현한다. Cell의 총 확률 질량은 1이고 검출기 전체가 이를 완전히 포함한다. Cell과 검출기를 함께 평행 이동해도 서로의 상대 배치, 면적, 순서, PSF는 같다.

- 원점 부근에서는 전체 질량이 1이다.
- 5 mm 위치로 옮기면 수정 전에는 총 확률이 1을 초과하여 `InputError`가 발생했다.
- 1000 mm 위치, cell 폭 0.0001 mm 사례에서는 수정 전 전체 포착률이 **0.8614733815193176**이었다. 검출기는 여전히 cell 전체를 포함한다.
- 새로운 테스트의 기울어진 직사각형도 실제 직교정규 basis의 가역 affine 투영으로 만든다. Self-intersection이 있는 polygon, 뒤집힌 detector edges, 렌즈 면 중첩은 사용하지 않는다.

원인은 global-coordinate shoelace 계산이다. 작은 면적을 구하면서 큰 절대 좌표들의 곱을 서로 빼므로 유효 숫자가 소실됐다. 겹침 면적을 FFT cell 면적으로 나누는 과정에서 그 오차가 확률 오차로 커졌다. 이 문제는 ray/lens geometry 오류와 별개의 수치적 좌표 의존성이다.

## 수정 내용

1. Detector edges에서 영상 reference origin을 먼저 뺀다. 작은 FFT cell 꼭짓점을 큰 전역 좌표에 더해 만들지 않는다.
2. 각 cell overlap clipping은 다시 그 cell 중심을 원점으로 하는 local frame에서 수행한다.
3. Half-plane 교점은 해당 clipping 경계 좌표를 정확히 재지정한다. 다음 경계 판정에 불필요한 교점 반올림을 전달하지 않는다.
4. Polygon area는 첫 꼭짓점을 뺀 triangle fan의 cross-product 합으로 계산한다. 절대 원점 크기의 곱을 빼지 않는다.
5. 최종 kernel을 1로 재정규화하는 보정은 제거했다. Crop된 power는 그대로 유지하고, 미세 양의 반올림은 기존 `SourceKernel`의 명시적인 확률 허용치로만 다룬다.
6. Detector grid를 계산 시작 시 검증한다. 잘못된 edges나 singular projection을 정상 geometry로 바꾸지 않는다.

변경하지 않은 것: 물리 source Q, radiance, fill, transmission, scalar diffraction PSF 정규화, 비축 detector projection의 명시적인 opt-in 및 국소 평행 투영 근사.

## 독립 기대값 및 실제 검증

수정 전 먼저 작성한 6개 테스트 메서드는 5 failures와 3 errors를 재현했다(한 메서드 내 여러 subcase 포함). 수정 후 실제 source producer 사례까지 더한 **7개 테스트가 모두 통과**했다.

| 검증 | 독립 기준 | 결과 |
|---|---|---|
| 0/5/1000 mm에서 정상 단위 cell 전체 포착 | 전체 면적이 detector 안에 있으므로 질량 1 | 통과, 절대 오차 허용 2e-14 |
| 비등방 직사각형 부분 crop | X/Y 구간의 `max(0, min(upper)-max(lower))` 곱 | 각 pixel 일치; 전체 포착 약 0.5625, 재정규화 없음 |
| 기울어진 비등방 cell | 직교정규 basis 검증, 가역 선형 변환의 중심 대칭에 따라 중심 반평면은 면적 1/2 | 전체 1, 반평면 0.5, 허용 2e-14 |
| 서로 다른 질량의 3×3 FFT cells | cell별 독립 구간 적분의 가중합 | 각 pixel 일치 |
| Polygon helper 및 순서 반전 | 실제 표현된 변 길이의 곱 | 일치 |
| 잘못된 edges/singular projection | 기존 명시 오류 계약 | 계속 거부 |
| 실제 `source_diffraction_kernel`의 전역/원점 배치 | 전체 FFT support 포함, 같은 PSF와 basis를 함께 이동 | 양쪽 전체 질량 1, pixel별 일치 |

기존 `test_source_wave.py` 7개, `test_review_source_wave.py` 4개도 모두 통과했다. 이번 담당자가 실행한 총 18개는 겹치지 않는 테스트 메서드 수다. 전체 프로젝트 suite와 독립 수정 후 리뷰는 통합 담당자가 별도로 실행한다.

실행 명령:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m unittest discover -s tests -p test_fix_c2.py -v
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m unittest discover -s tests -p test_source_wave.py -v
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m unittest discover -s tests -p test_review_source_wave.py -v
```

## 남는 수치 및 물리 한계

검출기 edges를 사용자가 전역 부동소수점 값으로 생성할 때 이미 잃은 비트는 복구할 수 없다. 테스트의 pixel별 기대값은 실제 전달된 edge 좌표를 기준으로 계산한다. 인접 edges가 같은 값으로 반올림되어 격자 폭이 0이 되면 기존 `PixelGrid` 검증이 거부한다.

실제 producer 검사는 면적 적분과 좌표 등록의 회귀 검증이다. 그 사례가 보고하는 phase sampling 상태를 물리 정확도 인증으로 승격하지 않았다. 비축 국소 투영의 광학 근사, scalar model과 상용툴 상관 검증은 이번 수정으로 달라지지 않는다.
