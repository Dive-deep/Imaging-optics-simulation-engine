# 50 mm 개구·24° 화각·3매 광학계 초기 후보 생성 과정

작성일: 2026-10-22 (Asia/Seoul). 이 문서는 후보 처방의 생성·수정·선정 과정과 계산 범위를 설명한다. 후보 데이터는 `candidates.json`, 그림과 곡선은 `reports/candidate_report.html` 및 `reports/candidate_layouts_and_mtf.pdf`에 있다. 데이터베이스 서버를 생성하거나 MySQL에 저장하지 않았다.

이 후보들은 이후 설계 공간 탐색을 시작하기 위한 실제 수치 처방이다. 전체 설계 공간의 최적해, 제조 가능한 완성품, Zemax와 일치가 검증된 설계 또는 수렴이 확인된 학습 정답이라는 뜻은 아니다. 기하 조건 통과와 MTF의 모델·수치 품질을 별도 상태로 저장한다.

## 1. 요구조건을 수치 입력으로 옮긴 방법

| 항목 | 구현한 정의 |
|---|---|
| 물체 | 무한 원거리의 평행 광선, 회전대칭 중심 광학계 |
| 화각 | Y 방향 0°, 6°, 9.6°, 12° = 0, 0.5, 0.8, 1F; 회전대칭에 의해 음의 화각 대칭 |
| 상면 | 모든 field와 파장이 공유하는 하나의 평면, 기준 상고 ±6.4 mm |
| 1F 위치 | 기준 파장 588 nm의 chief ray Y좌표 6.4 mm; 잔차와 수치 판정 허용값을 별도 저장 |
| 파장 | 사용자가 지정한 486, 588, 656 nm. Material의 nd는 별도로 587.5618 nm 기준 |
| 매수 | 공기로 분리한 3매: PMMA/PC 두 장 + FDS18W 또는 TAFD45 한 장 |
| 구경 대안 A | 앞 STOP 지름 50 mm = EPD 50 mm. STOP 앞에 광학적 파워가 없음 |
| 구경 대안 B | 첫 번째 렌즈의 전면 clear diameter 50 mm. STOP/EPD 지름은 별도 변수 |
| 두께·공기 간격 | 각 렌즈 축상 중심 두께 >2 mm, 두 렌즈 간 축상 air gap ≥1 mm |
| BFL | 마지막 면 꼭짓점에서 실제 고정 상면까지 ≥3 mm; 별도로 근축 BFL ≥3 mm도 확인 |
| NA | 실제 굴절 후 광선의 image-space NA. 최종 선정은 축상 전체 rim의 최소값이 486/588/656 nm 각각 ≥0.7인지 확인. 주변 field의 cone 값은 별도 보고 |
| STOP | 첫 렌즈 앞의 평면. 앞 꼭짓점까지 거리도 최적화. 전체 첫 clear surface보다 앞인지 검사하고 EFL보다 가까운지 기록 |
| 플라스틱 면 | standard 또는 even_asphere. 이번 비구면은 conic과 A4, A6, A8, A10, A12만 사용 |
| 유리 면 | standard 구면. k=0, 비구면 계수 모두 0을 추적 입력 단계에서 강제 |
| 열·흡수·코팅 | 열해석 없음. 벌크 흡수는 모르는 값으로 추정하지 않음. 광량 비교는 손실 없는 기하값과 무코팅 Fresnel 근사를 구분 |

첫 렌즈 지름 50 mm와 EPD 50 mm는 같은 조건이 아니다. 후보마다 `aperture_mode`, 실제 지원되는 EPD, STOP 지름, 첫 전면 clear diameter를 함께 저장한다. 일부 광선만 통과하는데 nominal STOP 지름을 유효 EPD라고 보고하지 않는다.

초기 최적화는 계산 효율을 위해 작은 pupil 표본에서 시작하고, 후보 선정은 STOP 16 rings ×48 arms에 chief와 rim을 추가한 817개 광선/field/파장으로 다시 평가한다. 구경 보정과 연속 시야 점검에는 0–12°의 13개 field ×3파장 ×817개 =31,863개 광선을 사용한다. 결과는 연속 pupil 전체에 대한 수학적 극값 증명은 아니다.

## 2. 먼저 검토한 광학적 난이도

Rectilinear 기준에서 필요한 EFL은

`f = 6.4 / tan(12°) = 30.1096327006621 mm`

이다. EPD 50 mm이면 근축 f-number는 약 0.6022다. 그러나 `1/(2 F-number)`를 실제 고NA로 사용하지 않고 추적한 출사 방향으로 NA를 구한다.

EPD 50 mm와 원형 24° 전체 시야를 손실 없이 image circle 12.8 mm로 전달한다는 가정에서는 étendue 용량상 등가 수용 NA가 `25 sin(12°)/6.4 ≈0.8122` 이상 필요하다. 이것은 모든 처방에 무조건 적용되는 축상 NA 공식이 아니라, 전 시야 광량·pupil 가정을 둔 난이도 검토다. 첫 렌즈 지름 50 mm를 고정하고 EPD를 작게 잡는 대안은 분리해 탐색할 이유가 있다.

양의 플라스틱 두 장과 음의 고굴절 유리 한 장은 색수차를 줄이기 좋은 출발점이지만 Petzval 곡률이 남는다. 반대로 양의 구면 유리와 음의 PC를 섞는 구조는 다른 색수차·상면 만곡 절충을 제공한다. 따라서 한 가지 렌즈 배열의 미세 변형만 만들지 않았다.

## 3. 재료 모델

| 재료 | 사용자 nd | 사용자 Vd |
|---|---:|---:|
| PMMA | 1.492 | 57.47 |
| PC | 1.586 | 32.00 |
| FDS18W | 1.94594 | 17.98 |
| TAFD45 | 1.95375 | 32.32 |

사용자가 채택한 normal-line 근사를 유지했다. `PgF = 0.6438 − 0.001682 Vd`와 `nF−nC=(nd−1)/Vd`를 사용하고, `n(λ)=a+b/λ²+c/λ⁴`의 세 계수를 nd·F−C·g−F 조건으로 결정했다. λ는 이 식에서 μm다. 표준선은 d=587.5618, F=486.1327, C=656.2725, g=435.8343 nm다. 요청 파장 486/588/656에서의 n을 이 모델로 평가하여 JSON 재료 표에 저장한다.

이 값은 실제 재료의 측정 Sellmeier 데이터나 catalog의 anomalous partial dispersion을 복원한 것이 아니다. 특히 유리 이름을 근거로 사용자 nd/Vd를 다른 catalog 값으로 교체하지 않았다. 단위와 공기/진공 기준은 기존 사용자 모델의 nominal absolute-index 해석을 그대로 유지한다.

## 4. 112개 출발점의 생성

### 4.1 수식 기반 96개

재료 조합, 유리의 앞·가운데·뒤 위치, 양/음 파워 배열, 면 bending과 구경 해석을 바꿨다. 근축 출발 파워는 `Σφ≈1/f`와 `Σφ/V≈0`의 색수차 균형을 참고했다. 이 식을 실제 광선 추적 대신 사용하지 않았다.

각 굴절면은 reduced ray `(y,n u)`에 대한 `R=[[1,0],[-(n_after−n_before)c,1]]`, 간격은 `T=[[1,t/n],[0,1]]`로 표현한다. 여섯 면과 다섯 내부 간격의 행렬 `M=[[A,B],[C,D]]`에서 `EFL=−1/C`, `BFL=−A/C`를 구하고, 세 렌즈의 초기 파워를 비선형 보정해 목표 EFL 근방으로 맞췄다. 렌즈 두께는 실제 sag를 고려해 중앙과 주변에서 면이 겹치지 않도록 준비했다.

초기 구조의 상세 파워·bending·보정 잔차와 반복 횟수는 각 seed의 `provenance`에 있다. 이 단계에서 얻은 EFL은 30.10921–30.11053 mm 범위지만, 근축 일치가 실제 NA·spot·MTF 합격을 의미하지 않는다. seed_0069의 음의 근축 BFL 같은 실패도 원자료에서 삭제하지 않았다.

### 4.2 공개 특허 기반 16개

[US3944337A — Three element objective lens with ellipsoidal, hyperbolic, and parabolic surfaces](https://patents.google.com/patent/US3944337A/en)의 4개 수치 예를 출발 구조로 사용했다. 이 특허는 앞조리개, 두 플라스틱 렌즈의 네 conic 면, 뒤의 구면 유리 렌즈를 제공한다. 원래 f/1.4–1.57로 이번 목표와 같지 않다.

수행한 변경은 다음과 같다.

1. 원문 표와 도면에서 반경/곡률·conic·중심 두께·공기 간격·STOP 위치를 전사했다.
2. 원래 focal length에서 목표 약 30.11 mm로 길이를 비례 조정했다.
3. 원래 polystyrene을 허용된 PC로, 원래 유리를 사용자 FDS18W 또는 TAFD45로 대체했다. 대체 후 원래 광학 성능이 유지된다고 가정하지 않았다.
4. 두께 >2 mm와 air gap ≥1 mm를 만족하도록 수정하고, exact ABCD 파워 보정을 다시 수행했다.
5. EPD 50 mm와 첫 렌즈 clear diameter 50 mm 대안으로 나눈 뒤 실제 광선 기반 구경 확대 최적화를 적용했다.

채택한 원문 예는 conic만 사용하므로 12차 초과 항을 삭제할 필요가 없었다. 최종 후보에는 A2·홀수차·14차 이상 항이 없다. 다른 공개 특허도 조사했지만 표면 종류나 목표 구경이 달라 직접 처방으로 채택하지 않은 이유를 `references/RESEARCH_KO.md`에 기록했다. 특허 원자료를 완성된 사용자 후보로 제시하지 않는다.

## 5. 실제 광선 최적화 과정

원본 v0.2.1 코어는 보존했다. 별도 실험의 FP64 C 커널로 6면의 실제 conic/다항 sag 교점, 벡터 Snell, OPL과 광선 이력을 구했다. 첫 전방 교점을 구간으로 분리하고 safeguarded Newton을 사용한다. clipping, TIR, 교점 실패, 뒤쪽 진행을 서로 다른 상태로 남긴다. 실패 광선을 버려서 작은 spot을 만들지 않는다.

최적화 변수는 6면 곡률, 플라스틱 4면 conic, 최대 20개 A4–A12 항, 세 중심 두께, 두 air gap, 상면 거리, STOP 위치 및 허용된 경우 STOP 반경이다. 수치 조건을 개선하기 위해 asphere는 최적화 중 `b_m(r/25)^m`으로 표현하고, export할 때 `A_m=b_m/25^m`로 변환한다. Radius는 curvature의 역수로 저장한다.

실제 box bound는 플라스틱 곡률 ±0.1/mm, 유리 곡률 ±0.98/clear-radius, conic −40…2, 각 정규화 b4…b12 −30…30mm, 중심 두께2.02…35mm, air gap1…20mm, 실제 상면 BFL3…50mm, 앞 STOP 거리0.1…28mm다. Opening50 대안의 STOP 반경은21.6…24.8mm다. 이 bounds는 이번 수치 탐색의 범위이며 허용 설계 공간 전체가 아니다. 예를 들어 더 작은 EPD나 더 두꺼운 렌즈를 탐색하지 않았다는 제한이 남는다.

1차 탐색은 112개 모두에 대해 STOP 반경의 0.45→0.60→0.72→0.82→0.90→0.96→1.00 비율로 구경을 확대했다. 처음에는 A4/A6를 사용하고 마지막 단계에서 A8/A10/A12도 풀었다. 실제 RGB 광선의 공통 평면 잔차, 588 nm full-field chief height, NA, 기하 간격과 실패 광선 penalty를 함께 최소화했다. 각 단계의 잔차 평가 횟수·종료 이유·시간·처방을 저장했다. optimizer의 종료 코드만으로 광학 합격을 판정하지 않았다.

후속 개선에서는 실제 pupil 면적 가중치, 같은 RGB centroid, 0/.5/.8/1F를 사용했다. 1.3, 2.9, 4.7, 6.8, 9.1, 12.3, 16.2, 20 lp/mm에서 **geometric complex OTF**를 빠른 merit의 일부로 사용하고, 낮은 real OTF와 음의 응답 및 RMS를 함께 벌점화했다. 이것은 회절 FFT를 최적화했다고 표현하지 않는다. 마지막 평가에서는 별도의 scalar 파동 곡선을 계산한다.

1차 단계는3rings×8arms와chief/rim에서 각stage 최대100번의least-squares 함수평가를 허용했고, 상대종료 tolerance는2e−5였다. 유한차분 Jacobian 때문에 실제 residual 호출 수는 `nfev`보다 크며 로그에 둘 다 남겼다. 후속 OTF 개선은6×16, 최대220 nfev, 종료 tolerance1e−6, finite-difference step5e−6을 사용했다. 최종 순위와 합격 판정은16×48에서 새로 계산했다. 낮은 표본의 merit에서 좋아진 값이 최종 표본에서도 자동으로 좋아진다고 가정하지 않았다.

추가로 상위8개 parent에서 모든 RGB/field가 공유하는 상면을−0.5…+0.7mm의49점으로 이동했다. 각점에서 전체 처방을 similarity scaling하여588nm/12° chief 상고6.4mm를 복원하고, 첫 전면 clear radius25mm를 다시 고정한 후 실제 광선을 재추적했다. 이 과정은 EPD50 고정 모드에는 적용하지 않았다. 평균/최악 AUC와 RI 편차의 보존 조건을 통과한 개선본2개만 남겼다. `focus_similarity_best/scan_summary.json`과 최종 후보의 history에 실제 shift·scale·전후 값을 기록했다.

MTF20 한 점만 최대화하지 않았다. OTF가 영점을 지나 부호가 바뀐 후 MTF 크기가 다시 커지는 처방은 유용한 연속 대역 해상도가 좋다고 볼 수 없기 때문이다. 곡선 면적·낮은 T/S·반전 지점과 RI 편차를 함께 저장한다.

다만 비대칭 spot의 OTF는 복소수이므로 real OTF가 음수라는 사실만으로 복소 OTF의 크기가0을 통과했다고 단정할 수 없다. JSON의 negative-real-OTF 점수는 위상·contrast response를 살피는 진단이다. 최종 magnitude-AUC 순위가 모든 주파수에서 단조롭게 높은 MTF를 보장하지 않는다.

## 6. Clear aperture 조정과 독립 검사

일부 시작점은 실제 광선 footprint보다 훨씬 큰 후방 clear aperture를 가지고 있었다. 모든 RGB와 13개 field를 다시 추적하고, 각 면의 최대 footprint에 `max(0.30mm,2%)` 여유를 둬 사용하지 않는 aperture만 줄였다. 처음 50 mm opening 조건의 첫 전면은 25 mm 반경을 고정했다.

보정 전후의 통과 광선·상태·전체 이력·방향·OPL을 비교한다. 이미 통과하던 광선을 하나라도 제거하면 해당 조정을 거부한다. 새 aperture에서 유효 광선만 다시 뽑아 통과율을 좋게 보이는 처리가 아니다. 구면 domain과 상면/앞 STOP 간섭도 다시 검사한다.

최종24개는 별도 물리 검토에서13field×RGB×817×24=764,712개 광선을 새로 추적해 실패0을 확인했고,4097개 반경 표본의 기하 간격과 export된 반경/계수/물성/두께의 역변환도 대조했다. 별도 소프트웨어 검토는24×300=7,200개의 chief·중간·rim 광선을 원본 Python 교점/Snell/OPL과 비교했다. 최대 경계 좌표 차2.93e−11mm, OPL차7.01e−11mm, 상면 좌표차3.66e−11mm였고 상태 불일치0이었다. 원본 방사 구간 geometry 검사도24/24 통과했다. 이는 검사하지 않은 모든 연속 광선이나 파동 모델 정확도의 증명은 아니다.

앞뒤 면의 clear semi-diameter는 독립 변수다. 이 데이터는 광학적으로 사용하는 면의 범위를 명시하며, 외경·bevel·고정구·사출 구배·중심 두께 공차나 완성된 taper CAD를 정의하지 않는다. layout의 면 끝을 잇는 선은 개략적인 단면 연결이다.

최종 조건 판정은 모든 렌즈의 재료/면 종류, 축상 두께·간격, 실제 상면 BFL와 근축 BFL, conic domain, radial 표본에서의 면 겹침, STOP/상면 간섭, 전체 rim의 NA, 상고 잔차, RGB/field 전체 광선 실패 수를 별도로 기록한다. 단순히 평균값이 좋아서 실패 후보를 합격시키지 않는다.

## 7. RI·왜곡·spot의 정의

RI는 균일한 원거리 radiance가 만드는 국소 영상면 irradiance를 근사한다. field 각도를 θ, chief image radius를 r(θ), 입사측 STOP의 유효 면적을 Aeff라 하면

`E/L = Aeff(θ) cosθ sinθ / [r(θ) dr/dθ]`

이고 축상 극한은 `Aeff(0)/f0²`다. chief map의 미분은 서로 다른 작은 간격으로 대조한다. Aeff에는 원래 pupil 면적 가중치와 실패/차단을 반영하며 surviving ray 수로 재정규화하지 않는다. 기하 RI와 여섯 면의 무코팅 scalar Fresnel을 적용한 RI를 둘 다 저장한다. 사용자의 coating·광원 spectrum·detector response가 없으므로 이를 임의의 제품 특성으로 추정하지 않았다.

Rectilinear `r=f tanθ`와 vignetting 없는 조건은 cos⁴ 법칙으로 돌아가며 12°에서 약 0.9154다. Sine mapping은 기하 RI가 평탄해지는 출발점이다. 그러나 실제 후보의 RI는 이 이론값을 복사하지 않고 실제 chief map·통과 pupil·Fresnel로 계산한다. 수차로 인한 전체 PSF의 검출기 재분포나 OpticStudio의 완전한 exit-pupil solid-angle 적분과 같다는 검증은 없다. [공식 RI 정의](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Relative_Illumination.html)를 참고했다.

Spot의 기본 숫자는 centroid 기준 **RMS radius**이며 단위는 mm 또는 그림의 μm다. RMS diameter도 따로 저장한다. RGB 합성 spot은 파장별 중심을 따로 맞추지 않고 공통 RGB centroid로 계산해 lateral color를 보존한다. 왜곡은 588 nm 근축 EFL의 rectilinear 기준과 실제 chief mapping을 명시하고, RI용 국소 chief slope와 혼동하지 않는다.

MTF와 RI는 lens-only 결과다. 12.8 mm 센서의 pixel pitch/fill factor, 검출기 crop, 입사각 response 또는 실제 촬영 장면을 포함한 영상 성능은 아니다.

## 8. 고NA MTF를 계산하고 표시한 방법

이번 NA는 일반적인 낮은 NA FFT 근사보다 높다. [Ansys의 FFT와 Huygens PSF 비교](https://optics.ansys.com/hc/en-us/articles/42661795535251-What-is-the-difference-between-the-FFT-and-Huygens-PSF)도 FFT 모델의 빠른 광학계·pupil distortion 적용 한계를 설명한다. 따라서 scalar 결과를 vector Debye/Maxwell 해나 Zemax 수치 일치로 표현하지 않는다.

주 곡선은 같은 scalar Fourier pupil 모델의 연속 주파수 autocorrelation이다.

`OTF(f)=∫P(χ−λf/2) P*(χ+λf/2) dχ / ∫|P(χ)|² dχ`

파장·field별 OPD를 실제 광선의 OPL과 기준 구면에서 구하고, OPD 자체를 보간한 후 위상차를 만든다. 복소 위상을 먼저 보간해 phase wrap을 섞지 않는다. 실제 FFT를 실행한 곡선은 별도 `actual_scalar_FFT` diagnostic으로 저장한다. 주 곡선을 실행하지 않은 FFT의 결과라고 부르지 않는다.

128 pupil grid는 NA~0.7에서 image window가 작아 0–20 lp/mm에 몇 개의 독립 주파수만 생길 수 있다. Padding만 늘리는 것으로 해결되지 않는다. 연속 주파수 autocorrelation은 요청 주파수를 직접 평가하고, 별도의 FFT 진단에는 1024 grid를 사용한다.

높은 NA에서 물리 boundary의 삼각형을 direction-cosine 평면에 직선으로 옮기면 실제 좌표 변환은 일대일이어도 직선 삼각형이 겹칠 수 있었다. 이를 해결하기 위해 원래 UV 삼각형에서 exact inverse `UV=Lχ/sqrt(1−|χ|²)`로 OPD를 조회한다. 적분 측도와 선언된 uniform-angular amplitude는 유지한다. 보간 좌표가 바뀌었으므로 예전 numerical profile과 bitwise 같은 계산이라고 주장하지 않는다.

고차항이 상쇄되는 후방 비구면은 항별 거친 상계가 실제보다 지나치게 큰 boundary를 만들 수 있다. 반경을 32,768개 구간으로 나누고 conic과 각 다항항의 단조 끝점 상계를 합해, 면 전체 뒤에 놓이는 보수적 평면을 구한다. 별도의 analytic/경계 control로 검토했으며 원본 코어 파일은 바꾸지 않았다.

RGB 합성은 같은 물리 좌표로 등록한 **complex OTF를 합한 뒤 절댓값**을 취한다. 파장별 MTF 크기를 평균하거나 RGB 중심을 각각 이동시키지 않는다. 현재 파동 합성은 파장별 정규화 image energy를 동일하게 둔다. RI의 equal-incident-radiance + Fresnel weighting과는 별도 정의이며 사용자 spectrum/detector response가 주어지면 바꿔야 한다.

파동 결과에는 서로 다른 적분 격자 간 점별 변화와 `conditional`, `underresolved_quadrature`, `not_applicable`를 그대로 저장한다. 0에 가까운 MTF에서 큰 상대 차이가 생겨도 평균으로 숨기지 않는다. 미수렴 값은 초기 탐색의 진단 곡선으로 볼 수 있지만, 바로 neural surrogate의 확정 학습 정답으로 넣으면 안 된다. Traced ray mesh의 별도 refinement도 대표 후보에 대해 기록한다.

최종 C01의0.8/1F,588nm에서 발견된 적용 불가는 별도로16×48→32×96→64×192 ray mesh까지 검사했다. Canonical boundary의 음의 방향 삼각형 수는0.8F에서8→24→80,1F에서44→168→662로 남았다. 마지막 S6면의 실제 XY map에도 가장 조밀한 표본에서 각각30/622개의 fold가 있었고, sampled exit footprint 바로 뒤로 평면을 옮겨도 사라지지 않았다. 따라서 이 사례는 단순 보간격자 부족이나 angular 직선 삼각형 문제로 보지 않는다. 단일값 spatial pupil에 대한 현재 모델의 적용 조건이 실패한 것이다. 가드를 해제하거나 접힌 triangle을 삭제하지 않았다. 다른 모든 실패 사례의 caustic 위치까지 같은 것으로 단정하지 않는다. 근거는`results/wave_recovery_C01_diagnosis.json`, 재현 코드는`wave_recovery.py`다.

후속 회절 엔진에는 출사 방향의 angular spectrum 또는 Huygens/multibranch 전파를 검토할 수 있다. 그 경우 위상 변환, 진폭 측도, caustic 위상 및 높은 NA의 vector/polarization 효과를 새로 검증해야 한다. 이번 데이터에서 미적용/미수렴 상태를 합격으로 바꾸는 대체식으로 즉시 사용하지 않았다.

## 9. 재현·데이터베이스로 옮기는 방법

후보 JSON의 surface table에는 radius, thickness-to-next-surface, material-after, clear semi-diameter, 표면 종류와 허용된 비구면 계수가 실제 mm 단위로 있다. STOP와 image plane은 따로 명시한다. 마지막 surface의 thickness는 실제 상면까지의 간격이다. 내부 `design`은 같은 처방을 계산 코드로 재현하기 위한 정규화 표현이다.

후보 간 보간은 재료 조합과 surface topology별로 나누는 편이 좋다. Radius의 무한대/부호를 직접 보간하기보다 curvature와 정규화한 asphere 값을 사용하고, 새 처방마다 geometry→통과 광선→NA/상고/BFL→wave 품질 순서로 검사한다. 두 합격 처방 사이를 보간했다고 중간 처방이 자동으로 합격하지 않는다. 실패 이유도 별도 레코드로 보존하면 DB에서 feasibility 모델과 성능 모델을 구별할 수 있다.

원본 엔진은 GitHub에서 가져온 v0.2.1이다. 이전 세션의 미업로드 속도 최적화 파일은 새 환경에 없었으며, 이번 결과를 그 파일에서 재현했다고 주장하지 않는다. 이 실험 폴더의 native build·seed 생성·최적화·평가 script와 source SHA를 보존했다. GPU는 사용하지 않았다.

### 실행 예

```bash
python experiments/candidates_20261022/native/build.py
python experiments/candidates_20261022/render_report.py experiments/candidates_20261022/candidates.json --output reports_regenerated
```

`seeds.py`의 생성 함수는 보존된 `seeds.json`을 만드는 원자료다. 기존 결과를 덮어쓰지 않는 새 output 디렉터리를 사용해 `optimize.py`, `refine.py`, `qualify.py`를 실행한다. 단일 후보를 다시 평가할 때는 저장된 `design`을 `radiometry.evaluate_design()`와 `wave_evaluation.evaluate_design()`에 전달한다. 새로운 성능 수치는 새 결과 파일로 저장한다. 현재 환경의 원본 코어277개 test도 모두 통과했으며, 새 기록은`results/restored_core_tests.json`에 별도로 보존했다.

## 10. 최종 선정·검증 요약

112개 출발점에서 탐색한 뒤 38개의 독립 seed가 기하 조건을 통과했다. 선정 시점의 76개 합격 checkpoint 중 같은 seed의 여러 수정본을 하나로 묶고, **서로 다른 seed의 24개**를 선택했다. 최종 재료 순서는 15종이며 유리는 앞·가운데·뒤 위치를 모두 포함한다.

선정 점수는 `0.6×평균 기하 MTF 정규화 AUC + 0.4×최악 field/axis AUC − 0.25×|RGB 무코팅 RI(1F)−1|`이다. 사용자 선호를 대신해 정한 초기 탐색용 휴리스틱이며 물리적으로 유일한 점수나 미수렴 파동 MTF의 우열 증명은 아니다. 각 후보의 `provenance.generation_history`에는 시작 처방, 7단계 최적화 로그, 재최적화 또는 공통 초점 보정 이력이 들어 있다.

**24개 모두 첫 렌즈 전면 clear diameter 50 mm 대안을 사용한다. EPD 50 mm 설계라고 부르면 안 된다.** 축상 NA≥0.7은 세 파장 각각 검사했지만 주변 시야의 모든 방향 cone NA≥0.7을 강제한 것은 아니다. 상고6.4 mm는588 nm chief 기준이다. 다른 파장의 lateral color와 spot이 센서 가장자리를 넘어가는 양은 데이터에 남기며, detector crop·pixel response를 적용하지 않았다.

| 집계 항목 | 최소 | 최대 |
|---|---:|---:|
| 실제 EPD [mm] | 43.2 | 43.95938 |
| RGB 축상 complete-rim NA | 0.7017721 | 0.7329018 |
| 실제 상면 BFL [mm] | 3.00067 | 9.027529 |
| 근축 BFL,588nm [mm] | 3.014 | 9.617752 |
| EFL,588nm [mm] | 29.65378 | 30.75316 |
| 중심 두께 [mm] | 2.069588 | 33.93733 |
| 축상 공기 간격 [mm] | 1.000064 | 7.804999 |
| 앞 STOP 거리 [mm] | 0.1012229 | 2.9688 |
| 첫 꼭짓점–상면 TTL [mm] | 53.52227 | 75.64402 |
| 1F 무코팅 RGB RI | 0.8562326 | 0.988907 |
| 1F 왜곡,588nm [%] | -2.092566 | 1.53725 |
| 축상 RGB spot RMS radius [mm] | 0.0709827 | 0.169949 |
| 1F RGB spot RMS radius [mm] | 0.1146283 | 0.4280942 |
| 1F chief 상고,486nm [mm] | 6.29290221 | 6.35915865 |
| 1F chief 상고,588nm [mm] | 6.39999006 | 6.4005068 |
| 1F chief 상고,656nm [mm] | 6.41729966 | 6.44311678 |

파동 평가는24×4field×3파장=288개의 T/S 쌍을 요청했고, 주 파동 곡선을 계산할 수 있었던 것은 236쌍, RGB 합성은 78field다. 나머지는 단일 boundary ray-map의 fold/caustic 조건으로 모델 적용 불가 상태를 보존했다. 128→256 적분 격자 점별0.5% 기준 통과는 11/288이다. 상태 집계: `{'underresolved_quadrature': 225, 'conditional': 11, 'not_applicable': 52}`. 1024 FFT 진단은588nm의0/.8F 두field씩48개를 계획했으나 primary 적용 조건 실패 시 실행을 생략했고, 실제 저장된 진단은 42개, 상태는 `{'underresolved': 33, 'not_applicable': 9}`다.

모든24개×4field×RGB에 대해 기하 ray-intercept OTF/MTF 진단 곡선도 따로 제공한다. 파동 모델 적용 불가 패널은 이를 명시하고 기하 곡선을 표시한다. 기하 MTF에는 diffraction이 없으므로 FFT나 Huygens MTF의 계산 성공으로 바꾸어 표시하지 않는다. 실패한 파동 값을 기하 값으로 덮어쓰지 않는다.

수렴 통과도 scalar 모델·균일 angular pupil amplitude 가정 안의 적분 비교다. 최종 처방 전체의 ray-mesh·vector 고NA 물리 정확도 검증으로 확대하지 않는다. 대표 처방의16×48→32×96 비교는 `results/wave_mesh_probe.json`에 있다. seed0012의588nm축상은 점별 기준을 통과했지만 다른 대표 주변field는 실패했다. 그 probe와 최종처방의 입력hash가 다르면 최종후보 인증으로 재사용하지 않는다.

| 후보 | 시작 seed / 방법 | 재료 순서 | EPD mm | NA588 | 실제 BFL mm | RGB RI1 | 1F RGB RMS radius μm |
|---|---|---|---:|---:|---:|---:|---:|
| C01 | seed_0059 / 수식·ABCD + OTF 재최적화 | PMMA / PMMA / TAFD45 | 43.663 | 0.71176 | 3.0076 | 0.9724 | 134.91 |
| C02 | seed_0100 / 특허 변환 + 공통초점/scale | PMMA / PC / TAFD45 | 43.441 | 0.70877 | 3.6731 | 0.9389 | 125.19 |
| C03 | seed_0011 / 수식·ABCD + 공통초점/scale | PMMA / PMMA / FDS18W | 43.383 | 0.70865 | 4.1363 | 0.9449 | 114.63 |
| C04 | seed_0036 / 수식·ABCD + OTF 재최적화 | PC / PMMA / FDS18W | 43.698 | 0.70711 | 3.0313 | 0.9692 | 171.26 |
| C05 | seed_0084 / 수식·ABCD + OTF 재최적화 | PC / PMMA / TAFD45 | 43.959 | 0.70771 | 3.0361 | 0.9889 | 166.32 |
| C06 | seed_0012 / 수식·ABCD + OTF 재최적화 | PMMA / PMMA / FDS18W | 43.200 | 0.71276 | 3.0396 | 0.9296 | 152.58 |
| C07 | seed_0108 / 특허 변환 + OTF 재최적화 | PMMA / PC / TAFD45 | 43.396 | 0.70821 | 3.0231 | 0.9712 | 152.68 |
| C08 | seed_0071 / 수식·ABCD + OTF 재최적화 | PMMA / PC / TAFD45 | 43.422 | 0.71911 | 3.6616 | 0.9239 | 164.77 |
| C09 | seed_0047 / 수식·ABCD + OTF 재최적화 | PC / PC / FDS18W | 43.373 | 0.70818 | 3.0310 | 0.9458 | 196.09 |
| C10 | seed_0096 / 수식·ABCD + OTF 재최적화 | PC / PC / TAFD45 | 43.850 | 0.71324 | 3.0038 | 0.9750 | 222.53 |
| C11 | seed_0095 / 수식·ABCD + OTF 재최적화 | PC / PC / TAFD45 | 43.540 | 0.70970 | 3.0135 | 0.9480 | 195.69 |
| C12 | seed_0048 / 수식·ABCD + OTF 재최적화 | PC / PC / FDS18W | 43.506 | 0.70907 | 3.0213 | 0.9539 | 195.73 |
| C13 | seed_0068 / 수식·ABCD + OTF 재최적화 | PMMA / TAFD45 / PC | 43.200 | 0.72308 | 4.3242 | 0.9132 | 180.28 |
| C14 | seed_0083 / 수식·ABCD + OTF 재최적화 | PC / PMMA / TAFD45 | 43.549 | 0.71039 | 3.0007 | 0.9282 | 428.09 |
| C15 | seed_0072 / 수식·ABCD | PMMA / PC / TAFD45 | 43.203 | 0.70917 | 9.0275 | 0.9169 | 167.83 |
| C16 | seed_0087 / 수식·ABCD + OTF 재최적화 | TAFD45 / PC / PC | 43.200 | 0.71919 | 4.1391 | 0.9414 | 223.57 |
| C17 | seed_0079 / 수식·ABCD + OTF 재최적화 | PC / TAFD45 / PMMA | 43.338 | 0.71413 | 7.2974 | 0.9022 | 187.87 |
| C18 | seed_0092 / 수식·ABCD + OTF 재최적화 | PC / TAFD45 / PC | 43.216 | 0.70976 | 4.0745 | 0.9229 | 225.86 |
| C19 | seed_0051 / 수식·ABCD + OTF 재최적화 | TAFD45 / PMMA / PMMA | 43.200 | 0.71464 | 7.6566 | 0.8798 | 148.16 |
| C20 | seed_0056 / 수식·ABCD + OTF 재최적화 | PMMA / TAFD45 / PMMA | 43.208 | 0.71498 | 5.1778 | 0.8914 | 181.18 |
| C21 | seed_0020 / 수식·ABCD | PMMA / FDS18W / PC | 43.200 | 0.72058 | 6.0375 | 0.8562 | 147.52 |
| C22 | seed_0016 / 수식·ABCD + OTF 재최적화 | FDS18W / PMMA / PC | 43.239 | 0.72133 | 6.7506 | 0.8614 | 180.21 |
| C23 | seed_0080 / 수식·ABCD + OTF 재최적화 | PC / TAFD45 / PMMA | 43.200 | 0.71212 | 5.6816 | 0.9215 | 203.08 |
| C24 | seed_0015 / 수식·ABCD | FDS18W / PMMA / PC | 43.200 | 0.71486 | 7.5346 | 0.8983 | 176.64 |

이 표의 순서는 위 기하 screening 점수 순서이며 최종 파동 성능 순위가 아니다. 좋은 축상 단색 성능을 가진 후보도 주변 시야·RGB 합성은 크게 낮아질 수 있다. 따라서 이번 결과를 “24개 모두 전 시야 고MTF 달성”으로 해석하지 않는다.
