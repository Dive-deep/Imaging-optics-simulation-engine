# 2차 독립 광학 물리 리뷰

검토자 역할: 최초 구현에 참여하지 않은 광학 물리학자. 먼저 릴리스 코드를 읽고 독립 반례를 작성했으며, 기존 리뷰나 183개 통과 결과를 정확성의 근거로 삼지 않았다.

검토 대상: `input/optical_engine`, 배포 archive SHA256 `d854cbaba525581088fe7b890ee54e2f1fb663c8ef8c6fb562b06a66a65284ef`. 원본 코드는 수정하지 않았다. 아래 경로는 이 프로젝트의 상대 경로이며 줄 번호는 동결된 릴리스 기준이다.

## 판단

스칼라 회절, OPL, 광량 전달의 기본 식에서 무작위로 식을 조합한 흔적은 발견하지 않았다. 특히 입체각 측도와 전기장 진폭을 구분하는 접근은 합리적이다. 그러나 **이 상태를 신뢰할 수 있는 일반 광학 데이터 생성기로 승인할 수는 없다.** 이번 검토에서 재현된 두 종류의 결함은 각각 원거리 MTF를 거의 0에서 거의 1로 바꾸거나, 설계를 바꿨는데도 예전 MTF를 정상적인 결과처럼 반환할 수 있다.

`conditional`은 이미 수렴 검증 전이라는 제한을 포함한다. 아래 지적은 이를 `validated`로 오해한 데서 나오지 않는다. 현재 구현에는 검사 가능한 요청 주파수의 해상도 부족과 명백한 입력 데이터 불일치를 감지하지 못하는 구체적인 경로가 있다.

| ID | 우선순위 | 분류 | 확실성 | 요약 |
|---|---|---|---|---|
| PH-01 | P1 | 구현·수치 검증 결함 | 재현 완료 | 원거리 angular OTF가 요청 주파수에 필요한 각도 해상도를 검사하지 않아 심각한 alias를 반환 |
| PH-02 | P1 | 데이터 일관성 결함 | 재현 완료 | 다른 처방에서 생성한 ray bundle을 현재 처방의 것으로 받아들임; STOP 변경은 기존 boundary cache 검사도 통과 |

P1은 이 경로를 설계 탐색·학습 데이터 생성에 사용하기 전에 수정해야 한다는 뜻이다. 모든 정상 입력의 모든 계산이 틀렸다는 뜻은 아니다.

## PH-01 — angular MTF 요청 주파수에 대한 quadrature 해상도 검사 누락

**위치**

- `src/sequential_optics/wave.py:204–219`, `angular_otf`
- 같은 파일 `559–580`, afocal 호출 및 최종 상태 판정
- `src/sequential_optics/source_wave.py:210–218,226–228`, 동일 함수를 공유하는 afocal source response

**원인**

코드는 다음 이산 적분을 한다.

\[
\operatorname{OTF}(f)=\sum_j p_j\exp(-2\pi i f\alpha_j),\qquad \sum_jp_j=1.
\]

이 식 자체는 적절한 midpoint quadrature이다. 문제는 요청한 `f`에 비해 `alpha_j` 간격이 너무 클 때도 그대로 계산한다는 점이다. pupil의 위상이 잘 샘플링되었는지와 이 적분의 진동 인자가 잘 샘플링되었는지는 서로 다른 조건이다.

작은 각도에서 FFT 각도 간격은 대략

\[
\Delta\alpha\simeq\Delta l=\frac{\lambda/n}{N\,\text{padding}\,\Delta u}.
\]

`f≈1/Δα`이면 이산 합의 위상들이 다시 정렬되어 거짓으로 MTF가 1 가까이 될 수 있다. 정확한 `atan2` 좌표를 사용하는 것만으로 이 문제가 없어지지 않는다. 특히 `analyze_wave`는 afocal 결과에 `linear_autocorrelation_supported`가 없으면 `True`로 간주하며, 요청 angular frequency의 sampling 조건은 확인하지 않는다.

**재현 A: 매끄러운 Gaussian 물리 경계장**

`repro_angular_otf_alias.py`와 `repro_angular_otf_alias.json`에 전체 코드와 출력이 있다.

- 동일한 64×64 경계장, 간격 0.01 mm, 진공 파장 550 nm, n=1
- 전기장 `A(u,v)=exp(-(u²+v²)/0.06²)`
- 요청 주파수 4654.54545 cycles/rad

| padding | 반환 MTF | Fourier 창 가장자리 에너지 비율 |
|---|---:|---:|
| 4 | 0.9999999931 | 2.60e-28 |
| 8 | 1.17e-17 | 1.12e-28 |
| 16 | 1.36e-16 | 5.33e-29 |

경계장 샘플은 그대로이고 angular output quadrature만 세분화했다. edge energy가 작은지 검사하는 것으로도 이 alias를 발견하지 못한다.

독립적인 **소각도** 연속 Gaussian 식은

\[
\operatorname{MTF}(f)=\exp\left[-\frac{\lambda^2 f^2}{2w^2}\right].
\]

해당 조건에서는 지수가 약 −910이다. 이를 exact-angle MTF의 엄밀한 정답이라고 주장하지 않는다. exact-angle 계산 자체에서도 padding 8과 16은 거의 0으로 일치하며, padding 4의 거의 1과는 명백히 다르다. 각도로 좌표를 바꾼 후에는 일반적인 유한 pupil의 cosine-space cutoff를 그대로 정확한 angular cutoff로 적용해서는 안 된다.

**재현 B: 실제 ray bundle → analyze_wave 전체 경로**

반지름 0.3 mm의 투명한 공기 평면 aperture, 축상 평행 광선, 표본 경계장 진폭을 입력했다. wave grid 64, traced pupil 8 rings ×32 arms를 사용했다. 경계 위상 차이는 0이다.

| 조건 | padding 4 | padding 8 | padding 16 |
|---|---:|---:|---:|
| 4581.81818 cycles/rad MTF | 0.999999466 | 4.19e-8 | 2.97e-8 |
| **기본 0–10000/121점 축에 포함된 4583.33333 cycles/rad MTF** | **0.999913461** | **9.40e-8** | **6.65e-8** |
| 반환 status | conditional | conditional | conditional |
| pupil phase max step | 0 | 0 | 0 |

이 사례는 사용자 원거리 광학계의 angular MTF 경로에 직접 관련된다. Gaussian을 선택한 이유는 aperture 불연속, 고차 수차, pupil phase alias와 분리하여 출력 적분 자체의 문제를 드러내기 위해서다.

**예상 동작 / 수정 방향**

1. 요청된 최대 angular frequency와 **실제 비균일 angular 좌표**에서의 최대 또는 허용 잔여 power를 고려한 위상 증가량 `2π f Δα`, `2π f Δβ`를 검사한다.
2. 기준을 만족하지 못하면 padding을 증가시키거나, 해당 주파수 결과를 `underresolved`로 표시하거나 명확하게 거절한다.
3. 비균일 각도 셀 적분을 더 정확히 수행하는 방법도 가능하다. 그러나 cell-center 값을 단순히 점 질량으로 놓고 임의의 높은 주파수까지 신뢰해서는 안 된다.
4. padding 4/8/16 비교, 입력 pupil grid 변화, angular window 변화는 서로 다른 검증 항목이다. 하나의 수렴 시험으로 대체하지 않는다.
5. source-wave afocal 경로도 같은 검사를 공유해야 한다.

단순한 `abs(f)>고정 상수` 검사로 해결할 문제는 아니다. `atan` 변환으로 angular grid가 비균일하며, 전체 angular window와 실제 에너지 지지영역에 따라 필요한 해상도가 달라진다.

## PH-02 — ray bundle과 처방의 일치 여부를 보장하지 못함

**위치**

- `src/sequential_optics/models.py:132–143`, `RayBundle`에 생성 처방 provenance 없음
- `src/sequential_optics/wave.py:320–332`, `_fingerprint`
- 같은 파일 `354–427`, `make_boundary_state`
- 같은 파일 `478–485`, boundary 재사용 검사
- `src/sequential_optics/source_wave.py:135–165`, supplied bundle 검사

### PH-02a: 새 boundary를 만들면 변경된 처방과 오래된 rays의 조합을 승인

`_fingerprint`는 trace가 생성될 때 기록한 처방을 검사하는 대신, **지금 전달된 system과 지금 전달된 bundle을 함께 해싱**한다. 따라서 old bundle + changed system으로 새 상태를 만들면 그 모순된 조합 자체에 새로운 정상 fingerprint가 붙는다.

릴리스의 synthetic triplet에서 마지막 surface curvature만 `0.995`배로 바꾸고, 명목 처방의 ray bundle을 그대로 전달했다. 이미지 위치는 고정했다.

| 60 lp/mm MTF(S) | 값 | status |
|---|---:|---|
| 원래 처방 + 원래 rays | 0.5081472574 | conditional |
| 변경 처방 + 원래 rays | **0.5081472574** | conditional |
| 변경 처방 + 재추적 rays | 0.0532438410 | conditional |

두 처방에서 나온 정밀한 실제 MTF의 절대 정확도를 이 비교가 보증하는 것은 아니다. 여기서 확인한 결함은 **렌즈가 바뀌어도 이전 ray histories를 현재 처방의 것으로 아무 검사 없이 수용한다는 것**이다. 재추적 결과는 이 입력 오류의 영향을 보여준다.

`source_diffraction_kernel(..., bundle=old_bundle)`도 받아들인다. 해당 함수는 첫 surface에서 출발점·Snell·OPL을 검사하지만, 이후 surface가 달라진 경우는 검출하지 못한다. 이 경우 Fresnel transmission은 바뀐 prescription의 normal/index를 사용하면서 direction/OPL은 이전 trace 값을 사용할 수 있어 물리적으로 일관되지 않은 혼합 데이터가 된다.

### PH-02b: STOP radius 변경은 기존 cached boundary fingerprint 검사도 우회

이 variant는 총괄 리뷰어가 코드에서 발견했고 본 리뷰에서 재현했다. `_fingerprint`에는 `stop_index`, `stop_radius`가 없다. image plane 이동을 허용하려고 `image_z`를 제외한 것은 의도적이지만, STOP 변경까지 제외하는 것은 다른 문제다.

- synthetic triplet STOP radius 4 mm →2 mm
- 기존 bundle + 기존 `boundary_state`를 그대로 전달
- fingerprint가 동일하므로 허용되고 이전 MTF가 그대로 나온다.

| 60 lp/mm MTF(S) | 값 | status |
|---|---:|---|
| STOP 4 mm의 cached 결과를 STOP 2 mm로 사용 | **0.5081472574** | conditional |
| STOP 2 mm에서 재추적 | 0.2265950279 | conditional |

`repro_stale_bundle.py`와 `repro_stale_bundle.json`은 두 variant 및 physical-source API 수용 여부를 재현한다.

**영향 / 수정 방향**

- 설계 탐색, 온도 sweep, cache 재사용 때 조용히 잘못된 레이블을 생성할 수 있다. 검토한 기본 `evaluate` 호출은 내부적으로 새 bundle을 만들므로 이 경로를 직접 밟지 않는다. 공개 low-level 모듈을 조합하거나 재사용할 때의 위험이다.
- ray tracing 시 생성 처방 fingerprint를 bundle에 저장하고, 해석 API에서 현재 처방과 비교해야 한다. STOP index/radius, aperture, geometry, 사용 파장, 해당 파장에서 평가된 매질 index, 실제 온도, OBJ/field와 필요한 trace 설정을 포함한다.
- IMA만 이동할 수 있다는 예외는 명시적으로 분리한다.
- 사용자가 직접 만든 synthetic wavefront oracle bundle은 일반 traced bundle과 구별하고, 별도 명시적 검증 생략 옵션이나 타입으로 처리할 수 있다. 테스트 편의를 위해 production 경로 전체를 무검증으로 만들 필요는 없다.
- bundle provenance만 검증해도 충분한 곳과 trace 내부 배열까지 불변으로 유지해야 하는 곳을 구분해야 한다. boundary의 mutable-array 위험은 소프트웨어 리뷰와 함께 검토해야 한다.

## 직접 확인한 수식과 계산 절차

아래는 **읽고 유도했을 때 새 수식 결함을 발견하지 않은 항목**이다. 모든 입력과 수치 범위에서 독립적으로 검증했다는 뜻은 아니다.

### 실광선과 OPL

`geometry.py`의 sag, gradient, Hessian은 conic + A2…A12 표현과 일관된다. uniform medium 구간의 `OPL += n*distance`는 실광선 길이를 사용한다. 무한 object의 launch eikonal `n d·r`와 유한 object point의 공통 초기 위상은 서로 다른 광원 조건에 맞다. Snell은 normal 방향을 incidence 반대쪽으로 정렬하고 실수 transmitted branch를 구한다. 전역적으로 모든 asphere 교점과 aiming branch가 올바르게 선택되는지는 별도 기하 검증 대상이며, 렌즈 리뷰어가 clipping/경계 반례를 발견했다.

### 유한 결상 회절의 boundary → Fourier 좌표

코드의 chief-normal boundary 좌표가 `(u,v)`, 기준점까지 축상 거리가 `L`이면

\[
R=\sqrt{L^2+u^2+v^2},\quad (\chi,\eta)=(u/R,v/R),\quad
\gamma=\sqrt{1-\chi^2-\eta^2}=L/R.
\]

\[
\frac{\partial(u,v)}{\partial(\chi,\eta)}=\frac{L^2}{\gamma^4},\qquad
B=\frac{A L}{\gamma^2}.
\]

후자는 Rayleigh–Sommerfeld 원거리 kernel의 `L/R²`와 면적 Jacobian을 결합한 결과다. `W=S−S_chief+n(R−L)`의 OPD 부호와 `R−L=(u²+v²)/(R+L)`의 안정적인 계산도 서로 일관된다. Fourier 영상 좌표가 `(λ_vac/n)×주파수`가 되는 단위 변환은 적절하다.

그러나 image 주변 전개가 들어가므로 이것은 임의의 큰 defocus/NA/영상 창에서 정확한 Maxwell 혹은 exact RS 결과가 아니다. 제공되는 direct RS 비교는 같은 재구성 boundary를 공유하기 때문에 ray/OPD 생성 자체의 독립 검증이 아니다.

### RS kernel

`direct_rs`는

\[
K=\frac{L}{2\pi R^3}(1-ikR)e^{ikR},\quad k=2\pi n/\lambda_{vac}
\]

를 사용한다. 이 식은 scalar RS-I의 exact normal-derivative kernel과 부합한다. `1/(kR)` 항을 남겨 둔 점도 확인했다. 정상적인 반공간 경계조건이라는 가정은 필요하며 벡터 편광 전파를 의미하지 않는다.

### 원거리 광량의 측도·전기장 단위

`l,m`이 direction cosines이고 `γ=sqrt(1−l²−m²)`이면

\[
d\Omega=dl\,dm/\gamma.
\]

scalar RS far field의 intensity per solid angle에는 `γ²|F|²`, power cell에는 `γ|F|² dl dm`이 들어간다. `afocal_fft`의 두 값은 이 차이를 유지한다. 물리 전기장 power 변환 `n/(2Z0 λ_medium²)`와 mm²→m²의 `1e−6` 인자도 차원이 맞는다. 이 검토는 PH-01의 **높은 angular frequency 적분 해상도 문제**를 해소하지 않는다.

### 광량 전달과 source amplitude

- `uniform_scene_ri`의 `L/n²` invariance 및 `E=∫L cosθ dΩ=∫L dl dm` 형태에서는 별도의 cos⁴나 source cosine을 또 곱하지 않는 것이 맞다.
- Fresnel s/p power 계수와 수면마다 basis를 운반하는 unpolarized Jones-power 계산을 확인했다. bulk absorption은 `exp(-μs)`이고 amplitude에서는 그 제곱근을 곱한다.
- finite source의 `Q_triangle=L·fill·∫cosθ_source dΩ·mean(τ)` 및 `E²=2Z0 I_normal/(n cosθ_boundary)`는 power packet의 제곱근을 그대로 field로 쓰는 오류를 피한다.
- source cell area와 fill을 한 번씩 적용하고 normalized kernel과 collected power를 별도로 보존하는 구조는 적절하다.

다만 triangle 방향 경계는 great-circle chord, transmission은 vertex mean, 물리 boundary amplitude는 triangle별 상수다. 고NA 편광의 pupil별 vector coherence를 Jones **power**만으로 보존할 수는 없으며, 현재는 scalar 모델의 제한이다. 임의의 빠른 렌즈에서 정확한 물리 PSF임을 주장하려면 추가 검증이 필요하다.

## 결함과 구별해야 하는 기존 모델 제한

1. **Zemax 동일 구현이라는 보증 없음.** 현재 profile은 독립 scalar FFT이고, 초점면은 chief-normal virtual plane이며 lens-only pupil illumination policy도 별도 정의다. Ansys 문서에 명시된 FFT의 적용 가정, exit-pupil sampling 및 image-surface 조건과 상호 검증해야 한다. 공개 scalar 식을 사용한다는 것만으로 상용툴과 동일한 지표가 되지는 않는다.
2. **같은 pupil에서 phase=0으로 만든 reference**는 적절한 same-pupil diffraction reference다. 균일 원형 aperture의 해석적 diffraction-limited curve와 항상 같은 것은 아니다.
3. **조건부 정확도.** `conditional`은 실제 MTF 오차가 얼마 이하라는 보증이 아니다. phase step, ray triangulation, aperture raster, angular/spatial output quadrature, local propagation approximation, raw power 각각의 수렴 검증이 필요하다.
4. **비축 영상면.** source-wave의 global-z detector 투영은 명시적 local affine opt-in이며 exact tilted-plane diffraction이 아니다.
5. **단일 경로.** ghost, scatter, vector diffraction, 다중 branch 및 일반 caustic 전파는 이 검토가 승인한 범위가 아니다.

## 출처와 증거 구분

아래 공식/저자 자료는 해당 방법과 적용 조건의 근거로 사용했다. 재현 버그의 직접 증거는 위 스크립트와 동결된 코드이다. 문서를 인용했다는 이유만으로 코드가 검증된 것은 아니다. 조회일: 2026-10-18.

- Ansys, **Methods for analyzing MTF in OpticStudio** — focal/afocal frequency 단위, Fraunhofer FFT 방법, scalar 적용·pupil distortion·sampling 조건.
  <https://optics.ansys.com/hc/en-us/articles/42661772706835-Methods-for-analyzing-MTF-in-OpticStudio>
- Ansys, **What is a Point Spread Function?** — pupil amplitude/OPD→FFT PSF 과정 및 chief ray의 image surface 입사 조건.
  <https://optics.ansys.com/hc/en-us/articles/42661723066515-What-is-a-Point-Spread-Function>
- Ansys OpticStudio help, **FFT MTF** — exit-pupil cosine-space sampling 및 해석 설정.
  <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html>
- Ansys OpticStudio help, **Relative Illumination** — field 중 최대 image-surface illumination으로 정규화하는 정의.
  <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Relative_Illumination.html>
- S. Konijnenberg, A. J. L. Adam, H. P. Urbach, TU Delft Open 저자 교재, **6.4: Rayleigh–Sommerfeld Diffraction Integral** — Helmholtz 기반 scalar 전파, obliquity factor, 생략 가능한 작은 항 및 angular-spectrum과의 관계. 아래는 저자 교재의 LibreTexts 공개본이다.
  <https://phys.libretexts.org/Bookshelves/Optics/BSc_Optics_%28Konijnenberg_Adam_and_Urbach%29/06%3A_Scalar_diffraction_optics/6.05%3A_Rayleigh-Sommerfeld_Diffraction_Integral>

## 실행

두 repro는 인접한 `../input/optical_engine`의 소스를 명시적으로 import한다. 시스템에 설치된 다른 버전이 테스트를 가로채지 않게 했다.

```bash
python physics/repro_stale_bundle.py
python physics/repro_angular_otf_alias.py
```

원본 엔진 수정은 하지 않았다. 이번 독립 반례를 먼저 regression test로 추가하고, 결함 수정 뒤 실제 처방 기반 수렴 검증을 이어가는 것이 적절하다.
