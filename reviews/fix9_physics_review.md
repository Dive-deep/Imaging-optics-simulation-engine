# 9개 결함 수정본 — 독립 광학 물리 관점 리뷰

검토일: 2026-10-18. 대상: 작업 루트 `src/sequential_optics/`. 비교 기준은 변경하지 않은 `audit_round2/input/optical_engine/`이다. 별도 AI 에이전트가 광학 물리학자 역할로 수행한 검토이며 외부 사람 전문가 인증은 아니다. 엔진 코드는 이 검토자가 변경하지 않았다.

## 결론

**검토 범위 R1/R2의 원인에 대응하는 수정은 타당하며, 미해상 angular MTF의 거짓 정상 승인과 변경 전 광선의 처방 재사용을 차단한다.** 기존 scalar 전파·진폭·광량 식에서 이 수정으로 새 물리 불일치가 생긴 증거는 찾지 못했다. 이는 전체 모델의 임의 처방 정확도, Zemax 상관, 수치 수렴을 보증하는 판정은 아니다.

R1/R2 회귀 테스트 최종 18개와 기존 독립 wave/source-wave 물리 대조 8개를 작업 `src`로 import 경로를 강제하여 실행했고 모두 통과했다. 아래 별도의 직접 대조도 수행했다. 전체 프로젝트 테스트 수와 결과는 루트 검증 보고서가 관리한다.

## 1. R1: angular OTF의 출력 적분 위상 검사

원래 결함은 경계 파면 위상 sampling을 충분하게 해도, 출력 적분 `sum(p_k exp(-i 2 pi f alpha_k))`의 각도 표본이 부족하면 alias가 발생한다는 것이다. 수정은 이 둘을 독립 검사한다.

- 실제 direction-cosine cell의 angular span을 `atan2` 좌표에서 계산한다. 각 cell의 양 끝과 cross-coordinate의 최소 절댓값을 포함하고, propagating disk 경계에서는 보수적인 수평 방향 한계를 사용한다. 각도와 direction cosine을 단순 동치로 두지 않는다.
- 주파수별 `2 pi |f| Delta alpha`가 설정 한계를 넘는 cell의 전력 비율을 계산한다. 한계와 허용 미해상 전력은 settings로 분리했다.
- 미해상 저수준 호출은 명시적으로 거부한다. 진단 경로는 해당 주파수에 NaN을 기록하며, 상위 lens 및 physical-source wave는 `underresolved`로 전달한다.
- 실제 pupil과 같은-pupil zero-phase reference에 동일한 검사를 적용한다. 고주파를 임의 hard cutoff로 0으로 만들지 않으며, 문제가 있는 cell을 삭제하거나 남은 전력을 다시 정규화하지 않는다. DC는 1로 보존한다.

이 검사는 **채택을 위한 보수적 sampling guard**다. 출력 함수 자체의 수렴이나 midpoint quadrature 오차 상한을 증명하지 않는다. 따라서 `convergence_verified=False`와 conditional 상태의 한계를 유지한 것이 적절하다. 수학적으로 모든 가능한 함수에 대한 절대적 필요충분조건이라고 해석해서는 안 된다.

회귀 검사는 기존 Gaussian alias를 거부하고 padding 16/24 결과가 1e-8 절대 오차 이내로 일치하는 것을 확인한다. 또한 high-angle fixture에서 direction-cosine Nyquist만으로 놓치는 각도 해상도 부족을 잡는다.

독립 직접 대조에서는 고정 난수 seed 401로 생성한 100개 후보 중 propagating 중심을 갖는 cell 각각을 301×301 내부 격자로 조사했다. X/Y 실제 angular range와 보수 bound 차이의 최대값은 부동소수점 수준인 약 1.1e-15 rad였다. 이 샘플 대조는 bound 구현의 지지 증거이며 모든 입력에 대한 증명은 아니다.

## 2. R2: trace provenance와 boundary 재사용

새 trace provenance는 소비 시점이 아니라 `trace_bundle`이 생성한 시점에 보존한다. 처방 identity와 기록 내용 identity를 분리한 것이 핵심이다.

- 처방 identity는 전체 면, STOP index/radius, 물점·공액, 파장, 온도, 실제 사용 굴절률, ray solver 설정을 포함한다.
- 내용 identity는 최종 ray만이 아니라 중간 위치·방향·OPL·power, 성공 상태, chief, field, wavelength, pupil 좌표와 적분 가중치를 포함한다.
- wave boundary 구성 전에 provenance를 확인하므로 오래된 광선과 새 처방을 함께 hash해서 다시 정상으로 인정하는 이전 경로가 사라졌다.
- IMA만 이동할 때는 광선 추적이 마지막 물리면에서 끝난다는 계약에 따라 동일 ray를 허용한다. FFT grid/padding/보고 주파수도 ray 입력과 구별한다.
- 미사용 RGB·CTE·흡수 데이터를 요구하지 않으며, 동일 유효 굴절률의 투과율 정보 갱신은 지오메트리 ray를 불필요하게 무효화하지 않는다.
- 직접 만든 analytic bundle은 명시적 opt-in으로만 lens-wave 경로에 진입하며 unverified 표시를 유지한다. 이 opt-in으로 mismatch가 있는 engine trace를 통과시킬 수 없다. Physical-source 경로는 검증된 engine trace를 요구한다.

추가된 BoundaryState는 bytes-backed 배열과 frozen dataclass로 방어적 snapshot을 소유한다. 별도 content seal은 OPL, 좌표, 굴절률, reference distance 등의 교체도 탐지한다. 검토자가 시도한 cached OPL의 직접 변경은 ValueError로 거부되었으며, dataclasses.replace를 이용한 값 교체도 회귀 테스트에서 거부된다. 이는 실수 방지용 무결성 계약이지 악의적 호출자 인증 장치는 아니다.

광원 투과율 갱신도 별도로 확인했다. n=1로 index matched된 1 mm plate에서 동일 ray bundle을 재사용하며 내부 흡수만 0에서 ln(2)/mm로 바꿨다. Physical source 수집 전력 비율은 **0.4999999968288031**이었다. 즉 ray identity에 필요 없는 loss 정보는 재추적을 강제하지 않으면서 실제 광량은 현재 물성으로 다시 계산한다. 첫 면의 Snell/OPL 검사만으로 뒤쪽 면 변경을 추정하던 이전 약점도 provenance가 보완한다.

## 3. 전파식과 전력 측도의 일관성

수정본의 핵심 식을 직접 미분·단위 대조했다.

1. 공통 physical boundary까지의 eikonal은 실제 전파 길이에 n을 곱해 더한다. Focal OPD는 `W = S - S_chief + n(R-L)`이고 위상은 `2 pi W / lambda_vacuum`이다. 굴절률을 중복 적용하지 않는다.
2. `chi=u/R`, `eta=v/R`, `gamma=L/R`에서 면적 Jacobian은 `L^2/gamma^4`다. RS local obliquity와 결합하면 Fourier amplitude가 `B=A L/gamma^2`가 된다. Focus 이동에서 physical A를 유지하고 B를 다시 구성하는 절차가 일관된다.
3. Afocal 전력 cell은 `gamma |F|^2 dl dm`에 비례하며 단위 입체각 density는 `gamma^2 |F|^2`에 비례한다. `dOmega=dl dm/gamma`를 곱하면 두 표현이 일치한다.
4. Physical-source triangle의 전력은 source projected solid angle, radiance, source 면적, 투과율로 구한다. 전기장 진폭은 실제 boundary 면적과 crossing cosine, n/(2 Z0) flux 관계로 변환한다. mm²→m²의 1e-6을 포함한다.
5. PSF 정규화 전에 ray/triangle power, raster boundary flux, propagating spectral power를 유지한다. 최종 Q×kernel 합만 맞춰 raw power 문제를 숨기는 방식이 아니다.

재실행한 기존 독립 물리 대조는 원형 aperture의 on-axis complex RS 해, 기울어진 단일 spectral mode의 raw flux cosine, 출력 굴절률에 따른 focal 축 배율, 실제 굴절 ray의 eikonal-gradient 관계, absolute n에 따른 전기장 진폭과 보존 전력, 흡수의 단일 적용, global detector 등록, 물점 거리 변경 거부를 포함한다.

## 4. 핵심 원리 문서 검토

`docs/core_principles_ko.md`의 식 (1)–(44)를 코드의 정의와 대조했다. 특히 다음을 확인했다.

- sag 미분의 A2 기여, reduced-ray ABCD 및 EFL/BFL 정의가 일관된다.
- 실제 광선과 paraxial 정의량, geometric T/S와 FFT 좌표축, global detector와 chief-normal wave plane을 구분한다.
- OPD reference sphere, Fourier chart Jacobian, medium wavelength, complex OTF와 intensity 합성의 순서가 명시되어 있다.
- RI의 radiance/n² 불변성과 direction-cosine 측도, Lambertian projected solid angle, 전력에서 경계장으로의 변환이 구별된다.
- 열팽창 asphere 계수의 `s^(1-2m)` 변환과 공동 RMS 초점 식은 해당 모델에서 직접 도출 가능하다.

원문 출처 중 Lucke 논문의 식 (5), (7), (8)과 MIT OTF 강의 자료를 직접 열어 원자료와 문서가 인용하는 범위를 확인했다. 해당 논문은 direction cosine Fourier 관계와 obliquity의 근거이며, 자체 엔진의 ray reconstruction이나 모든 구현 세부가 그 논문에 있다는 주장은 하지 않는다.

- Robert L. Lucke, *Rayleigh-Sommerfeld Fraunhofer Diffraction*: https://arxiv.org/pdf/physics/0604229
- MIT OCW 2.71, Lecture 22: https://ocw.mit.edu/courses/2-71-optics-spring-2009/bebe78c6a2ae2603a0d1a3ea85e0b4ca_MIT2_71S09_lec22.pdf

## 5. 남아 있는 모델 범위

Angular guard 통과는 출력 적분 수렴 인증이 아니다. 원래 경계장 grid의 undersampling, 제한된 propagating window, pupil support rasterization, source macro-cell 적분은 독립적으로 수렴시켜야 한다. 같은 scalar 경계장을 공유하는 direct RS 비교도 ray/OPD 재구성까지 독립 검증하지는 않는다.

실제 출사 경계장에 대한 scalar/단일 branch 모델, 초점 근방 국소 Fourier 전개, 비축 global detector의 선택적 국소 투영이 현재의 모델 범위다. Vector diffraction, 임의 큰 NA/비축 처방의 정확성, Zemax 내부 pupil 규약 일치와는 구별한다. 이번 R1/R2 수정은 이 범위를 넓히는 작업이 아니라 기존 범위에서 잘못된 수치 상태와 재사용을 차단하는 작업이다.
