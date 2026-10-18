# 3매 Sequential 광학 계산 엔진 — 이론·수식·계산 절차 명세

문서 버전: **1.0-draft-for-implementation** · 작성일: **2026-10-18 (Asia/Seoul)**

이 문서는 사용자의 광학 처방을 받아 MTF, spot, RI, 수차, 초점거리, 영상·grid, 균일 온도 성능을 계산하는 핵심 엔진의 구현 자료다. 기존 수식 자료 rev0.4, 8개 전문 감사와 광학 물리학자·카메라 설계자·소프트웨어 엔지니어의 독립 검토 및 교차토론을 반영했다. 기존 자료의 발췌 모음이 아니라 **계산 순서와 모듈 사이의 수학적 연결을 명시한 독립 기술 명세**다.

**검증 상태:** 본 문서의 이론 정리와 식의 제한적 검산은 실제 3매 엔진의 구현·검증 또는 Zemax 출력 대조를 뜻하지 않는다. 특히 공개된 FFT 이론으로 구현하는 독립 profile과 Zemax 내부의 모든 cosine-space 처리·옵션을 동일시하지 않는다. 계산 가능한 독립 규약은 구체적으로 지정하고, 제품 수치 호환이 미검증인 부분과 아직 선택이 필요한 물리 모델은 뒤의 쟁점 목록에 남긴다.



## 문서 목차

1. [범위·입력·전체 계산 순서](#part-1)
2. [재료·온도별 형상 전처리](#part-2)
3. [실광선 커널·초점거리·기하 수차](#part-3)
4. [OPD·FFT PSF/MTF·회절 기준](#part-4)
5. [전력·RI·원거리 분포·영상·grid](#part-5)
6. [균일 온도 성능과 공통 재초점](#part-6)
7. [결과·검증·미결 쟁점과 피드백 반영](#part-7)
8. [통합 참고문헌](#part-8)

처음 구현할 때는 1→2→3→4 순서로 한 field의 G 결과를 연결한다. 광량·영상은 5, 온도 스윕은 6을 적용하고, 7의 상태·검증 요구를 각 단계에서 함께 구현한다.

<a id="part-1"></a>

## 1. 범위·입력·전체 계산 순서

### A. 읽는 방법과 규범

문서의 문장을 다음 네 종류로 구분한다.

| 표기 | 의미 |
|---|---|
| **이론** | 명시한 가정 안에서 성립하는 수학·물리 관계. 출처 또는 직접 유도를 표시한다. |
| **본 구현 규약** | 여러 가능한 convention 중 재현 가능성을 위해 선택한 기준. Zemax의 모든 설정 기본값을 뜻하지 않는다. |
| **제품 정의** | 공식 OpticStudio 문서에서 확인한 분석의 의미·단위·동작. 특정 build의 수치 일치는 별도 대조가 필요하다. |
| **OPEN / 조건부** | 알고리즘 동등성, 적용성 또는 필요한 입력이 아직 닫히지 않은 항목. 값을 임의로 채우거나 다른 분석으로 조용히 대체하지 않는다. |

식의 파생은 본문에서 정의한 단위·좌표를 사용한다. 각 장의 별도 기호는 그 장의 정의를 따른다. 각 모듈은 입력, 출력, 가정, 실패 상태를 함께 갖는다. 뒤의 출처 목록은 원문 전문을 복사한 것이 아니라 구현에 필요한 사실과 유도의 근거를 연결한 것이다.

**범위:** 3개의 동축 회전대칭 렌즈와 6개 굴절면, 균질·등방성 매질, 순차 직접 투과 경로. Standard 또는 Even Asphere A₂…A₁₂. STOP 전용 무굴절 dummy 면을 추가해도 렌즈 수가 늘어나는 것은 아니다. 기본 OBJ/IMA는 평면이다. 임의 decenter/tilt, 곡면 센서, ghost·scatter·다중 반사, 제조 공차, 열전도·열구배·응력·GRIN·복굴절·mount 변형은 기본 지원에 포함하지 않는다. 해당 기능이 없다는 이유로 현재 범위의 계산을 중단하지 않는다.

**실광선의 의미:** 실제 sag의 교점과 벡터 Snell 법칙으로 유한 aperture 광선을 추적한다. 광선 경로를 근축 경로로 대체하지 않는다. 다만 EFL/BFL은 1차 광학량, field curvature는 chief 주변의 미분적 광선량, FFT PSF는 명시한 scalar diffraction 근사다. 모든 출력이 무근사 Maxwell 해라는 의미는 아니다.

### B. 전체 계산 의존 관계

1. 처방·분석 요청을 읽고 입력 convention, 단위, 좌표, 재료 기준 상태를 정규화한다.
2. 기준 온도 재료와 형상을 검증하고 aperture 변환에 필요한 1차 광학량을 구한 뒤, 입력 EPD/F/#를 물리 STOP으로 resolve한다. 온도별 재료 n·흡수와 렌즈 형상은 원 기준 상태에서 생성한다.
3. 전체 clear 영역에서 형상·면 순서·두께·교차·재료 유효 범위를 검사한다.
4. 필요한 first-order 계산을 수행해 초기 aiming과 EFL/BFL을 얻는다. 최종 실광선은 별도로 계산한다.
5. 각 field·파장에 대해 chief와 pupil 광선의 STOP aiming을 풀고 면별 실제 경로·OPL·전력을 저장한다.
6. 동일 광선에서 spot, ray fan, chief mapping, differential T/S curvature 등 기하 출력을 계산한다.
7. 같은 광선의 공통 입사 위상과 출력 경계장을 연결하여 reference OPD·복소 pupil을 만든다. 선택한 FFT profile로 PSF·OTF·MTF·동일 pupil의 회절 기준을 계산한다.
8. 일반 RI는 균일 장면에 대한 별도 radiance 적분, finite-source far-field는 발광면의 면적·방향 적분으로 계산한다. 두 출력의 분모와 측도를 구분한다.
9. Image/grid는 source별 수집량 Q, 해당 조명에서의 PSF, chief mapping을 결합한다. 각 field를 개별 재초점한 영상을 공통 센서 영상으로 저장하지 않는다.
10. 고정 IMA/정해진 angular 상태와 별도 공통 refocus 상태를 반복 평가한다. 기준 온도 대비 mapping·광량·성능 변화를 함께 저장한다.
11. 각 metric의 수치 수렴·물리 적용성·입력 상태를 평가하고 원 데이터와 함께 반환한다.

```text
evaluate(design, materials, profile, request):
    base = resolve_units_coordinates_and_reference_state(design, materials, profile)
    validate_geometry_and_material_ranges(base)
    first_order_at_T0 = paraxial_properties(base, profile.primary_wavelength)
    base.physical_stop = resolve_aperture_once_at_reference_temperature(base, first_order_at_T0)
    for T in request.temperatures:
        state = build_from_reference(base, T)  # 이전 온도 상태를 다시 팽창시키지 않음
        validate_geometry_and_material_ranges(state)
        first_order = paraxial_properties(state, requested_wavelengths)
        for field, wavelength in request.field_wavelength_pairs:
            metric_rays = aim_and_trace(state, field, wavelength, profile.lens_metric)
            geometry = geometric_metrics(metric_rays, state.image)
            wave = wave_metrics(metric_rays, state, profile.wave)
            save_raw_and_derived(geometry, wave)
        if requested: uniform_scene_RI(state, profile.radiometry)
        if requested: finite_source_transport_and_farfield(state, request.source)
        if requested: image_and_grid_from_matched_kernels(state, request.image)
        if requested: optimize_one_common_compensator_and_reevaluate(state)
        record_metric_states_and_convergence()
    return ResultBundle
```

위 의사코드는 기능의 의존 순서를 나타내며 구현 완료된 API를 뜻하지 않는다.

### C. 입력 명세와 내부 단위

입력은 `DesignSpec`, `MaterialSet`, `AnalysisProfile`, `RunRequest`로 나누고 결과는 `ResultBundle`로 반환한다. 이는 같은 값을 여러 위치에서 독립 편집하지 않게 하기 위한 구분이다. 기존 P001–P111은 추적 ID로 유지할 수 있으며 전부 설계마다 묻는 필수값은 아니다.

| 입력군 | 필요한 내용 | 구현에서 지킬 조건 |
|---|---|---|
| 6면 처방 | surface type, R 또는 c, conic k, A₂…A₁₂, 앞면 center thickness, 뒷면 다음 vertex까지 gap | Standard의 A항은 0. 마지막 면 뒤 IMA 위치는 별도 의미. 중심 두께와 gap을 혼동하지 않는다. |
| 구경·기구 | 각 면 active/clear diameter, mechanical diameter, 지정 무효부 높이/형상 | Diameter와 semi-diameter 명시. 광학 계산에서 비유효 영역은 기본 차단. 미정 rim을 광학 굴절면으로 외삽하지 않는다. |
| 매질 배치 | 각 구간의 재료 ID, 외부·내부 공기/진공 | 렌즈 front 뒤는 렌즈 재료, rear 뒤는 air 등 실제 구간으로 resolve. 물성의 absolute/relative 기준 보존. |
| OBJ·field | finite/infinite, OBJ 위치, object-height 또는 angle 좌표, 실제 2D field 값 | 정규화 field이면 역변환용 최대 field를 함께 지정. OBJ 위치와 광원면 위치의 관계를 명시. |
| IMA·출력 모드 | 위치, focal/afocal, chief-normal 분석면과 실제 IMA의 구분 | Afocal mode는 고유 power C=0이라는 조건이 아니다. 물리 IMA와 마지막 gap을 중복 강제하지 않는다. |
| STOP·system aperture | STOP surface ID, EPD/F-number/physical stop 등 한 입력 방식과 값 | STOP 위치는 별도 입력. F/# 종류를 명시하고 서로 모순되는 독립 입력을 거절. |
| 파장 | 실제 R/G/B 파장 수치, 입력 매질 기준, G 식별 | 내부 λ₀는 vacuum. 상대 RGB weights는 다파장 합성에 필요. G 단색 렌즈 MTF에는 광원 bitmap이 필요 없다. |
| 재료 DB | 실제 grade의 n(λ,T), 기준 T/P와 파장 convention, CTE·정의, internal transmission의 기준 두께/범위 | N-BK7/crown/flint/PMMA/PC 이름만으로 임의 계수를 채우지 않는다. 제공된 실제 값과 model ID를 사용한다. |
| RI·광원 | 12.8×3.2 mm, 중심 축상, 균일 Lambertian, fill=0.64 연속 근사, emitting side | 실제 pixel은 metadata. 절대 radiance가 없으면 상대량은 계산 가능. 절대 W/sr에는 절대 방출량이 필요. |
| 영상 입력 | sample image, source 좌표에 대한 크기·방향, 채널의 선형 방출량 의미 | 표시용 sRGB code를 물리 intensity로 바로 취급하지 않는다. 채널과 파장 대응은 선언한다. |
| 온도 요청 | 형상 기준 T₀와 평가 T, 위치 anchor, gap·STOP·OBJ·IMA 정책, 허용 refocus 자유도 | n fit의 기준 T와 형상 T₀는 별도. 파장은 기본 vacuum λ 고정. |
| 분석 규약 | pupil/illumination, OPD reference, PSF 방법, 좌표·T/S, 수치 표본, 정규화·reporting domain | 버전으로 고정한다. 실제 Zemax 초기값이라고 임의로 부르지 않는다. |

#### C.1 공통 단위와 기호

| 기호/값 | 내부 단위·의미 |
|---|---|
| x,y,z,R,t,OPL,W | mm. z는 기본 광축 진행 방향. W는 OPD 길이. |
| λ₀ | 진공 파장 mm. 입력 µm는 10⁻³을 곱해 mm로 변환. |
| n | 절대 굴절률. 매질 내 파장 λₘ=λ₀/n. |
| d | 단위 광선 방향. Optical direction n·d와 구분. |
| 각도 | rad. 출력 deg/mrad는 마지막 변환. |
| f | Focal: cycles/mm=lp/mm. Afocal: cycles/rad 또는 명시한 angular unit. |
| A₂m | mm^(1−2m). 수치 단위를 바꿀 때 계수도 변환. |
| Φ, I, E, L | W, W/sr, W/mm², W/(mm²·sr)를 기본 조합으로 사용할 수 있다. SI 변환을 명시. |
| τ, η | 각각 전력 투과율과 명시된 분모의 효율. 혼동하지 않음. |
| h(x;q) | Source q에 대한 intensity PSF, 영상면 적분 1이면 mm⁻². |
| hΩ(ω;q) | 입체각 적분 1인 angular kernel, sr⁻¹. 각도 평면 밀도와 구분. |
| Q(q) | Source 면적당 지정 출력으로 전달되는 전력 W/mm². |

의미가 다른 단위가 우연히 숫자상 비슷해도 같은 배열로 취급하지 않는다. 예를 들어 angular histogram의 bin당 W, W/sr, 정규화 intensity, angular MTF는 서로 다른 출력이다.

#### C.2 기본 조명·분석 profile

- `lens_metric`: 사용자 확정 G 렌즈 단독 FFT MTF. Source pixel/fill·bitmap의 전달함수를 곱하지 않는다. 균일 진폭을 어느 Fourier pupil 측도에서 정의하는지는 회절 장의 profile에 고정한다.
- `source_transport`: 공간적으로 비간섭인 연속 Lambertian 발광면. 각 물점·한 파장에서는 coherent pupil 전파, 물점·파장 사이에는 intensity 합산.
- `uniform_scene_RI`: 균일 평면 Lambertian 장면의 국소 영상 조도. 유한 사각 광원의 edge blur를 이 분석에 섞지 않는다.
- `thermal_fixed_geometry_policy`: 기준 상태에서 물리 STOP을 정한 후 온도별 위치·팽창 규칙 적용. 매 온도 목표 F/#를 재강제하는 경우는 별도 제어 시나리오로 기록한다.
- `focus_illumination_policy=hold_boundary_field_amplitude`: IMA만 이동하는 공통 초점 스윕에서는 같은 canonical 경계평면과 경계장 진폭을 유지한다. 변경된 reference에 맞춰 Fourier pupil을 다시 구성하되 매번 균일 B로 초기화하지 않는다. 대안 규약은 회절 장에서 별도 이름으로 정의한다.

다음 장부터 이 입력을 실제 수식과 계산 순서로 연결한다.

---

<a id="part-2"></a>

## 2. 재료·온도별 형상 전처리

이 장은 제공된 실제 재료 물성으로 각 온도의 광학 처방을 만든 다음, 같은 real-ray·OPD·FFT 엔진을 다시 실행하는 절차를 정의한다. 열전도·시간 응답·응력·mount 변형 해석은 수행하지 않는다. 따라서 여기서 온도란 각 렌즈와 지정한 주변 매질이 균일한 평형 온도에 도달한 상태다. 아래 위치 규칙은 그 제한 모델의 기하학적 경계조건이며 구조 FEA를 요구하지 않는다.

### 1. 재료 입력을 계산 가능한 함수로 변환

#### 1.1 공통 내부 규약

광선 엔진의 물성 인터페이스는 다음 세 함수를 제공한다.

$$
n_j^{\mathrm{abs}}(\lambda_v,T),\qquad
\mu_j(\lambda_v,T),\qquad
s_j(T;T_0).
$$

- $\lambda_v$: 진공 파장. 공통 엔진에서는 mm로 저장하고 제조사 식을 평가할 때만 해당 식의 단위로 바꾼다.
- $n_j^{\mathrm{abs}}$: 진공 기준 위상 굴절률. Snell과 OPL에 동일한 값을 사용한다. 군굴절률을 대신 사용하지 않는다.
- $\mu_j$: 내부 **전력**의 흡수 계수, 단위 mm⁻¹. 
- $s_j$: 기준 형상 대비 등방성 선형 치수의 비율. 재료별로 계산한다.

각 함수는 값 외에 `model_id`, 원자료 식별자/버전, 파장·온도 유효범위, reference 조건, 보간/외삽 정책을 반환해야 한다. 범위를 벗어난 값은 정상값과 구분한다. 현재 목표는 물성이 제공된 경우의 계산이며, 여기서 새 제조사 조사나 물성 확보를 선행조건으로 다시 요구하지 않는다.

N-BK7, 별도 crown/flint, PMMA, PC는 재료 선택지다. 실제 데이터는 `manufacturer + grade + dataset_version`을 키로 연결한다. Crown/flint 또는 PMMA/PC라는 분류명만으로 굴절률·CTE·투과율을 고정하지 않는다. 서로 다른 제품의 n, CTE, 투과율을 임의 조합한 가상 재료를 실제 grade처럼 저장하지 않는다.

#### 1.2 RGB의 n(T) 표만 있는 경우

계산 파장이 RGB 세 값뿐이면, 각 진공 파장의 $n^{\mathrm{abs}}(T)$ 표로 충분하다. 별도 Sellmeier fit은 필수 조건이 아니다. 온도 표의 사이를 선형 보간하는 프로파일을 택했다면 다음 식을 사용하고, 구간 밖은 외삽하지 않는 것이 기본이다.

$$
n(T)=(1-u)n(T_a)+u n(T_b),\qquad
u=\frac{T-T_a}{T_b-T_a},\quad T_a\le T\le T_b.
$$

이는 보간 모델이며 재료의 정확한 온도 법칙이라는 주장은 아니다. 표가 air-relative이면 **각 표본의 기준 공기 조건을 먼저 해석**해 공통 absolute 값으로 변환한다. 단순히 숫자에 `absolute`라는 label만 바꾸면 안 된다.

분산식 계수가 제공된 경우에는 해당 제조사/카탈로그의 식을 사용한다. 예를 들어 3항 Sellmeier 형식은 다음과 같다. $B_i,C_i$, 파장 단위와 공기/진공 파장 기준은 데이터와 함께 제공되어야 한다. [TH01](#ref-TH01)

$$
n^2(\lambda)-1=\sum_{i=1}^{3}\frac{B_i\lambda^2}{\lambda^2-C_i}.
$$

Abbe 수는 다음의 요약량이다. $n_d,V_d$만으로 임의 RGB 파장의 분산곡선을 유일하게 복원할 수 없다. 표 또는 전체 분산식이 있을 때 Abbe 수는 검증용 metadata로 보존한다. RGB의 G가 d선 또는 e선과 같다고 가정하지 않는다. d/F/C선의 굴절률과 파장·공기 기준은 catalog 원래 정의를 사용한다. 공통 내부 n_abs로 바꾼 값을 그대로 넣은 Abbe 수가 catalog의 reference-air Abbe 수와 동일하다고 강제하지 않는다. [TH01](#ref-TH01)

$$
V_d=\frac{n_d-1}{n_F-n_C}.
$$

#### 1.3 Absolute/relative index와 파장 기준의 변환

같은 물리적 파동을 비교할 때 다음 관계가 성립한다.

$$
n_{j,\mathrm{rel}}=\frac{n_j^{\mathrm{abs}}}{n_{\mathrm{air}}^{\mathrm{abs}}},\qquad
\lambda_{\mathrm{air}}=\frac{\lambda_v}{n_{\mathrm{air}}^{\mathrm{abs}}},\qquad
\frac{n_j^{\mathrm{abs}}}{\lambda_v}=
\frac{n_{j,\mathrm{rel}}}{\lambda_{\mathrm{air}}}.
$$

따라서 모든 굴절률과 파장을 일관되게 변환하면 Snell의 굴절률 비와 전파 위상은 같다. 반대로 `n_rel`과 `λ_vac`를 섞으면 공기/유리의 위상·주파수 환산이 일관되지 않는다. 직접 제공된 absolute 데이터의 기본 온도 실험에서는 **$\lambda_v$를 고정**하여 같은 광원 색을 유지한다.

OpticStudio의 입력 파장은 system T/P의 공기 중 파장이고, 그 상태의 공기 굴절률을 1로 정규화한다. 공통 엔진의 absolute convention과 다르므로 catalog adapter에서 변환한다. [S54](#ref-S54)[S55](#ref-S55)[S57](#ref-S57) 아래 두 thermal profile은 비슷한 계수를 사용하지만 **파장 인자의 규약이 서로 같다고 합치지 않는다**.

| 프로파일 | 분산·온도 처리 | 내부 엔진 출력 |
|---|---|---|
| `provided_absolute_n` | 제공한 진공 파장별 n(T)를 직접 사용 | absolute n |
| `manufacturer_schott_tie19` | 제조사 분산 기준을 해석한 뒤 TIE-19의 vacuum-λ 온도식을 사용 | absolute n |
| `zemax_agf_reference` | system-air λ → reference-air λ, catalog dispersion, catalog thermal 규약을 적용 | 계산된 absolute n; 비교용 relative n도 보존 |

공기의 압력은 재료의 catalog reference pressure와 운용 pressure를 구분한다. 이 과제에서는 균일 air 온도·압력 모델을 명시하면 충분하다. 밀폐 구조의 압력 변화나 재료의 압력 의존 굴절률을 자동 추정하지 않는다.

### 2. 굴절률의 온도 모델

#### 2.1 제공된 함수 또는 측정된 dn/dT

가장 직접적인 경로는 제공된 $n^{\mathrm{abs}}(\lambda_v,T)$를 평가하는 것이다. 일정한 absolute $\beta=dn/dT$의 선형 모델을 명시한 경우에만 다음 식을 사용한다.

$$
n^{\mathrm{abs}}(\lambda_v,T)=n^{\mathrm{abs}}(\lambda_v,T_c)
 +\beta^{\mathrm{abs}}(\lambda_v)(T-T_c).
$$

Relative $dn/dT$는 absolute $dn/dT$와 동일하지 않다. 공기의 온도 의존성도 포함되기 때문이다. CTE 또는 Abbe 수로 $dn/dT$를 도출하지 않는다. 이 선형식의 온도 범위는 원자료 또는 선택한 근사 프로파일의 범위다.

#### 2.2 SCHOTT 제조사 thermal profile

SCHOTT TIE-19 식에서 $\delta=T-T_c$, $\lambda_v$는 **μm**이고 $T_c$는 계수의 기준 온도다. 제공된 계수 $D_0,D_1,D_2,E_0,E_1,\lambda_{TK}$로 다음을 계산한다. [S29](#ref-S29)

$$
F(\lambda_v,\delta)=D_0\delta+D_1\delta^2+D_2\delta^3+
\frac{E_0\delta+E_1\delta^2}{\lambda_v^2-\lambda_{TK}^2},
$$

$$
\Delta n^{\mathrm{abs}}(\lambda_v,T)=
\frac{n_c^2-1}{2n_c}F(\lambda_v,T-T_c),\qquad
n^{\mathrm{abs}}(\lambda_v,T)=n_c^{\mathrm{abs}}+\Delta n^{\mathrm{abs}}.
$$

제조사 설명은 prefactor의 $n_c$로 기준 온도의 relative index를 사용하는 근사를 허용한다. 해당 convention을 `prefactor_index_basis`로 고정한다. TIE-19의 일반 계수 적용 범위는 −100…+140 °C, 0.3650…1.014 μm이며, 실제 실행은 **해당 grade 데이터의 더 좁은 범위가 있으면 그것을 따른다**. 이 계수를 PMMA/PC의 보편 법칙으로 적용하지 않는다. [S29](#ref-S29)

굴절률 기준 `index_basis`와 분산식의 파장 인수 기준 `wavelength_basis`는 독립 metadata다. 기준 분산 데이터가 reference-air index이고 인수도 reference-air wavelength인 경우, 선택한 vacuum-argument air model을 $a_c=n_{\mathrm{air}}^{\mathrm{abs}}(\lambda_v,T_c,P_c)$로 평가하고 $\lambda_c=\lambda_v/a_c$, $n_c^{\mathrm{abs}}=a_c f_{\mathrm{disp}}(\lambda_c)$로 변환한다. 반면 relative index가 vacuum wavelength의 함수로 주어졌다면 $n_c^{\mathrm{abs}}=a_c f_{\mathrm{disp}}(\lambda_v)$이며 λ를 다시 나누지 않는다. TIE-19의 air 식은 아래 adapter와 같은 계수 형식이지만 파장 인자를 vacuum μm로 정의한다. 이 차이도 `air_model_id`에 보존한다. [S29](#ref-S29)

형상 기준 $T_0$와 index fit 기준 $T_c$가 다르면 $T_c$를 임의로 $T_0$로 바꾸지 않는다. 제공된 기준 index가 $T_0$에 있을 경우:

$$
n^{\mathrm{abs}}(\lambda_v,T)=n^{\mathrm{abs}}(\lambda_v,T_0)
 + K_c\{F(\lambda_v,T-T_c)-F(\lambda_v,T_0-T_c)\},
\qquad K_c=\frac{n_c^2-1}{2n_c}.
$$

이 차분 형태는 원래 fit의 기준을 보존하는 직접 대수 변환이다. $K_c$도 fit 기준에서의 값을 유지한다.

#### 2.3 OpticStudio catalog 호환 adapter

다음은 공개 reference의 계산 규약을 구현하기 위한 경로다. 이 profile은 SCHOTT 식을 물리적으로 더 정확하게 만드는 모델이 아니라 **카탈로그 convention을 재현하는 별도 경로**다. $T_s,P_s$는 system 조건, $T_c,P_c=1\,\mathrm{atm}$은 catalog reference, $T_g$는 glass 온도다. 본 과제의 균일 조건에서는 $T_g=T_s=T$다. [S54](#ref-S54)[S57](#ref-S57)[S58](#ref-S58)

공식 reference의 air 식은 다음과 같으며, 이 adapter에서 $\lambda_s$는 system-air μm, $P$는 atm 기준 비율, $T$는 °C다. [S54](#ref-S54)

$$
a_{15}(\lambda_s)=1+10^{-8}\left[
6432.8+\frac{2949810\lambda_s^2}{146\lambda_s^2-1}
+\frac{25540\lambda_s^2}{41\lambda_s^2-1}\right],
$$

$$
a(\lambda_s,T,P)=1+
\frac{[a_{15}(\lambda_s)-1]P}{1+3.4785\times10^{-3}(T-15)}.
$$

실행 순서는 다음과 같다. [S57](#ref-S57)

$$
a_s=a(\lambda_s,T_s,P_s),\quad a_c=a(\lambda_s,T_c,1),\quad
\lambda_c=\lambda_s\frac{a_s}{a_c},
$$

$$
n_{c,\mathrm{rel}}=f_{\mathrm{catalog}}(\lambda_c),\qquad
n_c^{\mathrm{abs}}=a_c n_{c,\mathrm{rel}}.
$$

Thermal term에는 reference-air 파장 $\lambda_c$를 사용한다고 Ansys tutorial이 명시한다. Generic catalog는 signed $\lambda_{TK}$ 규약을 갖는다. [S54](#ref-S54)[S57](#ref-S57)

$$
\Delta n^{\mathrm{abs}}_{\mathrm{AGF}}=
\frac{n_{c,\mathrm{rel}}^2-1}{2n_{c,\mathrm{rel}}}
\left[D_0\delta+D_1\delta^2+D_2\delta^3+
\frac{E_0\delta+E_1\delta^2}{\lambda_c^2-\operatorname{sgn}(\lambda_{TK})\lambda_{TK}^2}\right],
\qquad\delta=T_g-T_c,
$$

$$
n_g^{\mathrm{abs}}=n_c^{\mathrm{abs}}+\Delta n^{\mathrm{abs}}_{\mathrm{AGF}},\qquad
n_{g,\mathrm{system\ rel}}=\frac{n_g^{\mathrm{abs}}}{a_s},\qquad
\lambda_v=\lambda_s a_s.
$$

이 adapter가 내보내는 $n_g^{\mathrm{abs}},a_s,\lambda_v$를 이후 ray/OPL에 사용한다. 고정 진공 파장 실험에서는 각 온도에서 다음 scalar equation을 풀어 $\lambda_s$를 결정한 뒤 adapter를 평가한다.

$$
\lambda_s a(\lambda_s,T_s,P_s)-\lambda_v^{\mathrm{target}}=0.
$$

반대로 `fixed_system_air_wavelength`를 선택하면 GUI 숫자를 고정하므로 진공 파장이 미세하게 변한다. 어느 실험인지 결과에 남긴다. 제공된 absolute 모델과 AGF 모델을 섞어 동일한 온도 보정을 두 번 적용하지 않는다. Thermal 계수가 없다는 사실을 전부 0이라는 물리 측정값으로 해석하지 않는다. `index_temperature_invariant`는 명시적으로 선택한 별도 근사다.

**남은 호환성 검증:** 공개 제조사식과 공개 제품식의 파장 규약 차이는 확인된 사실이다. 어느 쪽을 사용할지는 model ID로 결정할 수 있어 계산이 막히는 문제는 아니다. 다만 실제 목표 release/build의 catalog export와 기준 수치 대조 전에는 완전한 OpticStudio 수치 일치를 인증하지 않는다.

### 3. Transmission 자료에서 흡수 계수 생성

물성의 “transmission range/rate”는 파장 범위, 시험 두께, 내부/외부 투과율 구분이 있어야 한다. 내부 투과율 $\tau_{i,\mathrm{ref}}(\lambda,T)$가 물리 두께 $d_{\mathrm{ref}}>0$에 대해 주어진 경우 Beer–Lambert 모델로 변환한다. [TH02](#ref-TH02)

$$
\mu(\lambda,T)=-\frac{\ln\tau_{i,\mathrm{ref}}(\lambda,T)}{d_{\mathrm{ref}}},\qquad
\tau_{i,\mathrm{ray}}=\exp[-\mu(\lambda,T)\ell_{\mathrm{ray}}]
=\tau_{i,\mathrm{ref}}^{\ell_{\mathrm{ray}}/d_{\mathrm{ref}}}.
$$

여기서 $\ell_{\mathrm{ray}}$는 real ray가 재료 안에서 이동한 실제 길이다. Center thickness를 모든 ray의 흡수 길이로 사용하지 않는다. 각 굴절면 Fresnel/coating 전력 투과율은 별도로 곱한다. 외부 투과율에 포함된 표면 반사를 내부 흡수로 해석하고 다시 Fresnel 손실을 곱하면 이중 계상이다.

$0<\tau_i\le1$이어야 하며, $\tau_i=0$가 측정 검출한계인지 실제 모델의 완전 불투과인지 구분한다. $\tau_i=1$은 이 근사에서 $\mu=0$이다. 파장 보간은 $\tau$ 또는 $\mu$ 중 어느 양에 적용했는지 기록한다. $\mu$ 선형 보간은 비음성을 유지하면서 직접 경로 길이에 적용하기 편한 선택이지만 제조사 지정 보간이 있으면 그것을 따른다.

온도 의존 흡수 자료가 없다면 `mu_temperature_invariant`를 명시하고 $\mu(\lambda,T)=\mu(\lambda,T_0)$로 둔다. 이 경우도 팽창과 굴절 변화가 $\ell_{\mathrm{ray}}$를 바꾸므로 총 내부 투과율은 변할 수 있다. RI 정규화가 이 가정을 없애 주지 않는다.

### 4. 균일 등방성 팽창으로 표면 갱신

#### 4.1 CTE 정의에 따른 선형 치수 비율

$s_j(T;T_0)=L_j(T)/L_j(T_0)$로 정의한다. 데이터의 CTE 정의를 확인하여 한 가지 식을 선택한다. ppm/K 값에는 $10^{-6}$을 곱한다. OpticStudio TCE는 일정 계수의 선형 팽창 모델을 사용한다. [S27](#ref-S27)

| 제공된 계수의 정의 | 사용할 팽창식 |
|---|---|
| 일정한 선형/평균 계수 $\alpha$ | $s=1+\alpha(T-T_0)$ |
| 순간 계수 $\alpha(T)=L(T)^{-1}dL/dT$ | $s=\exp[\int_{T_0}^{T}\alpha(u)du]$ |
| 기준 길이 계수 $\alpha_0(T)=L(T_0)^{-1}dL/dT$ | $s=1+\int_{T_0}^{T}\alpha_0(u)du$ |
| 기준점이 다른 제공 길이 함수 $L(T)$ | $s=L(T)/L(T_0)$ |

위 적분식은 CTE 정의를 적분한 직접 유도다. 구간 평균 CTE를 순간 계수처럼 연속 적분하지 않는다. 재료 자료에 다른 기준 온도가 있으면 먼저 길이 함수를 구성한 뒤 비율로 바꾼다. $s\le0$이면 모델/입력 오류다.

#### 4.2 Standard/Even Asphere의 계수 변환

렌즈 고유 좌표에서 기준 sag를 다음과 같이 쓴다.

$$
z_0(r)=\frac{c_0r^2}{1+\sqrt{1-(1+k_0)c_0^2r^2}}
 +\sum_{m=1}^{6}A_{2m,0}r^{2m}.
$$

모든 좌표가 $s$배 되는 등방성 자유 팽창은 다음과 같다.

$$
z_T(r)=s z_0(r/s),\qquad
c_T=\frac{c_0}{s},\quad R_T=sR_0,\quad k_T=k_0,
\quad A_{2m,T}=A_{2m,0}s^{1-2m}.
$$

이는 원래 sag에 좌표 배율을 직접 대입한 결과다. $A_{2m}$의 단위가 길이$^{1-2m}$이므로 차원도 일치한다. $c_0=0$이면 갱신 후에도 c=0이지만 A항이 남으면 표면 전체가 평면인 것은 아니다. 모든 A항도 0인 실제 평면은 계속 평면이다. 숫자 infinity에 곱하는 구현 대신 curvature를 사용할 수 있다.

| 항목 | 온도 갱신 |
|---|---|
| Radius, center thickness, 렌즈에 속한 형상 직경·rim 높이 | $s$배 |
| Conic | 불변 |
| A2, A4, A6 | 각각 $s^{-1},s^{-3},s^{-5}$배 |
| A8, A10, A12 | 각각 $s^{-7},s^{-9},s^{-11}$배 |

앞·뒷면은 같은 lens material의 동일 $s_j$를 사용한다. Clear aperture가 렌즈 자체의 경계이면 함께 팽창시키되, 별도 고정 diaphragm의 경계이면 해당 diaphragm 정책을 적용한다. Mechanical diameter와 active diameter의 차이만으로 rim 형태·접촉 조건을 새로 추정하지 않는다. 제공된 무효부 형상은 같은 렌즈 고유 좌표에서 갱신한다.

각 온도의 표면은 항상 원래 $T_0$ 처방에서 생성한다. $T_1\to T_2\to T_3$의 누적 배율을 기본 구현으로 사용하지 않는다. 온도 sweep의 순서가 결과에 영향을 주어서는 안 된다.

### 5. Axial anchor로 global vertex와 air gap 계산

#### 5.1 기준 처방의 vertex 생성

표면 $i$의 기준 vertex 좌표를 $v_{i,0}$, 해당 surface의 LDE thickness를 $t_{i,0}$라 하면:

$$
v_{i+1,0}=v_{i,0}+t_{i,0}.
$$

렌즈 $j$는 front/back vertex $v_{f_j,0},v_{b_j,0}$와 중심 두께 $CT_{j,0}=v_{b_j,0}-v_{f_j,0}>0$를 가진다. Dummy/STOP/IMA surface는 렌즈의 굴절면 수와 별개다.

#### 5.2 일반 anchor 식

렌즈마다 0…1의 고정 비율 $a_j$를 사용해 anchor를 정의한다.

$$
b_{j,0}=v_{f_j,0}+a_jCT_{j,0},\qquad
v_{f_j}(T)=b_j(T)-a_js_jCT_{j,0},\qquad
v_{b_j}(T)=b_j(T)+(1-a_j)s_jCT_{j,0}.
$$

- $a_j=0$, $b_j(T)=v_{f_j,0}$: front vertex 고정.
- $a_j=1/2$, $b_j(T)=(v_{f_j,0}+v_{b_j,0})/2$: **축상 두 vertex의 중점** 고정. 질량 중심과 동일하다고 가정하지 않는다.
- $a_j=1$, $b_j(T)=v_{b_j,0}$: back vertex 고정.
- 명시된 $b_j(T)$: 사용자가 정의한 위치 함수. 본 과제에서 그 함수의 기계적 발생 원인을 FEA로 계산하지 않는다.

형상은 이 vertex 위에 갱신된 local sag를 놓아 생성한다. 렌즈 사이 공기 간격은 반드시 결과 vertex에서 유도한다.

$$
g_j(T)=v_{f_{j+1}}(T)-v_{b_j}(T).
$$

독립적인 고정 air-gap 모델을 원하면 anchor 식과 동시에 강제하지 말고, 첫 lens anchor를 정한 뒤 $g_j(T)$를 입력 규칙으로 사용해 후속 lens를 순차 배치한다. 즉 위치 제약 그래프에 독립 조건만 둔다. 전역 고정 IMA와 마지막 rear→IMA thickness 고정도 일반적으로 동시에 만족하지 않는다.

열 갱신 뒤에는 baseline과 동일한 전체 구경 검사, 앞·뒷면의 내부 교차, lens-lens 충돌, clear/mechanical aperture 관계를 다시 검사한다. Center/edge 두 점의 두께가 양수인 것만으로 유효 렌즈임을 보장하지 않는다.

#### 5.3 STOP·OBJ·IMA 정책

| 요소 | 실행 규칙 |
|---|---|
| 렌즈 surface에 지정된 STOP | 그 surface의 vertex/shape를 따른다. Aperture 반경은 렌즈 부착/독립 aperture 중 선택한 정책을 따른다. |
| 독립 STOP dummy surface | 별도 $v_{\mathrm{stop}}(T),a_{\mathrm{stop}}(T)$를 평가한다. 렌즈 CTE를 임의로 적용하지 않는다. |
| EPD/F-number 입력 | 우선 $T_0$에서 실제 STOP 크기로 resolve하고 온도별 물리 STOP 규칙을 적용한다. 매 온도 F/#를 강제로 유지하는 실험은 `recontrolled_aperture`로 분리한다. |
| 유한 OBJ | 기본 열 실험에서는 전역 위치와 12.8×3.2 mm 발광면 좌표를 고정한다. 이동 보상변수를 사용하면 별도 결과로 저장한다. |
| 무한 OBJ | 입력 방향과 진공 파장을 고정한다. 온도별 field를 다른 입사 방향으로 자동 재맞춤하지 않는다. |
| 고정 IMA | 전역 위치를 고정한다. 마지막 간격은 갱신된 rear vertex에서 다시 계산한다. |
| 재초점 IMA | 고정 IMA 결과를 먼저 계산한 뒤 하나의 공통 이동량을 구한다. |

EPD/F/#를 매 온도 재고정하거나 온도별 FOV에 맞춰 field를 다시 정규화하면 서로 다른 물리 실험이 된다. 의도적으로 선택할 수 있으나 고정 hardware 열 성능으로 혼동하지 않는다. Source radiance/spectrum의 온도 변화는 이번 모델에 주어지지 않았으므로 광원은 동일하게 둔다.

---

<a id="part-3"></a>

## 3. 실광선 커널·초점거리·기하 수차

이 장의 `제안 규약`은 독립 엔진이 실제로 실행할 알고리즘이다. `제품 정의`는 공식 자료에서 확인한 OpticStudio의 의미다. 두 표시를 구별하며, 여기의 알고리즘이 비공개 vendor solver와 동일하거나 수치 대조를 통과했다는 뜻은 아니다. 전제는 동축 3매, 균질·등방성 재질, 6개 Standard/Even Asphere 굴절면, 직접 순차 투과다. STOP용 무굴절 면과 OBJ/IMA는 굴절면 수에 포함하지 않는다. 제조 공차·반사 ghost·산란·응력·GRIN은 계산하지 않는다.

### G1. 입력을 실제 처방으로 해석하는 순서

내부 길이는 mm, 각도는 rad, 방향은 단위벡터, 재료 굴절률은 absolute index, 위상 파장 λ₀는 vacuum wavelength를 mm로 변환하여 통일한다. G라는 이름만 입력받지 않고 실제 수치 파장을 요구한다. 온도 모듈은 해당 온도의 완전한 처방과 재료값을 이 커널에 전달한다. `surface_id`와 `refracting_surface_index`를 구별한다.

1. 앞면 thickness는 그 렌즈의 두 vertex 간 중심 두께, 뒷면 thickness는 다음 렌즈 앞면 vertex까지의 간격으로 해석한다. 마지막 굴절면 뒤 거리와 IMA 위치는 하나의 편집 원천에서 유도한다.
2. 전역 +z는 OBJ에서 IMA로 향한다. 모든 면의 국소 z축도 +z이며 vertex는 $\mathbf v_j=(0,0,z_j)$다. signed radius R의 양수는 곡률 중심이 vertex의 +z 쪽이라는 뜻이다. 평면은 `curvature=0`으로 저장하고 JSON Infinity에 의존하지 않는다.
3. 면 j의 입력은 surface type, c=1/R, k, A₂…A₁₂, vertex, 광학 유효 구경, 기계 구경, rim 모델, 입사/출사 재질이다. Standard에서도 k를 허용하며 A 계수만 모두 0이다. [S01](#ref-S01)[S02](#ref-S02)
4. 하나의 system aperture 원천만 선택한다: `stop_size`, `entrance_pupil_diameter`, `nominal_image_f_number`. 나머지는 유도하거나 일관성 확인용이다. STOP 위치는 별도 필수 입력이다. [S44](#ref-S44)
5. canonical field는 무한 물체의 방향 또는 유한 OBJ의 실제 물점이다. 원래 GUI field type과 단위·normalization은 별도 메타데이터로 보존한다. 유한 물점과 그 물점의 chief-ray 방향을 동시에 독립 조건으로 고정하지 않는다. chief는 STOP 중심 조건으로 결정된다.
6. 영상 모드는 `focal_plane` 또는 `afocal_direction`이다. Afocal 설정은 광학계 power가 0이라는 선언이 아니다.

원형 active aperture의 반경은 a=D_active/2, mechanical 반경은 b=D_mech/2다. 광선 차단 판단과 형상 표시를 분리한다. **제안 기본 rim 규약은 active 영역 밖을 opaque로 처리**하는 것이다. 무효부 높이는 기하 표시·간섭 검사에 사용하지만 그 자체로 투명 rim의 굴절 형상을 유일하게 정하지 않는다. 투명 무효부까지 추적하려면 그 구간의 piecewise sag·재질을 추가해야 한다. 단순 높이값으로 그 면을 추정하지 않는다. [S37](#ref-S37)

### G2. Sag, gradient, Hessian

다음은 Even Asphere의 공개 sag와 그 직접 미분이다. $q=x^2+y^2$, $\gamma=(1+k)c^2$, $s=\sqrt{1-\gamma q}$로 둔다. [S01](#ref-S01)

$$
z(q)=\frac{cq}{1+s}+P(q),\qquad P(q)=\sum_{m=1}^{6}A_{2m}q^m.
$$

다항식과 도함수는 q에 대해 Horner 방식으로 계산한다. $(1-s)/[(1+k)c]$ 형태는 작은 q에서 cancellation이 생기므로 기본 sag 식으로 사용하지 않는다. c=0에서도 위 식은 그대로 유효하다. $1-\gamma q<0$이면 sag 정의역 밖이다. s=0에서 sag는 유한할 수 있지만 gradient가 발산하므로 매끄러운 surface 커널의 지원 경계로 취급한다.

$$
B(q)=\frac{c}{s}+\sum_{m=1}^{6}2m A_{2m}q^{m-1},
\qquad \nabla z=B(q)\begin{bmatrix}x\\y\end{bmatrix}.
$$

$$
C(q)=\frac{c\gamma}{s^3}+\sum_{m=2}^{6}4m(m-1)A_{2m}q^{m-2},
\qquad H_z=B I_2+C\begin{bmatrix}x^2&xy\\xy&y^2\end{bmatrix}.
$$

이 표현은 r=0에서 z′(r)/r를 직접 나누지 않아도 된다. 특히 다음 정점 곡률은 1차 광학에 사용한다.

$$
z(r)=\left(\frac c2+A_2\right)r^2+O(r^4),
\qquad c_v=z''(0)=c+2A_2.
$$

따라서 A₂가 있는 처방의 paraxial power를 c만으로 계산하지 않는다. 계수 A₂m의 단위는 mm^(1−2m)이며 렌즈 단위 변경 시 이 차원을 적용한다.

면의 implicit 함수와 미분은 다음과 같다.

$$
F_j(\mathbf r)=r_z-z_j-z(r_x,r_y),\quad
\mathbf g_j=\nabla F_j=(-z_x,-z_y,1)^T,\quad
H_{F_j}=\begin{bmatrix}-H_z&0\\0&0\end{bmatrix}.
$$

추적 전에는 유효/기계 모델 전체 구경에서 정의역·재료 순서·surface crossing을 확인한다. 렌즈 앞뒤 두 면의 두께 함수

$$
t(r)=t_c+z_{\rm back}(r)-z_{\rm front}(r)
$$

가 공통 지지 영역에서 양수인지 검사한다. 중심과 가장자리 두 점만으로 판단하지 않는다. q 구간을 세분하여 interval bounds로 t>0을 증명하거나, stationary point들을 분리해 최솟값을 확인한다. 증명/수치 결정 예산이 부족하면 `geometry_unresolved`로 남긴다. 이 장의 interval root 선택과 구간 검사 절차는 독립 구현 제안이며 vendor 내부 알고리즘이라는 출처는 없다.

### G3. 순차 면 교점: 첫 물리적 해를 선택한다

입사 ray는 $\mathbf r(t)=\mathbf p+t\mathbf d$, $\|\mathbf d\|=1$이고 t는 mm 단위 물리 경로 길이다. 다음 순차 면에서 풀 방정식은

$$
f(t)=F_j(\mathbf p+t\mathbf d)=0,\qquad
f'(t)=\nabla F_j\cdot\mathbf d.
$$

Newton 후보는

$$
t_{k+1}=t_k-\frac{f(t_k)}{f'(t_k)}.
$$

이는 공개된 일반 ray tracing 접근과 같다. 다만 수렴한 한 해가 최초의 물리적 해라는 보장은 없으므로 다음 **bracketed first-hit 절차**를 기본으로 한다. [S03](#ref-S03)

1. 해당 면의 mechanical cylinder와 sag z-range가 만드는 유한 bounding volume에 ray가 머무르는 t 구간을 얻는다. 가능한 forward t만 남긴다. Bounding range는 사전 interval sag bounds로 만든다. 정확히 연속된 dummy 면은 별도 zero-thickness 규칙을 사용한다.
2. 후보 구간을 t가 작은 순서로 처리한다. 구간 arithmetic로 f의 포함 범위를 계산해 0이 없으면 버린다. 실제 interval 연산은 outward rounding을 사용해야 한다. 일반 float min/max 샘플을 엄밀한 bound라고 부르지 않는다.
3. f′ interval에 0이 없으면 f는 단조다. 양 끝 f가 같은 부호이면 근이 없고, 부호가 다르면 유일한 근을 bracket한다. 근이 endpoint에 있는 경우는 endpoint와 다음 구간의 중복을 제거한다.
4. 단조 bracket에서는 Newton step이 bracket 안에 있고 유한하며 잔차를 줄일 때 사용하고, 아니면 midpoint로 이분한다. bracket 폭과 잔차 기반 위치 오차를 모두 수렴시킨다.
5. f′ interval에 0이 있으면 구간을 이분한다. 접선근/다중근은 부호 변화만으로 발견할 수 없다. 예산 내 분리하지 못하면 `intersection_unresolved` 또는 `grazing_degenerate`로 반환한다. 뒤의 쉬운 근으로 건너뛰지 않는다.
6. 선택한 근보다 앞선 모든 후보 구간이 제거됐을 때만 first hit로 확정한다. 원래 sag 식, 정의역, 순방향 거리, 교점 오차, 재료 반공간을 재확인한다. conic을 제곱한 다항식의 다른 branch를 채택하지 않는다.
7. 첫 교점에 active mask를 적용한다. 앞의 opaque/rim hit를 버리고 뒤의 active hit를 찾는 행위는 금지한다. 표면 바깥으로 빠져나간 ray는 물리적 escape이고, 근을 찾지 못한 ray는 unresolved다.

평면은 해석적 t를, 구면/conic은 안정적인 이차방정식 후보를 사용할 수 있다. 그러나 후보마다 원래 sag branch와 first-hit 조건을 확인한다. 정상적인 처방에서는 conic seed+Newton이 빠른 경로이며 위 검사를 fallback/검증으로 결합할 수 있다. 임의의 고차 비구면 계수 전 영역에 대해 제한 iteration만으로 해가 보장된다고 주장하지 않는다.

잔차만 보고 정확도를 정하지 않는다. 단순근에서 위치 오차의 1차 추정은

$$
\epsilon_t\simeq \frac{|f(t)|}{|\nabla F\cdot\mathbf d|}.
$$

분모가 작은 grazing에서는 작은 잔차도 큰 위치 오차일 수 있다. bracket 폭과 $|\nabla F\cdot\mathbf d|$, root isolation 상태를 함께 저장한다. 허용 오차는 전체 OPD/영상좌표 오차 예산에서 정하며 보편적인 임의 숫자를 고정하지 않는다.

### G4. Vector Snell, TIR, OPL 및 ray 기록

교점의 단위 normal을 입사 방향 반대쪽으로 정한다.

$$
\widehat{\mathbf g}=\frac{\mathbf g}{\|\mathbf g\|},\qquad
\mathbf N=-\operatorname{sign}(\widehat{\mathbf g}\cdot\mathbf d_i)\widehat{\mathbf g},
\qquad c_i=-\mathbf N\cdot\mathbf d_i\ge0.
$$

단, 법선 뒤집기가 재료 오류를 숨겨서는 안 된다. 입력한 입사·출사 매질 및 허용 순차 방향이 이 반공간과 일치해야 한다. $\widehat{\mathbf g}\cdot\mathbf d_i=0$은 grazing 별도 상태다. $\eta=n_i/n_t$일 때 물리적 투과 해는 [S06](#ref-S06)

$$
Q=1-\eta^2(1-c_i^2),\qquad
\mathbf d_t=\eta\mathbf d_i+(\eta c_i-\sqrt Q)\mathbf N.
$$

Q<0이면 TIR이며 이 직접 투과 경로는 종료한다. 임계각 부근의 roundoff bound에 의해 Q의 부호를 확정하지 못하면 `critical_angle_unresolved`로 남기거나 더 높은 정밀도로 재계산한다. 임의로 모든 음수 Q를 0으로 자르지 않는다. 다음 세 조건을 함께 확인한다.

$$
\|\mathbf d_t\|=1,\quad
n_i[\mathbf d_i-(\mathbf d_i\cdot\mathbf N)\mathbf N]
=n_t[\mathbf d_t-(\mathbf d_t\cdot\mathbf N)\mathbf N],\quad
\mathbf d_t\cdot\mathbf N=-\sqrt Q\le0.
$$

마지막 부호 조건이 잘못된 root를 배제한다. n_i=n_t인 dummy 면은 동일 방향을 유지한다. unit vector normalization은 미세한 부동소수점 오차를 교정할 때만 수행하고 큰 이상을 숨기지 않는다.

$$
\mathrm{OPL}_{\rm abs}=\sum_j n_{{\rm abs},j}\ell_j.
$$

물리 추적의 $\ell_j$는 실제 전진 길이이다. 파면 기준 구면/평면까지의 가상 signed propagation은 wave 모듈의 별도 단계이며 이 ray ledger와 혼합하지 않는다. 무한 물체의 유한 launch plane에서는 plane-wave 시작 위상 $S_{\rm in}=n_o\mathbf d\cdot(\mathbf p-\mathbf p_{\rm ref})$를 wave 모듈에 함께 전달한다. 유한 물점은 공통 시작 위상을 가진다. raw OPL 차이만으로 OPD를 만들지 않는다.

각 ray의 최소 기록은 `field_id, wavelength, sample_id, launch position/direction, input phase, sampling coordinates/measure/weight`, 각 면의 `hit, d_in, d_out, n_in/out, physical length, cumulative OPL, normal, incidence cosine, stop coordinate, aperture status, solver residual/bracket/conditioning`, 최종 `status, last valid surface`다. 파워/편광 모듈이 전달한 transmission은 ray의 기하 sampling weight와 다른 필드다.

#### G4.1 파면·radiometry 모듈로 전달하는 공통 평면 API

`trace_to_plane(field, λ₀, T, launch_parameters, plane)`은 동일 trace 뒤 최종 균질 매질에서 명시 평면 Π까지 전진한다. Π 원점 r_Π, 단위 법선 n_Π와 basis e₁,e₂가 주어지면

$$
t_\Pi=\frac{\mathbf n_\Pi\cdot(\mathbf r_\Pi-\mathbf p)}{\mathbf n_\Pi\cdot\mathbf d},\quad
\mathbf r_{hit}=\mathbf p+t_\Pi\mathbf d,\quad
S_{\Pi}=S_{in}+\sum_j n_{abs,j}\ell_j+n_{i,abs}t_\Pi.
$$

실제 wave boundary plane 경로는 t_Π≥0이고 마지막 모든 굴절면 뒤 같은 균질 매질에 있어야 한다. 반환값은 transverse (u,v), direction, S_Π(mm), physical power/transmission ledger, mask와 실패 상태, 가능한 경우 launch→(u,v) Jacobian이다. 무한 물체의 S_in과 마지막 propagation을 각각 한 번만 더한다.

`aim_to_exit_plane(target_u,target_v)`은 G6의 launch 변수를 그대로 쓰되 residual을 STOP 위치 대신 위 평면의 (u,v) 오차로 바꾼다. 같은 damped Newton/branch 검증을 사용하고 실제 STOP을 포함한 모든 physical aperture를 최종 trace에서 시험한다. wave 모듈이 direction-cosine grid를 지정하면 먼저 그 모듈의 선언된 mapping으로 (u,v) target을 구한다. STOP 좌표에 균일한 표본과 exit-plane/cosine grid에 균일한 표본의 weight가 같다고 가정하지 않는다. 평면 mapping의 fold·다중해·미결정 root는 mask=0이 아니라 wave 모델의 미해결/미지원 상태로 보낸다.

### G5. Ray Jacobian: aiming과 parabasal 분석에 공유

다음은 고정 처방에서 launch parameter $\xi$에 대한 직접 미분이다. finite difference를 기준 구현으로 먼저 사용해도 되지만 동일 branch·aperture 내부에서 step halving과 비교한다. production에서 미분 ray를 사용하면 intersection·Snell·aiming·field curvature를 하나의 일관된 경로로 연결할 수 있다.

교점 t와 $\mathbf h=\mathbf p+t\mathbf d$에 대해

$$
\frac{\partial t}{\partial\xi}
=-\frac{\mathbf g^T(\partial_\xi\mathbf p+t\partial_\xi\mathbf d)}{\mathbf g^T\mathbf d},
\qquad
\partial_\xi\mathbf h=\partial_\xi\mathbf p+t\partial_\xi\mathbf d+\mathbf d\,\partial_\xi t.
$$

normal 부호가 국소적으로 일정한 branch에서

$$
\partial_\xi\widehat{\mathbf g}
=\frac{I-\widehat{\mathbf g}\widehat{\mathbf g}^T}{\|\mathbf g\|}
H_F\,\partial_\xi\mathbf h.
$$

고정 n의 Snell 미분은

$$
\partial_\xi c_i=-(\partial_\xi\mathbf N)^T\mathbf d_i-\mathbf N^T\partial_\xi\mathbf d_i,
\qquad \partial_\xi Q=2\eta^2c_i\partial_\xi c_i,
$$

$$
\partial_\xi\mathbf d_t=\eta\partial_\xi\mathbf d_i
+\left(\eta\partial_\xi c_i-\frac{\partial_\xi Q}{2\sqrt Q}\right)\mathbf N
+(\eta c_i-\sqrt Q)\partial_\xi\mathbf N.
$$

이 미분은 Q>0인 매끄러운 transmitted branch에 한정한다. grazing, TIR 경계, aperture 경계, root 전환에서 Jacobian은 불연속/발산할 수 있다. 그곳에서 유한 숫자를 만든다고 민감도가 유효한 것은 아니다. 설계 파라미터 미분은 $F_\xi$, $\eta_\xi$, pose 미분 등 추가 항이 필요하며, 위 ray-parameter 미분을 그대로 재사용하지 않는다.

### G6. STOP을 맞추는 실제 ray aiming

#### G6.1 Field와 발사 변수

무한 물체의 rectangular field angle을 $(\theta_x,\theta_y)$로 받는 제안 표현은

$$
\mathbf d_o=\frac{(\tan\theta_x,\tan\theta_y,1)}{\sqrt{1+\tan^2\theta_x+\tan^2\theta_y}}.
$$

이 장은 각 성분이 forward chart에 있는 경우를 지원한다. field normalization은 이 실제 angle/point로 변환한 후 적용하지 않는다. 무한 물체 aiming의 미지수는 앞쪽 launch plane의 $\boldsymbol\xi=(x_L,y_L)$이고 $\mathbf d_o$는 고정이다.

유한 물체는 $\mathbf p_o=(x_o,y_o,z_o)$를 고정하고 launch slope $\boldsymbol\xi=(s_x,s_y)$를 미지수로 둔다.

$$
\mathbf d_o(\boldsymbol\xi)=\frac{(s_x,s_y,1)}{\sqrt{1+s_x^2+s_y^2}},\qquad
\frac{\partial\mathbf d_o}{\partial s_a}
=\frac{(I-\mathbf d_o\mathbf d_o^T)\mathbf e_a}{\sqrt{1+s_x^2+s_y^2}}.
$$

Finite `Angle` GUI 입력은 entrance pupil을 기준으로 물점을 먼저 정한 뒤 실제 aiming에서 방향이 바뀔 수 있다. 예를 들어 축상 entrance pupil vertex가 z_EP라면 $x_o=-(z_{EP}-z_o)\tan\theta_x$, y도 동일하다. 이 변환은 평면 OBJ·동축 paraxial pupil을 전제로 한다. 실제 Object Height를 기본 canonical 입력으로 쓰면 이 중간 규약을 피할 수 있다. 제품도 finite Angle에서 ray aiming이 실제 시작 각도를 바꾼다고 설명한다. [S48](#ref-S48)

#### G6.2 Aperture 입력을 STOP 목표 크기로 변환

STOP을 독립 물리 크기로 지정한 경우 a_stop을 그대로 사용한다. EPD 또는 nominal F/# 입력일 때에는 **제안 기본을 paraxial stop-sizing + actual ray aiming**으로 한다. 이름에 paraxial이 있어도 이후 여섯 면의 실제 추적은 근축 근사가 아니다. 제품의 Paraxial Ray Aiming도 STOP의 목표 높이 산정에 paraxial marginal ray를 사용할 수 있다. [S04](#ref-S04)

첫 굴절면 직전에서 STOP vertex까지의 reduced-ray 행렬을 $M_s=\left[\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right]$라 하자. STOP이 굴절면에 있으면 그 면에 도달한 높이는 굴절 전후 동일하다. a≠0일 때

$$
z_{EP}=z_1+\frac{n_o b}{a},\qquad
a_{stop}=|a|\frac{D_{EP}}2.
$$

이는 STOP의 paraxial entrance-pupil image를 구하는 직접 행렬 유도다. a=0/거의 0은 finite EPD 표현이 부적합하거나 ill-conditioned하므로 stop-size 입력 경로로 전환하도록 상태를 반환한다. 자동으로 임의 큰 pupil을 만들지 않는다. nominal F/#의 positive convention은 $N_{nom}=|F_{report}|/D_{EP}$이며, F_report의 index 기준은 G7에서 정의한다. 이 quantity는 지정 공액의 working F/#와 다르다. 제품 Image Space F/#도 infinite-conjugate paraxial EFL/EPD를 사용한다. [G01](#ref-G01)

이미 고정된 STOP clear aperture가 유도 a_stop과 다르면 어느 값을 수정/유도할지 입력 ownership 규칙으로 해결한다. 다른 면의 active aperture는 그대로 적용하므로 설정한 nominal cone이 모두 통과한다는 뜻은 아니다. 실제 통과 pupil과 nominal pupil은 결과에 따로 남긴다.

#### G6.3 2D damped Newton 절차

정규화 sample $\boldsymbol\rho=(\rho_x,\rho_y)$, $\|\boldsymbol\rho\|\le1$의 STOP 목표는 $\mathbf s=a_{stop}\boldsymbol\rho$다. chief는 $\boldsymbol\rho=0$, 네 기준 marginal ray는 (±1,0),(0,±1)이다. STOP은 실제 surface sag에 있으므로 x/y 목표를 정하고 z는 sag에서 얻는다.

$$
\mathbf R(\boldsymbol\xi)=\mathbf h_{stop,xy}(\boldsymbol\xi)-\mathbf s,
\qquad J=\frac{\partial\mathbf h_{stop,xy}}{\partial\boldsymbol\xi}.
$$

1. paraxial 선형식으로 chief seed를 만든다. OBJ→STOP 행렬이 A,B이면 finite seed는 축별 $u_o=(s-Ay_o)/B$, slope≈u_o/n_o다. 무한 물체에서는 launch→STOP 행렬의 a,b로 $y_L=(s-bu_o)/a$를 쓴다. 분모가 특이하면 alternative seed/명시 unsupported 상태를 사용한다.
2. 원점 sample을 먼저 풀고, 가까운 pupil sample을 이용하여 center→edge로 continuation한다. field도 작은 field→큰 field로 이어갈 수 있다. warm start는 수렴/branch 확인을 대체하지 않는다.
3. STOP까지 실제 surface 교점·Snell을 수행한다. aiming 중에는 active clipping을 적용하지 않는 `geometry_only_aim` 경로로 residual을 평가할 수 있으나 sag 정의역·TIR·기계적 branch는 유지한다. 최종 물리 trace에서 모든 mask를 다시 적용한다. 기계 구경 밖으로 무한 sag를 연장하려면 별도 명시 옵션이며 실제 통과 ray로 채택하지 않는다.
4. J를 G5 또는 central difference로 구한다. 잘 conditioned하면 $J\Delta\xi=-R$를 pivoted solve로 풀고, 나쁘면 $(J^TJ+\mu I)\Delta\xi=-J^TR$의 damped step을 사용한다. 서로 다른 차원의 미지수는 seed scale로 무차원화한다.
5. step을 1, 1/2, 1/4…로 줄여 같은 물리 branch에서 residual이 감소하는 첫 trial을 채택한다. 유효 trial이 없거나 iteration budget을 소진하면 `aim_unresolved`다.
6. STOP residual과 launch 변화가 목표 오차 예산 안에 들고 J 조건이 기록됐을 때 수렴으로 반환한다. 서로 다른 초기값에서 여러 해가 생기거나 map fold가 확인되면 `multi_branch_aim`이다. 기본 단일-branch pupil 프로필에서는 그 field의 wave 분석을 미지원/무효 처리하며 임의 해를 섞지 않는다.

`no physical ray`와 `solver failed to find ray`는 다르다. 실패한 aiming sample을 단순 vignetted ray로 세면 RI와 MTF가 함께 왜곡된다. 수동 vignetting factor는 pupil의 shift/compression을 바꾸는 설정이며 실제 surface clipping과 별개다. 기본은 0이고, 수입한 제품 factor를 지원할 때에는 정확한 sample transform을 별도 adapter에 둔다. [S04](#ref-S04)[S05](#ref-S05)

#### G6.4 Working F/# 진단

단순 축상 실제 cone은 $N_w=1/(2n_i\sin u_i)$다. n의 absolute/relative 기준을 이름에 포함한다. 제품과 비교할 때에는 system-air-relative n을 사용한다. 비축 product convention의 구조는 axis ray와 네 marginal ray 사이 각 $u_k$로

$$
\mathrm{NA}_{rms}^2=\frac14\sum_{k=1}^{4}(n_{i,rel}\sin u_k)^2,
\qquad N_{w,report}=\frac{1}{2\sqrt{\mathrm{NA}_{rms}^2}}.
$$

제품은 이 quantity에서 surface apertures를 무시하고 vignetting factors를 고려한다. 추적 불능 marginal의 축소-pupil 추정과 큰 F/# cap 등 제품 예외까지 기본 엔진이 모방한다고 보지는 않는다. [S49](#ref-S49) 따라서 `actual passed pupil`, throughput, T/S directional support를 별도로 제공하며 이 평균값 하나로 비축 FFT 모든 축의 scale을 정하지 않는다.

### G7. EFL, BFL 및 paraxial conjugate

EFL/BFL은 정의 자체가 1차 광학량이다. 다른 성능이 real ray 기반이라고 이 지표까지 full-pupil best focus로 바꾸지 않는다. 고정 파장에서 $(y,u=n_{abs}\theta)$를 상태로 사용한다.

$$
T(t,n)=\begin{bmatrix}1&t/n\\0&1\end{bmatrix},\qquad
S(n_1,n_2,c_v)=\begin{bmatrix}1&0\\-(n_2-n_1)c_v&1\end{bmatrix}.
$$

$$
M=S_6T_5S_5\cdots T_1S_1
=\begin{bmatrix}A&B\\C&D\end{bmatrix}.
$$

M은 첫 굴절면 직전부터 마지막 굴절면 직후까지이며 마지막→IMA 이동은 포함하지 않는다. intermediate dummy와 propagation은 필요한 위치에 삽입한다. 동일 refracting surface에 STOP marker를 붙였다고 S를 두 번 적용하지 않는다. 각 행렬의 determinant=1과 전체 determinant≈1은 구현 확인식이다.

$$
\Phi_{abs}=-C,\quad
f'_{physical}=-\frac{n_{i,abs}}{C},\quad
\mathrm{BFL}_{last}=-\frac{n_{i,abs}A}{C}.
$$

**Absolute index와 제품 EFL 보고 기준을 연결하는 보완.** OpticStudio는 system T/P의 공기 index를 1로 정규화한다. [S54](#ref-S54)[S57](#ref-S57) $n_{rel}=n_{abs}/n_a$라면 같은 광학계의 reduced state 변환에 의해

$$
M_{rel}=\begin{bmatrix}A&n_aB\\C/n_a&D\end{bmatrix},\qquad
F_{report}=-\frac{n_a}{C},\qquad
f'_{physical}=n_{i,rel}F_{report}.
$$

이 식은 제품의 reciprocal-power EFL 및 image-space EFL 구별 [S07](#ref-S07)과 index 정규화의 직접 결합 유도다. $-1/C_{abs}$를 무조건 Zemax EFL로 출력하지 않는다. $n_a=1$을 명시한 공기 근사에서는 기존 단순식 $-1/C$가 된다. 각 파장에서 n_a도 해당 파장값을 사용한다.

제품 BFL은 보통 마지막 glass surface에서 paraxial focus까지다. 비공기 image-space의 제품 예외에서는 IMA부터 측정한다. 마지막 굴절면→IMA 거리가 d_I이고 그 사이 동일 매질이면 [S08](#ref-S08)

$$
\mathrm{BFL}_{product}=\mathrm{BFL}_{last}-d_I.
$$

해당 비공기 예외일 때만 위 원점 변환을 적용하고 두 quantity를 모두 남긴다. 이는 문서에 공개된 원점 정의의 직접 변환이며 특정 build의 실제 보고값·마지막 glass surface 식별까지 수치 검증한 adapter는 아니다. C≈0은 power 및 uncertainty로 판정한다. 정확한 C=0이면 유한 cardinal focus 없음이 정상 결과다. 수치적으로 0 근처이면 매우 큰 숫자만 내지 말고 ill-conditioned 상태를 남긴다. 유한 물점에서 collimated output이 되는 설정은 C=0을 뜻하지 않는다.

유한 OBJ의 paraxial image는 $M_o=M T(z_1-z_o,n_o)=\left[\begin{smallmatrix}A_o&B_o\\C_o&D_o\end{smallmatrix}\right]$에서

$$
\ell_{image}=-\frac{n_iB_o}{D_o},\qquad
m=A_o+\frac{\ell_{image}}{n_i}C_o.
$$

D_o=0이면 finite object의 paraxial image가 무한대다. BFL은 무한 공액 지표이고 위 유한 공액 image 위치나 full-pupil best focus와 혼합하지 않는다. [S07](#ref-S07)[S08](#ref-S08)[S09](#ref-S09)

### G8. Spot diagram, T/S 투영, angular spot

Focal 분석은 모든 field를 동일 IMA에 교차시킨다. 평면 IMA의 원점 $\mathbf r_I$, 법선 $\mathbf n_I$, 직교 basis $\mathbf e_x,\mathbf e_y$에 대해

$$
t_I=\frac{\mathbf n_I\cdot(\mathbf r_I-\mathbf p)}{\mathbf n_I\cdot\mathbf d},\qquad
\mathbf x_i=\begin{bmatrix}\mathbf e_x^T\\\mathbf e_y^T\end{bmatrix}
(\mathbf p+t_I\mathbf d-\mathbf r_I).
$$

기본은 전역 z에 수직인 평면이다. $\mathbf n_I\cdot\mathbf d$≈0, 뒤쪽 가상 IMA hit, 물리 aperture hit 여부를 별도 상태로 저장한다. Extended plane의 geometric analysis가 virtual intercept를 허용한다면 `virtual_intercept`를 명시하며 물리 screen power와 혼합하지 않는다.

Sample weight w는 선택한 pupil sampling measure와 분석 정책이다. 균일 STOP sample, 균일 entrance-pupil sample, physical source power는 같은 분포가 아니다. Spot 결과에 그 이름과 실제 quadrature weight를 남긴다. [S10](#ref-S10)

$$
\bar{\mathbf x}=\frac{\sum_iw_i\mathbf x_i}{W},\quad W=\sum_iw_i,\qquad
R_{RMS}=\sqrt{\frac1W\sum_iw_i\|\mathbf x_i-\mathbf x_{ref}\|^2},\qquad
R_{GEO}=\max_i\|\mathbf x_i-\mathbf x_{ref}\|.
$$

`reference=centroid` 또는 `chief`를 명시한다. GEO는 추적한 표본의 최댓값이지 연속 pupil의 엄밀한 상한이 아니다. Centroid covariance는

$$
\Sigma=\frac1W\sum_iw_i(\mathbf x_i-\bar{\mathbf x})(\mathbf x_i-\bar{\mathbf x})^T,
\quad \sigma_T^2=\mathbf e_T^T\Sigma\mathbf e_T,\quad
\sigma_S^2=\mathbf e_S^T\Sigma\mathbf e_S.
$$

동축 계의 radial/meridional 표시는 $\mathbf e_T$를 field meridian의 영상면 방향, $\mathbf e_S$를 그 직교 방향으로 둔다. 축상에서는 T=y,S=x의 고정 규약을 쓴다. **이 radial spot 표시는 FFT MTF의 object-local y/x T/S label과 별개**다. 모든 출력에 basis vector를 남겨 그림 회전을 재현한다. [S60](#ref-S60)

Afocal spot은 출사 위치가 아니라 방향의 분산이다. 분석이 고정 vendor 축을 요구하지 않는 한 transverse basis는 global X를 chief 수직면에 투영·정규화한 e₁, e₂=d_c×e₁로 정한다. 투영이 특이하면 global Y를 사용하고 선택을 기록한다. 제안 angular chart는 chief 방향 $\mathbf d_c$와 그에 수직인 $\mathbf e_1,\mathbf e_2$를 사용한다.

$$
\alpha_i=\operatorname{atan2}(\mathbf d_i\cdot\mathbf e_1,\mathbf d_i\cdot\mathbf d_c),\qquad
\beta_i=\operatorname{atan2}(\mathbf d_i\cdot\mathbf e_2,\mathbf d_i\cdot\mathbf d_c).
$$

동일 식으로 angular centroid/covariance/RMS를 계산하고 rad/mrad/degree를 명시한다. 전방 hemisphere의 regular chart를 전제로 한다. 이는 명확한 독립 chart이며 모든 Zemax angular 표시의 동일 convention이라고 주장하지 않는다. 자유공간 IMA 이동은 $\mathbf d_i$를 바꾸지 않아 angular spot을 refocus하지 못한다.

### G9. Real-ray 수차: 팬, 종수차, coma와 색수차

각 field에서 실제 pupil diameter를 따라 q∈[−1,1]을 추적한다. meridional fan은 STOP T축, sagittal fan은 S축을 통과한다. 각 ray의 실제 IMA 좌표를 chief 기준으로 투영한다.

$$
\Delta_T(q)=\mathbf e_T^T[\mathbf x_T(q)-\mathbf x_c],\qquad
\Delta_S(q)=\mathbf e_S^T[\mathbf x_S(q)-\mathbf x_c].
$$

Fan별 두 transverse 성분을 모두 보존하면 cross-component와 2D shape도 검사할 수 있다. 원래 pupil 좌표, clipping, wavelength, common IMA를 기록한다. `chief reference`와 centroid reference를 자동 교체하지 않는다.

축상 meridional ray의 종방향 축 교차 위치는

$$
z_{cross}(q)=z-y\frac{d_z}{d_y},\qquad
\mathrm{LSA}(q)=z_{cross}(q)-z_{paraxial,infinity/finite}.
$$

출력 공액에 맞는 paraxial reference를 선택한다. q=0에서 0/0은 parabasal 극한으로 처리하며 $d_y=0,y\ne0$은 finite crossing이 없다. 이 LSA는 원래 ray의 정확한 기하학적 교차 계산이며 모든 비축 ray에 적용해 구면수차라고 부르지 않는다.

Coma는 일반적으로 단일 숫자가 아니다. 독립 진단량으로 meridional fan의 짝/홀 성분을 분리할 수 있다.

$$
\Delta_{even}(q)=\frac{\Delta_T(q)+\Delta_T(-q)}2,\qquad
\Delta_{odd}(q)=\frac{\Delta_T(q)-\Delta_T(-q)}2.
$$

저차 대칭 조건에서 짝 성분은 coma 계열, 홀 성분은 defocus·spherical 계열을 분리하는 데 유용하다. 그러나 비축 고차항을 포함하므로 이를 vendor Seidel coma coefficient라고 이름 붙이지 않는다. 세부 OPD 수차 descriptor는 다음 G9.1의 Zernike basis fit 및 잔차로 계산한다. Seidel은 저차 근사 지표이고 A₆…A₁₂ 전체 영향을 대체하지 않는다. [S13](#ref-S13)[S14](#ref-S14)

Astigmatism은 아래 parabasal T/S focus separation으로 정량화할 수 있다. 색수차는 고정 physical field·STOP·IMA에서 파장별 chief 또는 centroid map 차이(횡색수차)와 같은 정의의 초점 차이(종색수차)로 계산한다. 파장별 image를 따로 recenter/refocus한 뒤 색수차라고 보고하지 않는다.


#### G9.1 OPD의 Zernike descriptor: 수차 계수를 실제로 계산하는 규약

이 출력은 성능 계산에 사용한 raw OPD의 **진단용 분해**다. fixed IMA의 MTF에 들어가는 W를 이 fit으로 교체하거나 자동 defocus 제거하지 않는다. 입력은 동일 field·파장·온도·reference sphere/plane·OPD 부호로 얻은 $W_i$(절대 굴절률로 계산한 OPD mm), 유효 ray의 STOP 교점, pupil mask, quadrature measure다. 파면 모듈에서 이미 reference된 W를 ray ID로 받아야 하며 서로 다른 raw IMA OPL을 단순 차감해 만들어서는 안 된다.

**선언한 좌표.** 기본 descriptor chart는 물리 STOP의 고정 원점·반경 a_stop에 대한 정규화 좌표다. $(u,v)=(x_{stop}/a_{stop},y_{stop}/a_{stop})$, $\rho=\sqrt{u^2+v^2}$, $\varphi=\operatorname{atan2}(v,u)$로 둔다. u는 object/global X, v는 Y이고 field에 따라 자동 회전하거나 통과 pupil의 bounding box에 맞춰 다시 축척하지 않는다. $\rho\le1$인 원판의 W는 해당 STOP ray를 따라 파면 기준 surface로 옮겨 얻은 pullback이다. 따라서 이 descriptor를 실제 exit pupil 좌표의 vendor Zernike와 같은 값이라고 부르지 않는다. 다른 chart를 사용한다면 `chart_id, center, scale, basis, mapping`을 별도 profile로 저장한다. 다중 ray branch가 같은 STOP 좌표에 존재하면 기본 single-valued fit은 미지원이다.

**Basis와 index.** 단일 Noll 번호 대신 tuple $(n,m,t)$를 사용한다. n≥0, 0≤m≤n, n−m은 짝수다. m=0이면 t=0만, m>0이면 t=c와 s를 각각 사용한다. [S14](#ref-S14)의 일반 radial polynomial을 m≥0에 대해 사용하고 sin/cos를 명시적으로 분리한다.

$$
R_n^m(\rho)=\sum_{j=0}^{(n-m)/2}(-1)^j
\frac{(n-j)!}{j![(n+m)/2-j]![(n-m)/2-j]!}\rho^{n-2j}.
$$

$$
Z_{n,0,0}=\sqrt{n+1}\,R_n^0(\rho),\qquad
Z_{n,m,c}=\sqrt{2(n+1)}\,R_n^m(\rho)\cos(m\varphi),\qquad
Z_{n,m,s}=\sqrt{2(n+1)}\,R_n^m(\rho)\sin(m\varphi)\quad(m>0).
$$

이는 전체 단위원판의 균일 면적 측도에서

$$
\frac1\pi\int_{u^2+v^2\le1} Z_j Z_k\,du\,dv=\delta_{jk}
$$

가 되는 정규화다. n 순서, 각 n에서 m 오름차순, m>0에서 c 다음 s 순서로 직렬화하며 tuple을 항상 함께 저장한다. source의 signed-m/index 예시를 변환 없이 복사하지 않는다. factorial 식은 낮은 차수에서 정수 계수로 생성하고 Horner evaluation을 사용할 수 있다. 큰 n에서는 recurrence 또는 고정밀 계수로 cancellation을 확인한다. 실제 sag가 A₁₂에서 끝난다고 OPD가 Zernike n=12에서 끝나는 것은 아니다. `n_max`는 분석 profile에 명시하고 증가시킨 결과의 잔차·낮은 차수 계수 안정성을 확인한다.

| Tuple | 명칭 | 이 규약의 polynomial |
|---|---|---|
| (0,0,0) | piston | 1 |
| (1,1,c), (1,1,s) | x/y tilt | 2u, 2v |
| (2,0,0) | defocus | √3(2ρ²−1) |
| (2,2,c), (2,2,s) | primary astigmatism | √6 ρ² cos2φ, √6 ρ² sin2φ |
| (3,1,c), (3,1,s) | primary x/y coma | √8(3ρ³−2ρ)cosφ, √8(3ρ³−2ρ)sinφ |
| (3,3,c), (3,3,s) | trefoil | √8 ρ³ cos3φ, √8 ρ³ sin3φ |
| (4,0,0) | primary spherical | √5(6ρ⁴−6ρ²+1) |
| (5,1,c), (5,1,s) | higher-order coma descriptor | 위 일반식의 (5,1) |
| (6,0,0) | secondary spherical descriptor | √7(20ρ⁶−30ρ⁴+12ρ²−1) |

이 명칭은 선언된 wavefront basis 성분이다. Seidel coefficient, longitudinal aberration, vendor Zernike index/normalization과 동일하지 않다. coma/astig 두 성분을 합성 크기로 제공한다면 $a_{pair}=\sqrt{a_c^2+a_s^2}$와 원래 두 signed 계수를 함께 저장한다.

**가중 최소제곱.** 기본 fit은 통과한 STOP 영역의 균일 면적 가중이다. 샘플 셀 면적 또는 quadrature weight가 $\omega_i>0$다. Polar sampling의 측도는 $\rho\,d\rho\,d\varphi$이므로 모든 ring의 weight를 같게 두지 않는다. exit/cosine grid의 샘플을 재사용하면 STOP chart로의 면적 Jacobian을 반영하거나, 별도의 uniform STOP quadrature를 추적한다. 진폭·전력 가중 fit을 원하면 명시된 다른 profile로 취급한다.

$$
\widehat{\mathbf a}=\underset{\mathbf a}{\operatorname{argmin}}
\sum_{i\in valid}\omega_i\left[W_i-\sum_{j=1}^{K}a_jZ_j(u_i,v_i)\right]^2.
$$

실제로 $A_{ij}=\sqrt{\omega_i}\,Z_j(u_i,v_i)$, $b_i=\sqrt{\omega_i}\,W_i$를 구성하여 pivoted QR 또는 SVD로 푼다. 정규방정식 AᵀA의 역행렬을 기본 구현으로 사용하지 않는다. 샘플 수가 K 이상인 것은 필요조건일 뿐 rank를 보장하지 않는다. SVD의 singular values, 선택된 rank threshold, effective rank, condition number를 저장한다. rank가 부족하면 계수의 고유한 식별이 불가능하므로 `fit_rank_deficient`로 반환한다. truncated-SVD 최소 norm 값을 진단으로 남길 수 있지만 unique aberration coefficient라고 보고하지 않는다. clipping으로 빠진 ray를 W=0인 샘플로 채우지 않는다.

$$
e_i=W_i-\sum_ja_jZ_j(u_i,v_i),\qquad
RMS_{res}=\sqrt{\frac{\sum_i\omega_i e_i^2}{\sum_i\omega_i}},\qquad
PV_{res,sampled}=\max_i e_i-\min_i e_i.
$$

전체 원판이 아니거나 weight가 바뀌면 위 Z_j는 통과 영역에서 직교하지 않다. 따라서 계수 제곱합을 곧바로 측정 wavefront RMS라고 사용하지 않고 **실제 샘플에서 weighted RMS와 잔차를 다시 계산**한다. 독립 평가 샘플에서 잔차를 확인하면 overfit/누락된 pupil 구조를 검사할 수 있다. sampled PV는 연속 pupil의 보장된 extrema가 아니다.

**단위와 제거 정책.** a_j와 RMS_res는 OPD mm다. waves 단위는 같은 vacuum λ₀로 $a_j/\lambda_0$, radians는 $2\pi a_j/\lambda_0$다. 기본 fit은 piston·tilt·defocus를 모두 포함해 계수를 보고한다. 보고용으로 제거할 tuple 집합 R을 지정하면 immutable raw W에서 별도 배열

$$
W_{diagnostic,i}=W_i-\sum_{j\in R}\widehat a_j Z_j(u_i,v_i)
$$

를 만들고 `removed_terms`, 기준 좌표, reference geometry, aperture/weights, n_max를 기록한다. 비직교 clipped pupil에서 이 residual은 해당 selected modes만 따로 best-fit한 결과와 다를 수 있으므로 `joint_full_basis_subtraction`으로 이름 붙인다. nuisance terms만 최소제곱 제거하는 분석은 그 작은 basis를 별도로 fit하여 다른 label로 저장한다. **어느 diagnostic subtraction도 fixed-IMA FFT에 쓰는 W, PSF registration, 실제 refocus 위치를 변경하지 않는다.** Defocus 계수를 0으로 만드는 수치 조작을 실제 sensor 이동으로 해석하지 않는다.

Source QA 메모: S14의 일반 radial 식과 원판 정규화는 근거로 사용하되 웹 페이지의 단일-index 표·일부 예시 표기는 그대로 이식하지 않는다. 예를 들어 coma의 radial polynomial은 위 일반식에서 직접 R₃¹=3ρ³−2ρ로 산출한다. 본 문서의 tuple 방식은 index와 signed-m 혼동을 제거한 독립 규약이다.

### G10. Parabasal T/S field curvature

제품 Field Curvature는 field별 chief 주변 parabasal X/Y ray의 초점 z와 시스템 IMA z의 차이다. full-pupil RMS spot 최소점과 다른 quantity다. Y scan에서는 T=YZ,S=XZ이며 X scan에서는 반대다. 기본 적용 전제는 동축 회전대칭계와 평면 OBJ/IMA다. [S11](#ref-S11)

실제 구현은 STOP pupil parameter $\mathbf a=(a_T,a_S)$에 대해 chief 주변 ray family를 미분한다. IMA 또는 같은 최종 균질 매질의 기준 z₀ 평면에서 transverse 위치 $\mathbf X(\mathbf a)$, slope $\mathbf b(\mathbf a)=(d_x/d_z,d_y/d_z)$를 얻는다.

$$
J_X=\left.\frac{\partial\mathbf X}{\partial\mathbf a}\right|_{0},\qquad
J_b=\left.\frac{\partial\mathbf b}{\partial\mathbf a}\right|_{0},\qquad
\partial_\xi b_x=\frac{d_z\partial_\xi d_x-d_x\partial_\xi d_z}{d_z^2}.
$$

G5의 launch Jacobian을 STOP-coordinate Jacobian으로 바꾸려면 aiming Jacobian $J_{stop}$을 사용해 $d\xi/da=J_{stop}^{-1}$로 연결한다. 또는 물점을 고정한 채 STOP 목표 a±εe_T/S를 실제 재aim하여 central difference를 구한다. finite difference는 ε,ε/2,ε/4 및 다른 step 계열에서 수렴을 확인한다. clipping된 이웃 ray를 숫자 0으로 넣지 않는다.

Field의 T/S basis를 위치와 pupil 좌표 양쪽에 적용하여 $J_{X,TS}$, $J_{b,TS}$를 만든다. 회전대칭 meridional field에서 cross term이 수치 오차 이내라면

$$
\delta_T=-\frac{(J_{X,TS})_{TT}}{(J_{b,TS})_{TT}},\qquad
\delta_S=-\frac{(J_{X,TS})_{SS}}{(J_{b,TS})_{SS}},
$$

$$
z_T=z_0+\delta_T,\quad z_S=z_0+\delta_S,\quad
FC_T=z_T-z_{IMA},\quad FC_S=z_S-z_{IMA},\quad
\Delta z_{ast}=z_T-z_S.
$$

유도는 $\Delta\mathbf X(z)=J_X\Delta\mathbf a+(z-z_0)J_b\Delta\mathbf a$의 국소 초점 조건이다. 분모가 0이면 유한 T/S focus 없음이며 작은 분모는 ill-conditioned다. 축상은 T와 S가 같아야 하지만 IMA가 paraxial focus와 다르면 FC 값은 0이 아니다.

교차 coupling이 유의하면 단순 대각비로 값을 강제하지 않는다. 일반화된 line-focus는 $\det(J_X+\delta J_b)=0$으로 탐색할 수 있으나 이 결과는 별도 `generalized_parabasal_focus`이고 기존 T/S 제품 label이 아니다. 기본 동축 범위에서는 coupling이 나타난 원인(좌표·aiming·처방 위반)을 먼저 조사한다.

Afocal에서는 이 공간 z를 무한히 큰 숫자로 출력하지 않는다. 독립 출력으로 기준면의 converging-positive local vergence를 다음처럼 정의한다(δ는 mm, V는 m⁻¹=diopter).

$$
V_{T,ind}=1000n_{i,abs}/\delta_T,\qquad V_{S,ind}=1000n_{i,abs}/\delta_S.
$$

이는 광선의 downstream convergence를 양수로 정한 독립 paraxial local 진단이다. δ가 무한대인 regular bundle의 V는 0이고, δ≈0인 caustic의 V는 singular다. reference z₀를 결과에 반드시 남긴다. Zemax afocal diopter의 완전한 부호·reference 호환은 미확정이다. **본 규약의 주 확정 결과는 angular spot/MTF와 출력 direction map이며, vendor diopter 동등성은 미검증 항목으로 남긴다.**

### G11. Distortion, grid mapping 및 local frequency 변환

실제 chief map을 $\mathbf g(\mathbf q)$라 한다. focal은 IMA의 실제 chief 위치, afocal은 명시 angular chart의 chief 방향이다. ROI 크기만으로 무한 물체 angular field를 추정하지 않는다. 유한 12.8×3.2 mm source를 사용하는 경우 물체 좌표 q의 범위는 명확하지만 출력 FOV는 trace 결과다.

일반적인 2D 기준은 지정 reference field q₀에서 실제 small-field chief trace로 Jacobian을 얻는 것이다.

$$
J_0=\left.\frac{\partial\mathbf g}{\partial\boldsymbol\chi}\right|_{\chi_0},\quad
\mathbf g_{ref}(\chi)=\mathbf g(\chi_0)+J_0(\chi-\chi_0),\quad
\Delta\mathbf g=\mathbf g-\mathbf g_{ref}.
$$

여기서 χ는 finite-object height, rectilinear infinite-field의 (tanθx,tanθy), 또는 명시 F-theta angular coordinate다. 작은 field ε의 central difference와 step convergence로 J₀를 구한다. 제품도 실제 small-field ray를 기준으로 왜곡을 일반화하므로 defocused IMA에서 무조건 EFL×tanθ를 참조로 쓰지 않는다. [S11](#ref-S11)[S12](#ref-S12)

회전대칭의 한 meridian에서 signed 실제/기준 높이를 h,h_ref라 하면

$$
D(\chi)=100\frac{h-h_{ref}}{h_{ref}}.
$$

광축에서 비율 0/0을 직접 계산하지 않는다. uncalibrated local derivative 기준에서는 그 극한을 사용한다. calibrated best-fit scale을 별도 지원하면 축상 극한이 0이 아닐 수 있다. 2D 비대칭/불안정 branch에서는 임의 하나의 %보다 $\Delta g_x,\Delta g_y$를 보고한다.

Grid는 실제 물체 점들을 trace한 chief 위치의 연결이다. Grid별 spot/PSF는 같은 IMA에 overlay하고, 개별 focus map은 별도 결과로 제공한다. 이 세 결과를 구별한다. Rectangle field는 중심·장변 edge·단변 edge·corner·중간 field를 포함하고 clipping/강한 비선형 구간에서 adaptive refinement한다. 회전대칭 PSF의 재사용을 할 때에도 전체 2D kernel을 회전해야 하며 T/S 두 곡선만으로 이를 복원하지 않는다.

온도/설계 비교에는 각 조건에서 다시 fitting한 distortion뿐 아니라 **nominal 기준에 대한 raw map 차이**를 남긴다.

$$
\Delta\mathbf g_{nom}(\mathbf q,T)=\mathbf g(\mathbf q,T)-\mathbf g(\mathbf q,T_0).
$$

이 값과 local scale/Jacobian, FOV를 함께 저장해야 전체 배율 변화가 refit에 의해 사라지지 않는다. 예를 들어 $g_T=1.01g_0$인 선형 map은 온도별 distortion이 모두 0이어도 영상 크기는 1% 달라진다.

국소 map $\mathbf x=J\mathbf q$에서 물체/영상 주파수는 위상 보존으로

$$
\mathbf f_q=J^T\mathbf f_x,\qquad \mathbf f_x=J^{-T}\mathbf f_q.
$$

Angular map에도 동일하되 단위가 rad/mm이면 그 역전치 단위까지 적용한다. EFL 하나로 모든 비축 주파수를 바꾸거나 lp/mm를 같은 숫자의 cycles/mrad로 바꾸지 않는다. 이 좌표관계가 forward/reverse PSF·RI의 무조건적 동일성을 뜻하지는 않는다.

### G12. 공통 초점면과 through-focus 계산

최종 균질 매질의 평면 이동에서는 ray를 다시 굴절시킬 필요 없이 기준면 위치 aᵢ와 slope bᵢ로

$$
\mathbf x_i(\delta)=\mathbf a_i+\delta\mathbf b_i
$$

를 사용할 수 있다. 유효 ray 집합·weights·sensor clipping이 δ에 따라 바뀌지 않는 조건에서 centroid RMS²의 최소는

$$
\delta_{RMS}=-\frac{\sum_iw_i(\mathbf a_i-\bar{\mathbf a})\cdot(\mathbf b_i-\bar{\mathbf b})}
{\sum_iw_i\|\mathbf b_i-\bar{\mathbf b}\|^2}.
$$

공통 sensor 초점을 여러 field에 맞추려면 field weight v_h와 field별 정규화 moments를 먼저 정의한 뒤 분자/분모를 field에 걸쳐 합한다.

$$
A_h=\frac{\sum_iw_{hi}(\mathbf a_{hi}-\bar{\mathbf a}_h)\cdot(\mathbf b_{hi}-\bar{\mathbf b}_h)}{\sum_iw_{hi}},\quad
B_h=\frac{\sum_iw_{hi}\|\mathbf b_{hi}-\bar{\mathbf b}_h\|^2}{\sum_iw_{hi}},\qquad
\delta_{common,RMS}=-\frac{\sum_hv_hA_h}{\sum_hv_hB_h}.
$$

분모가 0이면 유일한 finite optimum이 없다. 이 식은 geometric RMS 기준이고 최대 MTF focus, wavefront RMS 최소, parabasal curvature와 같지 않다. sensor aperture·wavelength throughput·ray 집합이 달라지면 위 closed form을 사용하지 않고 각 δ에서 목적함수를 재평가한다.

MTF 기준 refocus는 허용 이동 범위를 sampling하여 동일 ordinary FFT API를 각 δ에서 실행하고, 선언한 field/frequency/T/S 가중 merit를 1차원 bounded search로 개선한다. IMA만 이동하면 회절 장의 canonical Π와 그 위 A를 유지하는 `hold_boundary_field_amplitude`를 적용하고, 변경된 reference로 B를 다시 계산한다. 매번 B=1로 재정의하는 대안과 혼합하지 않는다. MTF merit가 비단봉일 수 있으므로 최초 local optimizer 한 번으로 global best라 하지 않는다. 모든 field·주파수에 공통 δ 하나를 적용한다. field별 최댓값 envelope는 진단용이며 평면 센서 결과가 아니다. 제품 FFT Through Focus의 별도 빠른 계산 경로와 ordinary FFT 반복 sweep의 동등성은 따로 검증해야 한다. [G02](#ref-G02)

Afocal angular 성능은 IMA 평행이동으로 변하지 않는다. refocus를 원하면 OBJ나 허용 lens group의 한 명시 자유도를 이동시키고 전체 trace를 다시 수행하거나 refocus를 `not_applicable`로 반환한다. 고정 시스템 결과와 보상 후 결과, 필요한 보상량을 함께 보존한다.

### G13. 반환 상태, 검증 및 남는 쟁점

| 상태 | 의미 및 출력 처리 |
|---|---|
| `valid` | 선언한 모델/수치 조건에서 해당 quantity 계산됨 |
| `aperture_blocked`, `escaped`, `TIR` | 물리적인 직접 경로 손실. power 원장에 반영; solver 실패와 구분 |
| `intersection_unresolved`, `aim_unresolved`, `critical_angle_unresolved` | 수치적으로 미결정. 차단/0 성능으로 바꾸지 않음 |
| `invalid_geometry` | 정의역 오류, 면 교차, 재료/순서 모순 등 |
| `grazing_degenerate`, `multi_branch_aim` | 기본 매끄러운 단일 branch 모델의 지원 경계 |
| `ill_conditioned` | EFL/초점/Jacobian 등의 작은 분모로 안정성 부족 |
| `not_applicable`, `unsupported` | 정의 부적합 또는 구현하지 않은 convention; 가능한 다른 outputs까지 버리지 않음 |

검증은 실제 생산 API에 연결한다: plane/sphere와 conic branch; A₂≠0의 vertex curvature; Snell의 접선 운동량과 투과 반공간; 평행판의 OPL; small-pupil real-ray limit와 독립 ABCD; forward ray를 하나씩 정확히 뒤집은 path/OPL; STOP residual; Jacobian과 finite difference; 축상 T=S; common RMS closed form과 plane sweep. 파면 모듈과는 공통 평면의 매끄러운 단일 branch에서 $\nabla_\perp S=n\mathbf d_\perp$가 연결되는지 확인한다. 반대쪽에서 새 plane-wave bundle을 만든 것을 정확한 역추적이라고 부르지 않는다.

남는 쟁점은 숨기지 않는다.

- 제품의 모든 field/aperture/vignetting/working-F 예외를 모방하는 범위는 별도 compatibility profile이다. 여기의 지정 알고리즘은 독립 계산 기준이며 버전별 수치 일치는 미검증이다.
- A₂가 있는 asphere의 실제 작은 pupil 한계는 c+2A₂를 반드시 반영한다. 특정 vendor 보고 창이 모든 surface type에서 이 항을 같은 방식으로 처리하는지의 수치 확인은 별도 benchmark다.
- active 영역 밖 transparent rim을 원하면 높이 외에 실제 surface profile이 필요하다. opaque 기본 프로필은 그 물리 현상을 포함하지 않는다.
- 다중 ray-aiming 해·pupil fold·grazing은 수학적으로 존재할 수 있다. 단일 branch 해석에서 임의로 한 해를 숨겨 채택하지 않는다. 더 일반적인 처방 지원은 transport/wave 모듈과 함께 확장한다.
- Afocal vendor diopter의 구체적 부호/기준면, calibrated distortion의 제품 내부 fitting objective, 특수 virtual/intermediate surface 분석까지 수치적으로 같다고 선언하지 않는다. 일반 chief map·angular spot·독립 local vergence는 명시 규약으로 계산 가능하다.
- 모든 derivative·행렬·spot/초점 식은 본문 조건에서의 직접 유도다. 공식 자료의 기능 정의만으로 엔진의 실제 정확도나 vendor 대비 오차 수치가 증명되지는 않는다.

---

<a id="part-4"></a>

## 4. OPD·FFT PSF/MTF·회절 기준

이 절은 `independent_fft_profile_v1`의 계산 규약이다. 공개 scalar diffraction 이론에 따라 실행 가능한 좌표·측도·진폭·부호·주파수 변환을 정의한다. **Zemax의 비공개 cosine-space mapping을 역으로 알아냈다는 뜻이 아니다.** 제품 FFT와 같은 공개 계산 계열을 사용하지만 수치 동등성은 `vendor_cosine_space_equivalence=OPEN`이다. [S15](#ref-S15) [S40](#ref-S40) [S42](#ref-S42) [W01](#ref-W01)

일반 출력은 G의 lens-only incoherent FFT MTF다. 한 물점·한 파장의 pupil 내 복소장은 coherent하게 합산하고, PSF intensity의 Fourier 변환을 OTF로 정의한다. LED bitmap, 320×80 sampling, fill factor 및 pixel aperture MTF는 여기에 곱하지 않는다. 공간적으로 서로 비간섭인 LED 물점 및 서로 다른 파장은 이후 intensity로 합산한다.

### W.1 공통 부호·매질·경계장

시간 의존성은 `exp(−iωt)`, 진행파 위상은 `exp(+ik₀S)`로 고정한다. 길이 단위는 mm, 진공 파장 λ₀는 mm, absolute phase index는 n이다.

$$
k_0=\frac{2\pi}{\lambda_0},\qquad
\lambda_m=\frac{\lambda_0}{n_i},\qquad k_i=\frac{2\pi}{\lambda_m}.
$$

`n_i`는 마지막 굴절면 이후 출력 영역의 균질 굴절률이다. Group index를 사용하지 않는다. 분산/온도마다 n을 갱신한다. 굴절률을 이미 포함한 OPL에 다시 n을 곱하지 않는다.

광선마다 발사점 `a`, 공통 입력 위상의 eikonal `S_in(a)`, 각 실제 구간 길이 `ℓ_j`, 구간 n을 저장한다. 경계면까지의 위상 길이는 다음과 같다.

$$
S_\Pi=S_{\rm in}(\mathbf a)+\sum_j n_j\ell_j.
$$

- 유한 실물점 `o`에서 발사하면 첫 실제 구간을 `n_o|a−o|`로 포함한다. 이 항을 ray마다 0으로 두지 않는다.
- 무한 물체 평면파는 `S_in(a)=n_o d_o·(a−a_ref)`다. 발사점들이 같은 등위상면 위에 있지 않으면 동일 phase로 놓으면 안 된다.
- 지정된 수렴 입사 구면파/가상 물체를 지원할 때에는 `S_in(a)=C−n_o|a−o_v|`처럼 **파면과 진행 방향이 일치하는 부호**를 사용한다. 가상거리의 절댓값을 무조건 더하지 않는다.
- Piston은 임의 상수를 빼도 된다. Tilt와 defocus는 고정면 결과에서 임의 least-squares fit으로 제거하지 않는다.
- 복소 coating 위상은 기본 scalar intensity-transmission 모드에 포함하지 않는다. 지원 시 ray마다 phase를 별도 더하고, amplitude만으로 대체하지 않는다.

마지막 광선점 `p`와 진행 방향 `d`로부터 공통 평면 `Π: (r−r₀)·e₃=0`까지는

$$
t_\Pi=\frac{(\mathbf r_0-\mathbf p)\cdot\mathbf e_3}
{\mathbf d\cdot\mathbf e_3},\qquad
\mathbf r=\mathbf p+t_\Pi\mathbf d,\qquad
S_\Pi=S_p+n_i t_\Pi.
$$

OPL의 이 단계는 signed다. 기본 profile은 **모든 유효 ray가 마지막 물리면을 지난 후의 실제 균질 영역에 Π를 둔다**. 역연장 `t_Π<0`을 쓸 경우 같은 균질 구간에서만 허용하며, 렌즈·aperture를 가로지르는 가상 평면을 물리 경계장으로 취급하지 않는다. 실제 출사동공의 근축 위치가 렌즈 내부라고 해서 그 위치에 무조건 `Π`를 두어서는 안 된다.

`Π`의 각 점에는 단일한 매끄러운 phase branch가 있어야 한다. Ray map fold, caustic, grazing, 서로 다른 위상의 여러 branch, 통과하지 못하는 chief ray가 있으면 다른 regular plane을 선택하거나 wave metric에 `unsupported`를 반환한다. 이때 geometric spot은 별도로 유효할 수 있다. 다중 branch의 Maslov phase·vector 전파는 초기 profile의 범위 밖이다.

Π의 경계장은 기하광학 eikonal·수송 진폭으로 구성한다. 내부 STOP/렌즈 edge에서 생긴 회절을 다음 렌즈까지 단계마다 전파하는 다면 wave-optics 모델은 아니다. 마지막 적분을 정확한 scalar kernel로 계산해도 이 boundary 근사의 오차가 사라지는 것은 아니다.

Π의 선택도 재현 가능한 규약이어야 한다. 기본 후보는 유효한 마지막 surface 영역을 e₃ 방향으로 투영한 최대값보다 `δ_clear>0`만큼 downstream인 chief-normal 평면이다. δ_clear는 geometry tolerance를 넘는 계산 설정으로 저장하고 실제 경계 전체에 대해 검사한다. Sparse ray의 최대값만으로 모든 ray가 이미 마지막면을 지났다고 보장하지 않는다. Focal 모드에서는 이 후보가 q₀보다 앞에 있고, ray map이 regular한지 확인한다. 가능한 Π가 없으면 이 profile을 거부한다. 물리 source kernel을 다른 regular Π에서 다시 계산하는 시험은 boundary/transport 일관성 진단이지만, 임의 Π에서 B를 다시 균일하게 정하는 lens metric은 조명이 바뀔 수 있으므로 자동 불변성을 요구하지 않는다.

선택한 Π의 origin·basis·선택 기준·δ_clear·유효성 범위를 **canonical immutable profile state**에 넣는다. 보통의 IMA-only focus sweep은 이 Π와 ray-boundary 데이터를 유지하고 출력면과 reference만 바꾼다. Sweep 중 Π가 유효하지 않아 새 평면이 필요하면 기존 결과를 조용히 이어 붙이지 않는다. 새로운 `profile_instance_id`로 재계산하여 평면 변경을 기록하고, 동일 비교군 전체를 새 기준에서 다시 계산하거나 해당 sweep 점을 unsupported로 남긴다. 렌즈 이동·온도 변화로 geometry 자체가 달라지는 경우에도 Π와 illumination 정책의 변경 이력을 명시한다.

### W.2 Focal 출력면과 reference sphere

G chief ray가 사용자가 지정한 물리 IMA와 만나는 점을 `q₀`라 하자. 단색 계산에서는 해당 파장의 chief를 사용한다. `e₃`는 그 점에서의 chief 진행 방향이다. `e₁`은 OBJ local X 방향을 `e₃`에 수직인 평면에 투영해 정규화하고, `e₂=e₃×e₁`로 잡는다. 투영이 퇴화하면 고정된 대체 축과 상태를 기록한다.

출력 좌표와 boundary plane은

$$
\mathbf q(x,y)=\mathbf q_0+x\mathbf e_1+y\mathbf e_2,
\qquad
\mathbf r(u,v)=\mathbf r_0+u\mathbf e_1+v\mathbf e_2,
\qquad
\mathbf q_0=\mathbf r_0+L\mathbf e_3,\quad L>0.
$$

이 출력은 **chief-normal 가상 평면의 mm 좌표**다. IMA의 위치는 기준점을 결정하지만 IMA가 기울었다고 이 결과가 실제 센서면 irradiance로 바뀌지는 않는다. Zemax FFT PSF도 이 가상면 성격을 명시한다. [S42](#ref-S42)

실제 IMA로 변환하려면 별도 평면 전파를 사용하거나, 유효성이 확인된 작은 PSF에 한해서 국소 투영을 적용한다. IMA의 직교 기저를 `b₁,b₂`라 할 때 chief 방향의 국소 평행 투영은

$$
\begin{bmatrix}x\\y\end{bmatrix}
=J\begin{bmatrix}X\\Y\end{bmatrix},\qquad
J_{ab}=\mathbf e_a\cdot\mathbf b_b,\quad a,b\in\{1,2\},
\qquad
h_{\rm IMA}(X,Y)=h_\perp(J[X,Y]^T)|\det J|.
$$

이는 광선들이 국소적으로 chief와 거의 평행하고 tilted plane에 걸친 defocus 변화를 무시할 수 있는 **추가 근사**다. 큰 NA·큰 tilt에서 단순 `cos(chief angle)` 보정으로 실제 센서 PSF를 복구하지 않는다. 실제 IMA 결과와 FFT 지표를 다른 result key에 저장한다.

Reference sphere는 중심이 `q₀`인 이상 수렴파다. Π 위에서

$$
R(u,v)=\sqrt{L^2+u^2+v^2},\qquad
S_{\rm ref}(u,v)=S_c+n_i[L-R(u,v)],
$$
$$
W(u,v)=S_\Pi(u,v)-S_c+n_i[R(u,v)-L].
\tag{W-1}
$$

`S_c`는 chief의 Π 위 eikonal이다. 이상 수렴 구면파이면 W=0이다. 이는 **reference sphere까지 실제 ray를 이동시킨 raw OPL 차이**를 요구하지 않고, 실제 Π 위 경계 phase에서 이상 구면 phase를 빼는 정의다. `R−L=(u²+v²)/(R+L)`로 계산하면 원점 근처 cancellation을 줄인다.

고정 IMA의 defocus는 W 안에 남는다. IMA를 이동한 focus sweep은 `q₀,L,R,W`를 그 위치에 맞춰 다시 만든다. 각 field/주파수/T·S의 최적 W를 따로 선택해 하나의 센서 성능으로 합치지 않는다.

### W.3 왜 이 좌표의 FFT인가: 독립 scalar profile의 유도

평면 경계장에서 downstream 평면으로의 Rayleigh–Sommerfeld I 전파를 사용한다. 아래 첫 식은 scalar Helmholtz half-space Green function의 법선 미분을 풀어 쓴 것으로, 파면 경계값이 정확히 주어졌을 때의 scalar 식이다. 실제 렌즈의 full Maxwell 해라는 뜻은 아니다. `R_q=|q−r|`, `q−r`의 longitudinal 성분은 L이다.

$$
U(\mathbf q)=\iint_\Pi U_\Pi(u,v)
\frac{L}{2\pi R_q^3}(1-i k_iR_q)e^{i k_iR_q}\,du\,dv.
\tag{W-2}
$$

`k_iR_q≫1`에서는 작은 `1/(k_iR_q)` 항을 생략하여

$$
U(\mathbf q)\simeq\frac{1}{i\lambda_m}
\iint U_\Pi(u,v)\frac{L}{R_q^2}e^{i k_iR_q}\,du\,dv.
\tag{W-3}
$$

표준 RS/inverse-Fraunhofer 출발식은 [W01](#ref-W01) §§2–3이다. 이하 임의 ray-derived `W,A`를 연결하고 유효성 잔차를 정의한 부분은 이 문서의 구현용 유도다. **모든 sinθ를 θ로 치환하는 근축 광선추적은 하지 않는다.** 다만 다음의 **출력 영역에 대한 국소 전개**는 남는다.

$$
\chi=\frac{u}{R},\qquad\eta=\frac{v}{R},\qquad
\gamma=\sqrt{1-\chi^2-\eta^2}=\frac{L}{R},
$$
$$
R_q=R-\chi x-\eta y+
\frac{x^2+y^2-(\chi x+\eta y)^2}{2R}
+O\!\left(\frac{\|(x,y)\|^3}{R^2}\right).
$$

기준 profile은 exponential에서 1차까지, kernel amplitude에서는 0차까지 유지한다. 역변환과 면적 Jacobian은 정확하다.

$$
u=\frac{L\chi}{\gamma},\qquad v=\frac{L\eta}{\gamma},\qquad
du\,dv=\frac{L^2}{\gamma^4}\,d\chi\,d\eta.
$$

`U_Π=A_Π exp(ik₀S_Π)`를 넣으면, x·y에 무관한 위상·상수를 제외한 다음의 Fourier integral을 얻는다.

$$
B(\chi,\eta)=A_\Pi(u,v)\frac{L}{R^2}\frac{L^2}{\gamma^4}
=A_\Pi(u,v)\frac{L}{\gamma^2},
$$
$$
P(\chi,\eta)=M(\chi,\eta)B(\chi,\eta)e^{i k_0W(\chi,\eta)},
$$
$$
V(x,y)=\iint P(\chi,\eta)
\exp\!\left[-i2\pi\frac{x\chi+y\eta}{\lambda_m}\right]d\chi\,d\eta.
\tag{W-4}
$$

`M`은 모든 실제 clear aperture·STOP의 통과 마스크다. χ,η는 여기서 **이상 reference sphere를 향한 방향의 위치 성분**이며 수차가 있는 실제 ray의 방향 cosine와 무조건 같지 않다. 따라서 실제 `d_x,d_y`를 χ,η 자리에 그대로 대입하면 다른 알고리즘이 된다.

`A_Π∝f(θ,φ)/R`인 이상 수렴파에서는 `B∝f/γ`가 되어 [W01](#ref-W01) 식 (8)과 일치한다. 이 관계가 Jacobian 및 obliquity를 이중 적용하지 않았는지 확인하는 제한 검산이다. 이 차트에서의 균일 B와 평면 Π의 균일 전력밀도는 서로 다른 조명이다.

**적용성 잔차를 계산한다.** 실제 출력 ROI의 q 표본과 pupil r 표본에서

$$
\epsilon_{\phi}=k_i\max_{q,r}|R_q-R+\chi x+\eta y|,
\qquad
\epsilon_A=\max_{q,r}\left|\frac{R^2}{R_q^2}-1\right|.
\tag{W-5}
$$

이들은 국소 전파 근사 진단값이지 바로 MTF 오차 상한은 아니다. `k_i min(R_q)`와 scalar/vector 적용성도 따로 기록한다. ROI가 커지거나 횡수차가 매우 커서 잔차가 목적 오차에 비해 무시할 수 없으면 더 많은 FFT pixel만으로 해결되지 않는다. 일부 q에서 (W-2)의 직접 적분과 대조하고, 통과하지 못한 경우 `model_validity=unresolved/unsupported`로 둔다. 그 값을 Huygens MTF로 조용히 바꾸어 FFT label로 내보내지 않는다.

### W.4 진폭을 정하는 두 가지 서로 다른 계약

**(a) `lens_metric_pupil` — 일반 렌즈 FFT 지표.** 기본은 위 Fourier chart에서 `B=1`, physical clipping의 `M`, polarization off다. 지정 apodization을 `a(χ,η)`라 하면 `B=a`이고, power transmission τ를 포함하는 옵션에서는 `B=a√τ`다. τ를 다시 PSF intensity에 곱하지 않는다. 이 기본은 **독립 규약**이며 Zemax가 임의 처방·옵션에서 같은 B를 만든다는 확인은 남아 있다.

같은 physical illumination으로 through-focus를 비교하는 경우에는 reference focus의 B에서 `A_Π=Bγ²/L`를 한 번 만들고 Π 위 A를 고정한다. 다음 IMA 위치의 B는 바뀐 `L,γ`를 이용해 (W-4)로 다시 구한다. 반면 매 IMA에서 B=1을 다시 정의하는 `uniform_fourier_pupil_each_focus`는 명시적인 지표 convention이며, 고정된 물리 조명의 sweep과 같다고 선언하지 않는다. 기본 공통 through-focus의 illumination 정책은 `hold_boundary_field_amplitude`로 저장한다. Geometry와 τ가 변하는 실제 compensator/열 sweep에서는 해당 transport를 다시 계산한다.

**(b) `source_transport_pupil` — 물리 LED 영상/각도 kernel.** 각 source cell에서 전송한 전력을 physical Π의 면적으로 나눈 `p_Π=dΦ/dA_Π`를 먼저 구한다. Ray power 하나의 제곱근을 grid field에 곧바로 대입하지 않는다. 비자성 균질 매질의 국소 평면파 convention에서

$$
p_\Pi=\frac{n_i}{2Z_0}|A_{\Pi,E}|^2(\mathbf d\cdot\mathbf e_3),
\qquad
A_{\Pi,E}=\sqrt{\frac{2Z_0p_\Pi}{n_i(\mathbf d\cdot\mathbf e_3)}}.
\tag{W-6}
$$

전기장 V/m를 실제 사용할 때 `p_Π`는 W/m²로 변환한다 (`1 W/mm²=10⁶ W/m²`). 정규화 PSF만 필요하면 모든 곳에서 같은 상수 `2Z₀`와 단위 변환을 생략한 scaled amplitude를 사용할 수 있다. `n_i`와 방향 cosine가 공간/파장마다 달라지는 부분까지 임의 상수로 지우면 안 된다. 여기의 cosine는 **실제 ray**의 cosine다. (W-4)의 γ는 reference geometry이므로 두 요인을 혼동하지 않는다.

`p_Π`에 이미 surface/bulk τ가 들어갔다면 (W-6) 뒤에 `√τ`를 재적용하지 않는다. 이 `A_Π`를 (W-4)의 B로 변환한다. Lens metric PSF를 총 Q만 곱해서 source kernel로 쓰는 것은, 두 B의 모양이 상수배로 같다는 검증이 있을 때만 허용한다. Kernel의 총 직접 투과량 Q는 radiometry 모듈이 따로 전달하고, h는 아래처럼 단위면적/입체각에 대해 정규화한다.

### W.5 Uniform grid, aiming, DFT와 물리 단위

χ,η의 균일 격자를 만든 후 (W-3)의 역변환으로 목표 (u,v)를 구한다. 각 점에서 launch의 두 자유도를 조정하여 **실제 6면을 지난 ray가 목표 Π 점에 도달하도록** inverse aiming한다. 물점/입사 평면파 조건과 모든 aperture는 유지한다. Pupil 바깥의 physical miss는 M=0이지만 수치 aiming 실패는 masked zero가 아니다. Regular map에 대한 continuation/Newton 및 branch 확인을 사용한다. 비균일 traced grid로부터 보간한다면 unwrapped W와 amplitude를 separately 보간하고 support 경계를 보존한다. Wrapped phase angle을 선형 보간하지 않는다.

Support는 유효 STOP/clear-aperture 영역의 forward map 및 경계 refinement로 먼저 확인하거나, 물리 aperture 판정을 기록하는 독립 geometric continuation으로 얻는다. “Newton이 실패했다”는 사실만으로 바깥 pupil이라고 판정하지 않는다. Fractional edge cell은 복소 적분을 subcell quadrature하는 방식으로 처리한다. 셀의 통과 면적비를 복소 amplitude에 적용하는 것과 그 제곱근을 곱하는 것은 서로 다른 모델이다.

짝수 `Nχ,Nη`에 정수 index `j=−Nχ/2,…,Nχ/2−1`, `k=−Nη/2,…,Nη/2−1`를 사용하고 `χ_j=jΔχ`, `η_k=kΔη`로 고정한다. Chief chart 원점은 배열의 zero index에 정확히 들어간다. χ²+η²≥1은 사용할 수 없다. 배열 physical support와 zero padding 영역을 구분한다.

$$
V_{ab}=\Delta\chi\Delta\eta\sum_{j,k}P_{jk}
\exp\left[-i2\pi\left(\frac{aj}{N_\chi}+\frac{bk}{N_\eta}\right)\right],
$$
$$
x_a=\lambda_m\frac{a}{N_\chi\Delta\chi},\qquad
y_b=\lambda_m\frac{b}{N_\eta\Delta\eta}.
\tag{W-7}
$$

구현은 명시한 indexing에 맞춰 `fftshift(fft2(ifftshift(P))) ΔχΔη`다. 이는 NumPy의 forward-negative-exponent convention을 가정한 수학적 정의이며 특정 library 버전 의존 동작을 숨기지 않는다. `x,y,Δx,Δy`는 mm다.

$$
I_{ab}=|V_{ab}|^2,\qquad
h_{ab}=\frac{I_{ab}}{\sum_{a,b}I_{ab}\Delta x\Delta y},
\qquad \sum h_{ab}\Delta x\Delta y=1.
$$

분모 0은 null/실패 상태다. 이는 **전체 계산창 내 단위 적분**이다. 출력 display crop을 다시 1로 정규화하면 tail을 숨긴다. 계산창 확장 시 총 적분·central PSF·MTF 수렴, boundary strip의 에너지, 원래 crop의 captured fraction을 기록한다. 충분히 큰 창에서 검증된 h만 Q와 곱하여 물리 영상을 합성한다. Peak=1 그림은 별도 display normalization이다.

Forward FFT의 부호를 OTF에도 그대로 사용한다.

$$
\mathrm{OTF}(f_x,f_y)=\iint h(x,y)e^{-i2\pi(f_xx+f_yy)}dx\,dy,
\qquad \mathrm{MTF}=|\mathrm{OTF}|,
$$
$$
f_{x,c}=\frac{c}{N_\chi\Delta x},\qquad
f_{y,d}=\frac{d}{N_\eta\Delta y},\qquad
\mathrm{OTF}(0,0)=1.
\tag{W-8}
$$

MTF 주파수는 image-space cycles/mm = lp/mm다. 요청한 0–120은 결과를 읽는 대역이며 계산 pupil/PSF의 충분성을 정하지 않는다. Complex OTF를 보존하고 MTF를 임의 단조화하지 않는다. 주파수 slice의 독립 적분 또는 pupil autocorrelation으로 일부 f를 검산한다. 이 부호에서는

$$
\mathrm{OTF}(\mathbf f)=
\frac{\int P(\mathbf s)P^*(\mathbf s+\lambda_m\mathbf f)d^2s}
{\int|P(\mathbf s)|^2d^2s}.
\tag{W-9}
$$

`h(x−x₀)`의 OTF가 `exp(−i2π f·x₀)OTF(h)`가 되는지 검사하면 부호·등록을 확인할 수 있다. Autocorrelation을 library 함수로 계산할 때는 인수 순서에 따라 complex conjugation이 달라지므로 magnitude만 확인하지 않는다.

T/S 표시는 본 profile에서 T=f_y slice, S=f_x slice로 정의하고 해당 basis를 기록한다. Zemax FFT의 명칭은 object X에 평행한 줄무늬의 modulation이 T, object Y에 평행한 줄무늬가 S다. [S40](#ref-S40) 임의 2D field에서 이 basis와 제품의 entrance-pupil orientation mapping이 완전히 같은지는 별도 adapter 검증 항목이다. 자동 radial/azimuthal 회전을 제품 정의로 선언하지 않는다. 2D PSF/OTF를 보존하면 이후 적절한 slice를 다시 평가할 수 있다.

### W.6 Sampling을 네 가지 오차로 나눈다

1. **실제 phase/edge sampling:** 인접 unwrapped phase 차이 `|Δ(k₀W)|<π`는 유용한 필요 진단일 뿐 충분조건이 아니다. 중간점/offset grid로 급격한 위상 및 clipping 경계를 확인한다. Pupil의 세밀한 구조가 coarse grid에서 모두 빠질 수 있다.
2. **Zero padding:** 고정 Δχ,Δη에서 N을 늘리면 PSF pixel pitch가 작아지지만 실제 ray phase 정보는 늘지 않는다. OTF의 circular wrap을 막으려면 pupil의 linear autocorrelation support가 padded domain 안에 들어가도록 한다. 통상 각 축 비영점 support 폭의 약 2배 이상이 시작점이며 실제 support를 검사한다.
3. **계산창:** 작은 Δχ는 `X_period=λ_m/Δχ`를 키운다. Padding만 증가시키면 이 주기/물리창 크기는 그대로다. Tail truncation·주기 replica와 displayed crop을 구분한다.
4. **보간/출력 주파수:** FFT native frequency spacing이 요구한 f를 충분히 표현하는지 확인한다. 필요하면 (W-9)를 selected frequency에서 평가한다. 보간한 MTF만 보지 말고 complex OTF를 보간한다.

단색 PSF의 Nyquist는 `1/(2Δx)`다. 120을 포함하는 것 외에도 pupil autocorrelation의 유효 cutoff를 alias 없이 수용해야 한다. 두 배 grid의 반복 일치만으로 수렴을 증명하지 않는다. Non-nested/offset grid와 selected-point direct integration을 추가한다. 제품의 “64×64 sampling”과 독립 profile의 같은 N이 동일한 effective rays·padding·창을 뜻하지 않는다. 제품 sampling 관련 공개 설명은 [W02](#ref-W02)에 있다.

### W.7 Diffraction limit

**동일 pupil의 zero-phase 기준:** 해당 B와 M을 그대로 두고 W=0으로 만든 `P_DL=MB`를 같은 FFT·정규화로 계산한다. 이는 위 정의에 따른 무수차 기준이다. 실제 pupil이 field마다 다르면 field별 DL도 함께 저장할 수 있다. Zemax 표시선의 reference field 정책과는 구별한다. 제품의 통상 reference는 (0,0)이 있으면 그 field, 없으면 첫 field다. [S41](#ref-S41)

**균일 원형 Fourier pupil의 해석 검산:** χη chart에서 반경 a의 원형이고 B가 일정하면

$$
f_c=\frac{2a}{\lambda_m}=\frac{2\mathrm{NA}_{\rm ref}}{\lambda_0},
\qquad
\mathrm{NA}_{\rm ref}=n_i a,
$$
$$
\mathrm{MTF}_{\rm circ}(q)=\frac{2}{\pi}
\left[\arccos q-q\sqrt{1-q^2}\right],\quad q=f/f_c\in[0,1],
\qquad \mathrm{MTF}=0\ (q>1).
\tag{W-10}
$$

이 식은 원판의 정규화 overlap area에서 직접 유도한다. `a`는 reference chart의 실제 support 반경이다. 심하게 수차가 있는 계에서 이를 실제 marginal-ray NA 또는 vendor working F/#와 무조건 같다고 놓지 않는다. 근축 원형 공기계에서는 `f_c≈1/(λ₀F_work)`로 연결되지만, nominal EFL/EPD를 유한 공액/비축 모든 경우의 정확한 scaling으로 쓰지 않는다. 차폐·비네팅·apodization pupil에 원형 공식을 적용하지 않는다. DL을 pupil이 다른 모든 field의 절대 상한이라고 하지 않는다.

### W.8 Afocal: plane OPD와 angular FFT

Afocal은 출력 plane wave를 reference로 삼는다. 명목 collimator의 finite object도 이 branch에 포함하며, 렌즈 시스템 행렬 C=0을 요구하지 않는다. [S39](#ref-S39) Physical plane Π와 single branch 조건은 동일하다. 좌표 e₃는 기준 G chief 방향이며, monochromatic 개별 계산이면 해당 chief 방향을 쓸 수 있다.

$$
W_\infty(u,v)=S_\Pi(u,v)-S_0-n_i\mathbf d_{\rm ref}\cdot[\mathbf r(u,v)-\mathbf r_0],
\qquad U_0=M A_\Pi e^{ik_0 W_\infty}.
$$

e₃=d_ref이면 Π 위 마지막 항은 0이다. 잔여 quadratic/defocus phase를 제거하면 안 된다. Physical IMA 평면을 같은 균질 구간에서 옮기는 것만으로 collimation을 보정했다고 판단하지 않는다.

Uniform physical-plane `(u,v)` grid에서 `F(ν_u,ν_v)=∫U₀ exp[−i2π(ν_u u+ν_v v)]du dv`를 계산한다. 출력 local direction cosine를 `(l,m,γ_o)`라 하면

$$
l=\lambda_m\nu_u,\qquad m=\lambda_m\nu_v,\qquad
\gamma_o=\sqrt{1-l^2-m^2},\qquad
\mathbf d=l\mathbf e_1+m\mathbf e_2+\gamma_o\mathbf e_3.
\tag{W-11}
$$

l²+m²<1의 propagating 방향만 사용한다. Baseline은 chief 축을 demodulation한 경우다. 공통 axis에 대해 다른 plane-wave carrier `(l₀,m₀)`를 제거했다면 실제 방향은 `l=l₀+λ_mν_u`, `m=m₀+λ_mν_v`로 복구해야 한다. 서로 다른 색의 chief tilt를 개별 제거한 채 합산하면 lateral angular color가 사라진다.

Afocal `lens_metric_pupil`의 독립 기본은 위에서 고정한 canonical physical Π의 통과 pupil에서 `A_Π=1`이다. 더 일반적으로 `A_Π(u,v)=a_Π(u,v)`로 지정하고 power transmission을 포함하는 옵션에서만 `A_Π=a_Π√τ`로 둔다. Mask M은 실제 모든 aperture를 통과한 support에 한 번 적용한다. 이 규약에는 focal의 B나 거리 L이 필요하지 않다. Afocal `source_transport_pupil`은 (W-6)의 실제 source plane-flux 진폭이다. 평면을 임의로 옮긴 뒤 A를 다시 균일하게 놓으면 다른 조명·지표가 되므로 동일 계산으로 취급하지 않는다. Π 위치·illumination measure·고정/갱신 정책을 결과에 저장한다. FFT pitch는 `Δl=λ_m/(N_uΔu)`, `Δm=λ_m/(N_vΔv)`다. Direction window와 angular PSF tail을 확인하며, FFT가 표현하지 않은 방향의 power를 계산창 정규화로 감추지 않는다.

RS far field의 obliquity를 포함하면 unit solid angle intensity shape는 다음과 같다. [W01](#ref-W01) §2. 이 profile의 measure 변환은 자체 유도다.

$$
J_\Omega(l,m)\propto\gamma_o^2
\left|F\!\left(\frac l{\lambda_m},\frac m{\lambda_m}\right)\right|^2,
\qquad
d\Omega=\frac{dl\,dm}{\gamma_o},
$$
$$
h_\Omega=\frac{\gamma_o^2|F|^2}
{\int\gamma_o|F|^2 dl\,dm},\qquad \int h_\Omega d\Omega=1.
\tag{W-12}
$$

절대 W/sr가 필요할 때만 전기장 convention의 `n_i/(2Z₀λ_m²)`와 SI 단위를 복원한다. 정규화 h의 형상과 source throughput Q는 따로 둔다. High angle에서 γ²를 빼면 다른 모델이 된다. 이 단일 물점 응답에 Q를 곱하고 source-area 적분해야 유한 LED 전체의 far-field intensity가 된다.

보고 좌표를 `α=atan2(l,γ_o)`, `β=atan2(m,γ_o)`로 선택하면

$$
d\Omega=J_{\alpha\beta}\,d\alpha\,d\beta,\qquad
J_{\alpha\beta}=\frac{\sec^2\alpha\sec^2\beta}
{(1+\tan^2\alpha+\tan^2\beta)^{3/2}},
\qquad h_{\alpha\beta}=h_\Omega J_{\alpha\beta}.
\tag{W-13}
$$

Uniform l,m에서 uniform α,β로 옮길 때 cell power를 보존한다. 실제 angular OTF는 `h_αβ`를 α,β에 대해 Fourier 변환한다. 작은 PSF에서 l≈α,m≈β 및 γ≈1로 계산하는 단순 FFT angular MTF는 별도 `small_angle` 근사다. 그 조건 없이 FFT array의 l index를 rad라고 부르지 않는다. 주파수는 cycles/rad이며 `f_cycles/mrad=f_cycles/rad/1000`이다. Focal의 120 lp/mm를 자동으로 동일 숫자의 angular band로 바꾸지 않는다.

위 α,β는 단일 field의 **local chief frame** 좌표다. 여러 LED 물점의 hΩ를 합칠 때는 (W-11)의 3D 방향 d를 시스템의 공통 frame으로 회전하고 global angle bin에 power를 보존하여 다시 넣는다. 서로 다른 chief 원점의 angular 배열을 index가 같다는 이유로 더하지 않는다. Rotation 자체는 dΩ를 보존하지만 각도 chart의 cell 면적은 다시 계산해야 한다.

Afocal DL도 같은 amplitude·mask에서 plane OPD=0으로 만든다. Uniform circular physical pupil 지름 D에 대한 `f_c≈D/λ_m` cycles/rad는 작은 angular coordinate/γ≈1 검산식이다. Exact angular coordinate, obliquity 및 measure를 유지한 (W-12)–(W-13) 전체 결과에 동일 원형 closed form을 무조건 강제하지 않는다. 공간 주파수와 각도 주파수를 연결할 필요가 있으면 실제 국소 map `x=Jq`의 `f_x=J^{-T}f_q`를 쓴다. Global EFL 하나로 모든 field의 mapping을 대체하지 않는다.

### W.9 RGB 합성과 위치 등록

주요 결과는 G 단색이다. RGB lens-only 합성 또는 source image가 필요한 경우, 모든 λ의 PSF를 **같은 물리 출력면·같은 좌표·같은 pixel area**에 놓고 intensity로 합한다. Focal에서 G의 q₀와 e basis를 공통으로 쓰고 `W_λ=S_λ+n_i(λ)R−C_λ`로 평가하면 lateral color와 chromatic defocus가 남는다. `C_λ`는 임의 piston이다. λ마다 chief/centroid를 원점으로 맞춘 뒤 단순 평균하지 않는다.

$$
H(\mathbf x)=\sum_\lambda w_\lambda Q_\lambda h_\lambda(\mathbf x),\qquad
h_{\rm poly}=\frac{H}{\int H dA},\qquad
\mathrm{MTF}_{\rm poly}=\left|\mathcal F\{h_{\rm poly}\}\right|.
$$

Lens-only 모드에서 Q를 포함하지 않기로 했다면 유효 spectral weights에 그 정책을 명시하고 Q=1로 둔다. Source 모드에서는 실제 Q와 source/detector spectral weights를 사용한다. 파장별 MTF magnitude를 직접 평균하지 않는다. 서로 incoherent인 RGB의 **복소 field를 더하지 않는다**. 필드별 PSF를 다른 격자로 계산했다면 flux-conserving interpolation과 registration 오차를 확인한다.

### W.10 연결 검증과 아직 열린 문제

| 시험/쟁점 | 구현 규칙과 기대 관계 |
|---|---|
| Boundary eikonal gradient | Smooth interior에서 `∂u S=n_i d·e₁`, `∂v S=n_i d·e₂`. 유한차분 step을 줄여 수렴 확인. 경계·caustic 제외. |
| OPD gradient | (W-1)에서 `∂u W=n_i(d·e₁+u/R)`, `∂v W=n_i(d·e₂+v/R)`. 이상 수렴 ray이면 0. 이는 raw phase/reference/sign을 함께 검사한다. |
| Piston/tilt | Piston은 PSF 불변. `W→W+aχ+bη`는 (W-4)의 PSF를 `(a/n_i,b/n_i)`만큼 이동시켜야 한다. MTF는 유지되지만 complex OTF phase는 바뀐다. |
| Defocus | 물리 IMA/허용 compensator 변화 후 W가 변한다. Best-fit sphere로 defocus를 삭제한 결과는 fixed IMA 결과로 저장하지 않는다. |
| 독립 scalar 적분 | 일부 field·PSF 점에서 (W-2)를 직접 quadrature하여 국소 FFT 근사·보간을 검사. 같은 ray-boundary 입력의 오류까지 독립적으로 검증한 것은 아니다. |
| Geometry+λ similarity | 모든 길이·STOP·OBJ·IMA·λ₀를 s배, n과 대응 조명을 고정하면 `h_s(x)=s⁻²h₀(x/s)`, `OTF_s(f)=OTF₀(sf)`. 정규화 angular 응답은 불변. 실제 온도 변화 법칙이 아닌 합성 시험. |
| Raw source-power closure | 아래 (W-14)를 PSF 정규화 전에 계산하여 physical Π flux 및 ray source-cell power와 비교한다. Q×정규화 h의 적분은 정의상 Q이므로 그 검사만으로 진폭 전달의 정확성을 입증하지 않는다. |
| 계산 수렴과 모델 유효성 | Ray/phase grid 수렴, padding/창 수렴, scalar/near-focus 유효성, vendor 일치를 별도 상태로 기록. 한 종류의 통과로 다른 종류를 대체하지 않는다. |
| 실제 IMA와 FFT 가상면 | 기본 결과는 chief-normal. Tilted physical detector의 정확한 전파는 별도 계산. 국소 J 투영은 조건부 근사. |
| 제품 cosine-space/진폭/축 | 독립 profile은 닫혔지만 제품 mapping·working F/#별 T/S scale·polychromatic registration·sampling indexing 전체의 수치 동등성은 미확인. Reference output 확보 후 차례로 대조한다. |
| Pupil distortion/branch | 정확한 ray map을 사용해도 regular scalar boundary 가정은 남는다. Fold/caustic을 큰 amplitude로 덮지 않는다. |
| 고 NA/복소 coating | Scalar RS도 vector polarization·종방향 field·coating phase를 자동 포함하지 않는다. 필요하면 확장하며 FFT metric과 이름을 구분한다. |

**정규화 전 source transport의 전력 진단.** (W-11)–(W-12)의 `F`가 실제 전기장 `A_{Π,E}`를 사용한 physical-plane Fourier 변환이면, 모든 길이를 m로 통일하여 propagating hemisphere에서

$$
\Phi_{\rm wave}=\frac{n_i}{2Z_0\lambda_m^2}
\int_{l^2+m^2<1}\gamma_o
\left|F\!\left(\frac l{\lambda_m},\frac m{\lambda_m}\right)\right|^2\,dl\,dm,
\qquad
\Phi_\Pi=\int_\Pi p_\Pi\,dA_\Pi,
\qquad
\Phi_{\rm ray}=Q(\mathbf q_s)\Delta A_s.
\tag{W-14}
$$

여기서 λₘ 및 F의 입력 u,v는 m, pΠ는 W/m², Q는 W/m²(source), ΔAₛ는 m²다. 모두 같은 source cell·파장·직접 통과 경로·τ 정책을 사용해야 한다. `Φwave/ΦΠ`, `ΦΠ/Φray`, 비정규화 차이, 적분 angular domain, 그 밖의 미평가 power 추정을 저장한다. 0 분모는 ratio를 null로 둔다. Angular window를 제한했다면 Φwave는 그 창의 부분 전력임을 표시한다.

이 검사는 scalar angular spectrum의 전력과 GO ray-density 경계장이 일관되는지 확인한다. GO 기반의 국소 plane-wave 전력 식, 유한 창, sampling, edge 근사 및 scalar 모델의 차이 때문에 모든 처방에서 정확한 등식을 강제하는 것은 아니다. 진폭을 임의 상수로 재조정하여 mismatch를 지우거나 Q×h의 정규화 후 적분으로 대체하지 않는다. 큰 차이는 단위·cell measure·진폭·투과율 중복·window/aliasing·모델 한계 중 원인을 조사할 신호다. 이것이 full vector Maxwell 에너지 검증을 대신하지는 않는다.

Jones 추적으로 얻은 **ray별 총 power만** scalar 진폭으로 바꾸는 경우에는 pupil 위치별 polarization 방향 차이가 만드는 coherent interference를 포함하지 못한다. 일반적으로 `|Σ_j E_j exp(iφ_j)|²`는 `|Σ_j |E_j| exp(iφ_j)|²`와 다르다. 비간섭인 두 입력 편광 상태는 각각 pupil의 vector field를 coherent하게 전파한 후 intensity를 더해야 한다. 현재 scalar profile에서 Jones power를 transmission weight로만 사용하는 것은 명시적 근사이며, 이를 완전한 vector PSF로 명명하지 않는다.

최소 result provenance는 `profile_id/version`, λ₀/nᵢ, field/object mode, Π/q₀/basis/L, OPD reference와 제거 항, mask/B/illumination 정책, actual support/실제 traced count/padded shape/spacing, x·f coordinate arrays, raw normalization, truncation/convergence, (W-5), 실패 상태다. 어떤 숫자가 아직 실제 engine에서 검증되지 않았는지도 저장한다.

---

<a id="part-5"></a>

## 5. 전력·RI·원거리 분포·영상·grid

이 장은 공통 실광선 추적기에서 얻는 교점·방향·경로 길이·매질·aperture 통과 상태를 사용한다. 기하광학적 직접 경로만 추적하며, 표면 반사로 빠지는 전력은 손실 장부에 기록한다. 반사 ghost·산란·형광·mount 영향은 이 모델에 포함하지 않는다. 회절 영상은 별도의 파면 전파 장에서 정의한 **조명 및 적분 측도가 일치하는** PSF를 받는다. 광선 히스토그램을 회절 PSF로 부르지 않는다.

아래 식의 `[유도]`는 광량 정의·좌표변환으로부터 이 문서에서 정리한 계산식이다. 제품 문서가 같은 코드나 동일한 수치 경로를 사용한다는 뜻이 아니다. 출처 번호 S는 기존 출처 목록과 동일하다.

### R.1 계산량·단위·광원 계약

세 종류의 출력은 구분하여 저장한다.

| 출력 | 정의·단위 | 분모 또는 기준 |
|---|---|---|
| 영상면 irradiance E | 면적당 입사 전력, W/mm² | 실제 IMA 또는 명시한 검출면의 실제 면적 |
| Sequential형 RI | 균일 장면에 대한 E의 상대값 | 선택한 field 영역에서의 peak; 축상 정규화는 다른 출력 |
| Far-field radiant intensity I | 입체각당 출사 전력, W/sr | 전체 유한 광원을 합한 출력 각분포 |
| 원거리 상대 각분포 | I의 상대값 | 선언한 각도 영역의 peak 또는 축상 값 |
| Φ | 전력, W | source·수집·출사·검출 ROI 중 기준을 명시 |

기하 단위는 mm, solid angle은 sr, 내부 각도는 rad이다. 이 장에서 L은 대표 선스펙트럼에 배정된 radiance [W/(mm²·sr)]이며 nm당 spectral density가 아니다. 연속 spectrum 자료를 받는 경우 대표 파장의 구간 적분 또는 가중치로 먼저 변환한다. W/(m²·sr) 입력을 내부 단위로 바꿀 때 10⁻⁶을 곱한다. 결과 E[W/mm²]를 W/m²로 출력할 때 10⁶을 곱한다. μ[mm⁻¹]와 실제 광선 길이 ℓ[mm]를 사용한다. 상대 입력 L=1도 계산할 수 있지만 이때 W 단위의 실제 밝기라고 표시하지 않는다.

유한 LED 모델을 다음으로 고정한다.

- source 평면의 로컬 좌표 q=(qₓ,qᵧ), 영역 Aₛ=[−6.4,6.4]×[−1.6,1.6] mm², 중심은 광축 위이다.
- source normal을 nₛ라 하고, 전방 반구에서 Lₛ(q,ω,λ)=L₀,λ로 일정한 Lambertian radiance를 사용한다. 각방향 방출 전력에는 cosθₛ가 붙지만 radiance 자체에 cosθₛ를 다시 곱하지 않는다.
- 320×80, pitch 0.040 mm, fill factor f=0.64는 균질화 광원의 메타데이터이다. 발광부 radiance가 입력이면 L₀,λ=fLactive,λ. 이미 패널 평균 radiance 또는 총광량이 입력이면 fill을 다시 적용하지 않는다.
- 전체 패널 면적은 40.96 mm², active area 합은 26.2144 mm²이다. 정사각 발광부로 해석할 때 변 길이는 0.032 mm지만, 현재 모델은 개별 pixel 경계를 추적하지 않는다.
- 공간적으로 서로 다른 source cell 및 서로 다른 대표 파장은 비간섭이다. 한 source point·한 대표 파장의 pupil 성분은 회절 PSF 계산에서 coherent하게 합한다. 각 점에서 출발한 spherical-wave eikonal과 Lambertian angular amplitude를 사용하는 scalar point-emitter 모델을 선언한다. Lambertian이라는 말만으로 이 coherence 가정이 자동으로 정해지는 것은 아니므로 별도 모델 항목으로 기록한다.

반구 방출 전력은 다음과 같다. [유도; 기본 radiometric integral S21]

$$
\Phi_{s,\lambda}=\int_{A_s}\int_{\Omega_+}L_{s,\lambda}(\mathbf q,\boldsymbol\omega)\cos\theta_s\,d\Omega_s\,dA_s
=\pi A_sL_{0,\lambda}.
$$

절대 광량이 미정이면 RI나 정규화 원거리 각분포는 임의의 공통 양의 L₀로 계산할 수 있다. RGB별 상대 스펙트럼이 미정이면 RGB 합성 밝기는 정해지지 않는다. G의 렌즈 단독 FFT MTF에는 LED pitch·fill·bitmap·반구광량을 넣지 않는다.

### R.2 표면·bulk 투과와 polarization 경계

#### R.2.1 무코팅 유전체의 전력 계수

각 표면의 실광선 입사각 i, 굴절각 t와 **해당 λ,T의 실제 매질 굴절률** n₁,n₂를 사용한다. 양의 cos i, cos t를 갖는 전파 해에 대해 다음 전력 계수를 계산한다. [S06](#ref-S06)

$$
R_s=\left|\frac{n_1\cos i-n_2\cos t}{n_1\cos i+n_2\cos t}\right|^2,
\qquad
R_p=\left|\frac{n_2\cos i-n_1\cos t}{n_2\cos i+n_1\cos t}\right|^2.
$$

$$
T_s=1-R_s,\qquad T_p=1-R_p.
$$

이 식은 경계 자체가 흡수하지 않는 매끈한 무코팅 유전체에 해당한다. TIR이면 직접 투과 계수는 0, 해당 표면의 반사 손실은 입사 전력 전부다. TIR을 수치 실패와 혼합하지 않는다.

- `lossless_scalar`: 표면·bulk 감쇠를 끄고 aperture 손실만 계산한다. 이상적 비교용 모델이다.
- `scalar_interface_average`: 표면마다 T=(Tₛ+Tₚ)/2를 곱한다. 최초 무편광 광의 **한 경계**에서는 정확하지만 여러 경계가 만드는 편광 상태를 유지하지 않는 근사다. 매 경계 전에 다시 무편광으로 만드는 가정과 동등하며 일반 유리 렌즈가 실제로 그렇게 depolarize한다는 뜻은 아니다.
- `unpolarized_jones_power`: 두 직교 초기 편광을 독립적으로 수송하고, 마지막 전력을 평균한다. 여러 경계에서 생성되는 편광 차이와 S/P basis 회전을 보존한다. 이 결과를 단순 scalar FFT의 전력 가중치로 쓸지, vector 성분별 파동 전파로 쓸지는 별도 설정이다.

두 번째 모드를 세 번째와 동일하다고 표기하지 않는다. Zemax의 `Use Polarization`도 분석 기능마다 intensity 가중과 vector 파동 계산의 범위가 다르다. [S51](#ref-S51)

실제 구현의 한 방법은 광선에 수직인 직교 basis의 **power-normalized Jones vector**를 사용하여, 표면마다 local S/P basis로 회전하고 diag(gₛ,gₚ)를 적용하는 것이다. 무코팅 경계에서 전기장 투과 계수와 power-normalized 계수는 다음이다. [Fresnel 관계의 유도]

$$
t_s=\frac{2n_1\cos i}{n_1\cos i+n_2\cos t},\quad
 t_p=\frac{2n_1\cos i}{n_2\cos i+n_1\cos t},\quad
 g_{s,p}=\sqrt{\frac{n_2\cos t}{n_1\cos i}}\,t_{s,p}.
$$

기하 basis 회전·표면 전파를 포함한 전체 2×2 행렬을 J라 하면 무편광의 직접 투과율은 τpol=½Tr(J†J)이다. 각 surface에서 임의의 basis 부호가 바뀌어 phase 불연속이 생기지 않도록 연속 basis를 사용한다. 정상 입사처럼 S/P plane이 유일하지 않은 경우 이전 basis를 연속적으로 수송한다. `T`는 전력비이므로 scalar amplitude 감쇠에는 √T를 사용한다. 이미 power-normalized Jones 행렬에 T가 들어갔다면 scalar T를 또 곱하지 않는다.

코팅은 복소 진폭 투과계수 또는 각도·파장별 전력투과 데이터와 정확한 convention이 있어야 한다. 단순 n·Abbe 값만으로 AR coating, coating phase, 흡수 및 산란을 복원하지 않는다. intensity 표만 있다면 전력 전달에는 사용할 수 있으나 coating phase가 필요한 coherent wave 계산은 지원 범위가 제한된다.

#### R.2.2 벌크 흡수

균질한 재료의 직접 경로 길이 ℓ에 대해

$$
\tau_{bulk}(\lambda,T,\ell)=e^{-\mu(\lambda,T)\ell},\qquad
\mu=-\frac{\ln\tau_{int}(\ell_0)}{\ell_0}.
$$

τint는 표면 Fresnel을 제외한 **내부 투과율**이다. 측정 두께 ℓ₀와 transmission 종류를 재료 입력에 저장한다. 총 transmission 자료를 내부 투과율로 잘못 사용하고 표면 손실을 한 번 더 적용하지 않는다. “transmission range”만 있고 두께당 transmission 또는 μ가 없으면 절대 bulk absorption은 결정되지 않는다. 실제 ray의 ℓ를 사용하며 렌즈 center thickness로 일괄 치환하지 않는다.

#### R.2.3 직접 전력 장부

각 입사 packet Φin의 상태를 propagation 순서로 갱신한다.

1. 균질 구간 후 Φafter=τbulk Φbefore, 차이는 bulk absorbed에 더한다.
2. 유효 aperture 밖이면 잔여 packet 전체를 blocked에 기록하고 종료한다.
3. 표면에서 Φtrans=TΦincident, Φrefl=(1−T)Φincident를 기록하고 직접 투과 가지를 계속한다. 복소 coating에 흡수가 있으면 T+R+A=1의 별도 계수를 쓴다.
4. 교점 미수렴·상태 판별 실패이면 마지막으로 확정된 잔여 전력을 unresolved에 보존한다.
5. 마지막 경로에서 실제 ROI 밖으로 나간 광량은 out-of-ROI이며, optical absorption이 아니다.

$$
\Phi_{launched}=\Phi_{direct}+\Phi_{reflected}+\Phi_{absorbed}+\Phi_{blocked}+\Phi_{missed}+\Phi_{unresolved}
$$

이 등식의 항들이 서로 배타적이어야 한다. 앞에서 Fresnel 반사로 제거한 광량을 뒤 surface의 차단 전력에 다시 포함하지 않는다. 방향 sampling을 반구 일부로 제한했다면 Φlaunched는 그 부분의 광량이다. 전체 반구 방출 Φs 대비 효율을 보고할 때 누락된 방향을 `not_collected/not_launched`로 별도 처리한다. Solver 실패를 손실 0 또는 정상 차단으로 바꾸면 전력 보존 시험의 의미가 사라진다.

### R.3 광선의 초기 전력: sampling과 물리량의 연결

#### R.3.1 source 면적·출사 입체각을 직접 적분

한 source point와 한 출사 방향이 결정하면 직접 Sequential 경로는 공통 trace로 결정된다. (q,ω)의 결합 확률밀도 p(q,ω)는 **dAₛ dΩₛ 기준**이며 ∫p dAₛdΩₛ=1이다. N개 표본을 쓰면 packet 초기 전력은 [유도]

$$
w_k=\frac{L_s(\mathbf q_k,\boldsymbol\omega_k)\cos\theta_{s,k}}{N\,p(\mathbf q_k,\boldsymbol\omega_k)}.
$$

독립 추출이 아니면 pA(q)pΩ(ω)로 임의 분해하지 않는다. 조건부 추출이면 p=pA(q)pΩ(ω|q)를 사용한다. deterministic quadrature에서는

$$
w_{ab}=L_s(\mathbf q_a,\boldsymbol\omega_b)\cos\theta_{s,b}\,w_{A,a}w_{\Omega,b}
$$

이며 wA,wΩ가 각각 면적·solid-angle 적분 가중치다. 좌표가 θ,φ이면 dΩ=sinθ dθdφ를 포함한다. 균일 θ,φ 표본을 균일 solid angle로 간주하지 않는다.

유용한 fixture는 균일 panel와 cosine-weighted hemisphere sampling이다. pA=1/Aₛ, pΩ=cosθ/π이면 모든 초기 packet의 w=πAₛL₀/N가 같다. 이때 Φs를 모든 packet이 정확히 나눠 갖는다. 광학계가 작아 대부분이 빗나가면 아래와 같은 entrance/reference-plane importance sampling을 사용할 수 있다. 확률밀도와 수집 효율 기준은 그대로 남겨야 한다.

#### R.3.2 source에서 보이는 앞쪽 평면을 sampling

source와 첫 굴절면 사이에 놓인 기준 평면 B의 점 b를 택한다. 두 점 사이에 굴절면이 없고, 거리 r=|b−q|, 그 평면 normal과 광선의 각도를 θB라 하면 [S21; 좌표변환]

$$
d\Omega_s=\frac{|\cos\theta_B|}{r^2}\,dA_B,
\quad
 d\Phi=L_s\frac{\cos\theta_s|\cos\theta_B|}{r^2}\,dA_s\,dA_B.
$$

두 cos는 서로 다른 투영에 해당한다. source Lambertian radiance에 cos를 미리 곱한 뒤 위 식의 cosθs까지 다시 쓰는 것은 중복이다. 기준 평면 좌표 대신 STOP 좌표를 역으로 aiming하여 발사한다면, STOP 이전 굴절 때문에 이 자유공간 r⁻² 식을 그대로 적용할 수 없다. 그때는 aiming map의 Jacobian

$$
d\Omega_s=\left|\frac{\partial\boldsymbol\omega_s}{\partial a}\times\frac{\partial\boldsymbol\omega_s}{\partial b}\right|\,da\,db
$$

를 계산하거나, 출사 방향 PDF로 돌아간다. 미분은 unit-sphere 방향벡터의 미분이다. 국소 folded map, 다중 aiming root, 누락 branch가 있으면 단순한 하나의 Jacobian으로 전력을 보장할 수 없다.

### R.4 Radiance 보존과 영상면 RI

#### R.4.1 보존 법칙

손실 없는 isotropic 경로에서 L/n²는 보존되고, 직접 경로의 전력투과율을 τdirect로 모으면 [S59; étendue 관계의 재정리]

$$
L_i=\tau_{direct}\left(\frac{n_i}{n_s}\right)^2L_s,
\qquad dG=n^2\cos\theta\,dA\,d\Omega.
$$

이는 광량이 굴절 후 증폭된다는 뜻이 아니다. 공간·각도 측도의 변화이다. 이미 W 단위 전력을 갖는 ray packet에는 n²를 또 곱하지 않는다. air-to-air 시스템에서 시작·끝 n의 비는 상쇄되지만 내부 ray 굴절과 surface T는 여전히 계산한다.

#### R.4.2 균일 장면에서의 실제 irradiance 적분

영상면 점 x의 로컬 normal과 입사 방향 사이의 양의 cosine을 cᵢ라 하면

$$
E_\lambda(\mathbf x)=\int_{\Omega_{visible}}L_{i,\lambda}(\mathbf x,\boldsymbol\omega_i)c_i\,d\Omega_i.
$$

영상면의 로컬 직교 접평면 방향 cosines를 l,m이라 두면 cᵢ=√(1−l²−m²), cᵢdΩᵢ=dl dm이다. 따라서 [S21; Jacobian 유도]

$$
E_\lambda(\mathbf x)=\left(\frac{n_i}{n_s}\right)^2\int_{D(\mathbf x)}\tau_{direct}(\mathbf x,l,m,\lambda)L_{s,\lambda}(\mathbf q_s,\boldsymbol\omega_s)\,dl\,dm.
$$

구현 절차는 다음과 같다.

1. 요청 field q의 평가 위치 x=gλ(q)를 선언한다. 일반적으로 해당 field의 chief-ray IMA 교점이다. Chief가 차단되면 임의의 0을 넣지 말고, 기하 mapping 정의의 대체 기준 또는 평가 불가를 기록한다.
2. x에서 광학계 방향으로 역광선을 발사하여 visibility와 앞쪽 source/domain 연결을 찾는다. 광학적 역방향 법칙으로 경로를 찾되 물리 광량은 forward source→image convention으로 평가한다.
3. image-direction-cosine disk에 일정 간격 격자 또는 적응형 적분을 놓는다. 각 cell의 dl dm 면적을 사용하고 추가 cosθi를 곱하지 않는다.
4. 실제 clear aperture/STOP/edge에서 차단된 방향의 integrand는 0이다. 유효 경로의 τdirect와 n비·source radiance를 넣는다.
5. 이미지 평면에서 여러 x를 적분하여 E map을 얻는다. image surface가 곡면이면 각 점의 실제 normal과 실제 면적을 사용한다. XY 투영 이미지의 mm²와 곡면 mm²는 다르다.

이것은 **정의에 기반한 독립적인 reverse radiance 적분기**다. 임의 pupil 표본을 uniform cosine 표본으로 재분류하는 방법이 아니다. Forward pupil 좌표 (u,v)로 적분하는 경우에는 |det ∂(l,m)/∂(u,v)|가 필요하다. 원점이 x가 아닌 기존 bundle에서 방향만 모은 것은 일반적으로 위 적분과 같은 경로 집합이 아니다.

Sequential형 RI는 source가 충분히 넓은 균일 평면 또는 정의한 무한 장면이라고 가정한다. 반면 source를 정확히 12.8×3.2 mm로 자르면 가장자리 field의 irradiance에는 유한 source 경계가 들어간다. 두 장면을 구분하여 `uniform_scene_ri`와 `finite_source_irradiance`를 저장한다. 동일 source 경계로 계산하지 않은 forward/reverse 결과를 오차로 비교하지 않는다.

#### R.4.3 RI 정규화 및 Zemax 대응

평가할 field/domain H를 기록하고

$$
E_{peak,\lambda}=\max_{\mathbf q\in H}E_\lambda(g_\lambda(\mathbf q)),\qquad
RI_{peak}(\mathbf q,\lambda)=\frac{E_\lambda(g_\lambda(\mathbf q))}{E_{peak,\lambda}}
$$

로 계산한다. 실제 이산 계산에서는 `sampled_peak`와 추정 peak 위치·field refinement 정보를 저장한다. 전역 최대를 증명하지 않은 것을 `exact_global_peak`라고 표시하지 않는다.

Zemax 공식 RI는 **축상 값이 아니라 field 안의 최대 illumination**으로 정규화하며 단색 분석이다. 문서는 uniform Lambertian scene, projected exit-pupil/cosine-space 적분과 Rimmer 방법의 확장을 설명한다. 충분히 국소화된 공액 영상과 angle-space caustic이 없는 상황 등 적용 조건이 있다. [S20](#ref-S20) 독립 reverse 적분기의 물리적 목적은 같지만, field 위치·sampling·apodization·product 옵션까지 대조하기 전에는 동일 수치를 보장하지 않는다.

별도로 RIaxis=E/Eaxis를 제공할 수 있으나 이름과 분모를 구분한다. Epeak>0인 상태에서 어떤 field가 완전 차단되면 RI=0은 정상이다. 전체 Epeak=0이면 상대값은 undefined이다. Solver 미수렴으로 얻은 낮은 Epeak는 정규화 기준으로 확정하지 않는다. 온도/설계 간 비교를 위해 raw E 또는 공통 nominal 기준의 상대값도 보존한다.

cos⁴θ는 thin·slow·무수차 렌즈, 렌즈 위치 STOP 등의 제한된 fixture에서만 점검한다. 일반 비구면·telecentric·비네팅 시스템에 전체 조도의 계산식으로 넣지 않는다. [S20](#ref-S20) 이 장의 projected-solid-angle integral에 cos⁴를 추가로 곱하면 중복이다.

### R.5 유한 Lambertian 광원의 far-field

#### R.5.1 직접 ray-transport 각분포

각 packet의 최종 방향 dₖ와 직접 출사 전력 wₖτₖ를 얻는다. 각도 bin B의 정확한 solid angle ΔΩB를 사용하여

$$
\Phi_B=\sum_{k:\mathbf d_k\in B}w_k\tau_k,\qquad
I_B=\frac{\Phi_B}{\Delta\Omega_B}.
$$

I는 bin 평균 radiant intensity [W/sr]이다. 출사 방향 histogram에는 detector incidence cosine을 다시 곱하지 않는다. `per_field_angular_spot`은 하나의 source point에서 출사하는 방향 분포이고, 위 I는 **전체 source area와 emission angle을 합한** 값이다. 둘을 같은 출력으로 저장하지 않는다. 원거리 I의 계산은 방향 분포를 다루므로 가상의 IMA를 매우 큰 거리로 이동시켜 좌표를 빼는 방식이 필요하지 않다.

구면 좌표 θ,φ이면 dΩ=sinθdθdφ이고 직사각 θ–φ bin의 solid angle은

$$
\Delta\Omega=(\phi_2-\phi_1)(\cos\theta_1-\cos\theta_2).
$$

전방 direction-cosine chart d=(l,m,c), c=√(1−l²−m²)>0이면 dΩ=dl dm/c이다. 단순한 ΔlΔm를 ΔΩ로 쓰지 않는다. 투사 화면식 좌표 p=dₓ/dz, q=dᵧ/dz에서는

$$
d\Omega=\frac{dp\,dq}{(1+p^2+q^2)^{3/2}}.
$$

수평·수직각을 α=atan p, β=atan q로 정의한 경우에만

$$
d\Omega=\frac{\sec^2\alpha\sec^2\beta}{(1+\tan^2\alpha+\tan^2\beta)^{3/2}}\,d\alpha\,d\beta.
$$

β를 통상 구면 elevation으로 정의한 convention과 혼합하지 않는다. 넓은 bin은 중앙 Jacobian 곱 대신 bin의 solid angle을 적분한다. 위 관계는 unit-sphere 면적에서 직접 유도한 좌표변환이며 Zemax의 각도 좌표를 암묵적으로 복제한 것이 아니다.

출사 효율과 상대 분포는 다음처럼 분리한다.

$$
\eta_{direct}=\frac{\Phi_{exit}}{\Phi_s},\qquad
RI^{far}_{peak}(\boldsymbol\omega)=\frac{I(\boldsymbol\omega)}{\max_{\Omega_{report}}I},\qquad
RI^{far}_{axis}(\boldsymbol\omega)=\frac{I(\boldsymbol\omega)}{I(\hat{\mathbf z})}.
$$

출력 이름은 `relative_far_field_intensity`로 한다. 이는 focal Sequential RI 분석과 다른 물리량이다. Peak bin은 angular resolution에 따라 변할 수 있으므로 raw ΦB, ΔΩB, sampling, smoothing 유무를 보존한다. Smoothing 결과만 데이터로 남기지 않는다. [Radiant intensity의 제품 정의 S47]

광학계 없는 panel fixture는 I(θ)=AₛL₀cosθ, ∫I dΩ=πAₛL₀를 만족한다. Ray 수·source quadrature·angle bin refinement를 각각 확인한다. 설계 간에는 같은 source radiance·정의한 각도 영역·bin을 사용해야 한다.

#### R.5.2 회절을 포함하는 각분포

이전 절은 기하광학적 각분포다. 회절을 포함하려면 각 source point q의 물리 조명에 맞는 angular PSF hΩ,λ(ω|q)를 구한다. hΩ는 global sphere에서 ∫hΩdΩ=1이다. Source point의 수집·직접 투과 전력밀도는

$$
Q_\lambda(\mathbf q)=\int_{\Omega_{s,+}}L_{s,\lambda}(\mathbf q,\boldsymbol\omega_s)\cos\theta_s\,\tau_{direct}(\mathbf q,\boldsymbol\omega_s)\,d\Omega_s.
$$

Q의 단위는 W/mm² source이다. 통과 불가능한 방향은 τdirect=0이다. 따라서

$$
I_\lambda(\boldsymbol\omega)=\int_{A_s}Q_\lambda(\mathbf q)\,h_{\Omega,\lambda}(\boldsymbol\omega\mid\mathbf q)\,dA_s.
$$

이는 비간섭 source cell의 intensity 합이다. 파면 복소 amplitude를 서로 다른 q 사이에서 먼저 합하지 않는다. Local angular PSF를 global sphere에 배치할 때 chart 회전과 solid-angle Jacobian을 반영한다. 각도가 커져도 모든 local PSF를 같은 2D pixel 배열에서 평행 이동만 하는 것은 일반적으로 정확하지 않다.

각 chart 좌표 ξ에 대해 dΩ=JΩ(ξ)d²ξ라면 histogram/PSF density는 hξ=hΩJΩ이다. ∫hξd²ξ=1인 FFT 출력과 sr⁻¹인 hΩ를 혼합하지 않는다. 특히 direction-cosine chart에서는 JΩ=1/c이다. Afocal FFT 규약의 amplitude·flux weighting은 파동 전파 장의 정의를 따른다.

Ray histogram과 wave PSF 합성을 함께 사용하여 같은 회절 blur를 두 번 넣지 않는다. 무회절 결과와 회절 결과를 각각 명시하여 비교할 수 있다. 균질화 source가 허용하는 모델은 pixel-gap 고주파 효과를 생략한 연속 emission 모델이며, 이 한계는 RI 계산 가능 여부와 별개다.

#### R.5.3 필요할 때만 유한 스크린으로 변환

출사 pupil 크기를 무시할 만큼 원거리이며 원점에서 방출하는 angular distribution으로 근사할 수 있을 때, z=Z 평면에서 r=Z/cosθ이고

$$
d\Omega=\frac{\cos^3\theta}{Z^2}dA,\qquad E=I\frac{\cos^3\theta}{Z^2}.
$$

Z를 mm로 쓰면 E는 W/mm²이다. 이것은 RI의 cos⁴ 법칙이 아니라 angular distribution→평면 면적의 Jacobian이다. 가까운 screen에서는 실제 출사 교점과 방향으로 각 ray의 screen hit를 구한다. 이 보조 출력을 자동으로 기본 MTF 측정면에 적용하지 않는다.

### R.6 영상 합성: Q와 정규화 PSF의 책임 분리

#### R.6.1 물리적 커널 계약

Finite-source 이미지에서 source q의 실제 chief mapping을 gλ(q), source-conditioned 정규화 intensity kernel을 hA,λ(x|q)라 한다. hA는 전체 실제 검출면에서 ∫hAdAi=1이고 단위는 mm⁻² image이다. Kernel은 centroid 정렬 PSF+g인지, chief-ray 정렬 PSF+g인지 명시한다. 둘을 혼합하면 coma와 lateral color의 위치가 바뀐다.

$$
E_\lambda(\mathbf x)=\int_{A_s}Q_\lambda(\mathbf q)\,h_{A,\lambda}(\mathbf x\mid\mathbf q)\,dA_s.
$$

자주 쓰는 표기는 hA(x−gλ(q);q)이지만 이는 각 field에서의 중심 규약을 고정했을 때만 같은 의미다. 연속 source 적분을 이산화하면

$$
E_\lambda(\mathbf x_j)\simeq\sum_a Q_\lambda(\mathbf q_a)\,h_{A,\lambda}(\mathbf x_j\mid\mathbf q_a)\,w_{A,a}.
$$

실제 detector pixel 평균 또는 cell power는 PSF의 pixel 면적 적분을 사용한다. Pixel 중심값×면적은 충분히 작을 때의 quadrature다. Sample image의 pixel 크기, LED의 물리적 pixel pitch, 수치 detector cell 크기를 구분한다.

Q는 angular collection, surface/bulk 손실, aperture clipping을 포함한다. h에는 그 손실이 만드는 **pupil 내 상대 amplitude 형태**가 반영되어야 한다. Q를 곱했다고 h의 조명 형태가 자동으로 맞아지지는 않는다.

#### R.6.2 `lens_metric_pupil`과 `source_transport_pupil`

렌즈 FFT MTF의 uniform apodization pupil과 실제 Lambertian 광원의 pupil을 구별한다. 동일 ray geometry/OPL은 재사용할 수 있지만 다음이 다르면 amplitude, Q, PSF는 다시 계산한다.

- Source 방출 방향 분포와 수집각.
- Reference plane에 도달하는 ray tube 면적 및 Jacobian.
- Fresnel·bulk·coating의 field/pupil별 감쇠.
- 평가하는 출력면과 intensity의 면적/각도 측도.

물리 광원에서 작은 source cell ΔAs가 downstream reference plane Π로 보내는 flux density는 regular single-branch map에서 [유도]

$$
p_\Pi(\mathbf u)\,dA_\Pi=L_s(\mathbf q,\boldsymbol\omega_s)\cos\theta_s\,\Delta A_s\,\tau_{direct}\,d\Omega_s.
$$

즉 pΠ=LₛcosθₛΔAₛτdirect|dΩₛ/dAΠ|이다. 이 mapping을 ray tube 또는 수렴한 미분으로 계산한다. Caustic에서 Jacobian이 특이해지면 단순 보간 amplitude로 넘기지 않는다.

해당 plane에서 전기장 amplitude U를 쓰는 비자성 매질의 local plane-wave convention이라면

$$
p_\Pi^{SI}=\frac{n}{2Z_0}|U|^2(\mathbf d\cdot\hat{\mathbf n}_\Pi),\qquad
|U|=\sqrt{\frac{2Z_0p_\Pi^{SI}}{n(\mathbf d\cdot\hat{\mathbf n}_\Pi)}}.
$$

여기서는 pΠSI가 W/m²이며 내부 pΠ[W/mm²]에 10⁶을 곱한다. 이를 다른 의미의 flux amplitude나 이미 obliquity를 포함한 Fourier pupil에 다시 곱하지 않는다. Fourier pupil로 변환하는 kernel·Jacobian·reference phase는 파동 전파 장의 명시된 규약이 담당한다. `sqrt(packet_power)`를 비균일 plane 격자에 그대로 놓고 FFT하는 방법은 면적 측도가 빠져 있다.

h의 정규화는 PSF 계산 오차를 지우는 장치가 아니다. 유한 window에 저장된 energy fraction, tail/truncation 및 pupil sampling 수렴을 먼저 점검한다. 작은 window에 남은 일부 power를 임의로 1로 만들어 Q 전체를 분배하면 wing 손실이 숨겨진다.

#### R.6.3 왜곡 warp와 Jacobian 중복 방지

위 직접 source-coordinate 적분에서는 별도의 |det Dg|를 곱하지 않는다. Zero-blur limit에서 h=δ²(x−g(q))를 대입하면 [유도]

$$
E(g(\mathbf q))=\frac{Q(\mathbf q)}{|\det Dg(\mathbf q)|}
$$

가 자동으로 나온다. 먼저 source를 image grid에 warp하는 구현은 이 면적 변환을 포함해 flux를 보존해야 한다. 이미 warp된 irradiance 또는 vendor식 RI map에 Jacobian을 한 번 더 적용하지 않는다. Mapping이 접히면 여러 source branch의 기여를 합한다. Singular caustic에서는 점 density 대신 유한 detector-cell power로 적분한다.

물체가 무한대이면 q는 mm² source가 아니라 field angle chart다. 그때 입력은 radiance scene와 angular measure이며 `Q dA_source` 식에 mm²를 억지로 배정하지 않는다. 동일한 일반형을

$$
d\Phi_{collected}(\mathbf q)=Q_\mu(\mathbf q)d\mu_s(\mathbf q),\quad
E(\mathbf x)=\int Q_\mu(\mathbf q)h_A(\mathbf x\mid\mathbf q)d\mu_s
$$

로 쓰고 μs의 solid-angle Jacobian과 entrance projected area의 수집 적분을 포함한다. 유한 12.8×3.2 mm source를 무한 장면의 field angle로 변환하려면 실제 광학 mapping이 필요하다. 수치가 같다는 이유로 mm와 degree를 교환하지 않는다.

#### R.6.4 Field-varying kernel과 source bitmap

1. Bitmap의 orientation·source 좌표 범위·linear intensity 또는 radiance 해석을 고정한다. 파일이 감마 부호화되어 있으면 명시한 color profile로 linearize한다. 알 수 없는 encoding을 임의로 물리 radiance라고 간주하지 않는다.
2. 각 source sample에서 Ls,λ(q)=Bλ(q)Lbase,λ를 지정한다. 균질화 모델에서는 B가 macro-scale emission modulation이며 개별 pixel gap을 복원하지 않는다.
3. 축상·수평/수직 끝·corner와 중간 field에 source-conditioned h, Q, g를 구한다. Aperture clipping onset 또는 빠르게 변하는 수차 부근은 field grid를 세분한다.
4. 인접 field의 h를 동일한 물리적 좌표 간격과 중심 규약에 맞춘 후 nonnegative, energy-aware 방식으로 보간한다. g, Q의 보간과 PSF 모양 보간을 구분한다. Coherent amplitude를 field 사이에서 임의 보간·합산하여 source coherence를 만들어내지 않는다.
5. 소스 면적 quadrature로 h 기여를 합한다. 전체 영상에 한 PSF를 convolution하는 것은 검증된 isoplanatic 구간에서의 근사다.
6. Detector ROI밖으로 퍼진 광량을 따로 저장한다. 영상 edge에서 잘린 convolution을 자동 peak normalize하여 숨기지 않는다.

중심정렬 회전대칭 렌즈는 field 반경 방향의 재사용이 가능하지만 2D PSF/OTF를 회전해야 한다. T/S MTF 두 개만으로 일반적인 비대칭 2D PSF를 복원할 수 없다. 직사각 source 전체의 radiometric 적분은 반경 하나의 함수로 줄어들지 않는다.

#### R.6.5 RGB 및 표시 영상

동일한 global image/angle grid에 λR,λG,λB별 Eλ 또는 Iλ를 합한다. Source spectral weights가 이미 Lbase,λ에 들어 있으면 같은 weight를 또 적용하지 않는다. 센서/시각 반응을 모델링하는 경우 채널 c의 결과는 명시한 선형 응답 D(c,λ)에 대해 Cc(x)=ΣλD(c,λ)Eλ(x)다. 물리 RGB radiometric 영상과 display RGB는 다른 데이터다. White balance, tone mapping, exposure normalization은 시각화 단계이며 raw 결과를 덮어쓰지 않는다.

파장별 chief 위치를 독립적으로 중앙에 맞춘 뒤 합하면 lateral color를 지운다. 공통 물리적 위치를 유지한다. 세 대표 파장의 합은 세 선스펙트럼 모델이지 임의 LED 연속 spectrum의 완전한 색 예측이 아니다. Sample bitmap RGB와 λ 순서 mapping은 명시적으로 받는다.

Zemax Image Simulation의 공개 설명도 field별 PSF grid·보간·왜곡/색수차 처리를 사용한다. [S23](#ref-S23) 그러나 제품의 `Diffraction` 선택은 Huygens PSF 경로이며 자체 fallback·채널 처리·표시 정규화가 있다. [S52](#ref-S52) 본 엔진은 적용성이 확인된 FFT 물리-source kernel로 독립 영상 합성을 제공할 수 있다. 이를 모든 제품 Image Simulation 옵션과 수치적으로 동일하다고 부르지 않는다. Huygens는 모든 영상의 필수 선행 조건은 아니지만 FFT 적용 밖에서 아무 대체 모델 없이 값을 확정할 수도 없다.

### R.7 Grid simulation과 각 grid의 초점

Grid 좌표는 distortion 기준이 명확한 source 좌표에서 만든다. 유한 object는 mm, 무한 object는 선언한 field-angle chart다. Real image height로 grid를 역정의하면 왜곡을 입력 단계에서 상쇄할 수 있으므로 원본 좌표를 유지한다.

각 grid node qk에 대하여 다음을 저장한다.

1. 실제 chief mapping gλ(qk), reference mapping gref(qk), vector distortion Δg=g−gref.
2. 공통 IMA 또는 공통 angular chart에 놓인 실제 spot samples와 중심/RMS, 필요한 경우 source-conditioned PSF.
3. 선언한 common-focus sweep의 해당 위치에서의 성능. Field마다 별도로 얻은 최적 초점은 `local_best_focus_diagnostic`로 저장하고 동시에 달성되는 평면 영상이라고 표시하지 않는다.
4. Q(qk), clipping fraction, throughput, trace/applicability 상태.

직선을 촘촘히 표본화하여 g로 옮기면 왜곡 grid를 얻는다. 각 node에 spot/PSF를 함께 표시하면 local blur를 보여줄 수 있다. 선 자체의 물리적 blurred image가 필요하면 line/grid bitmap을 R.6의 동일 합성기로 계산한다. 점 표시의 크기 확대율은 legend에 쓰며 실제 길이와 혼동하지 않는다. Afocal grid는 각도 mapping과 angular PSF를 사용하고 임의 screen 거리로 spatial focus라고 바꾸지 않는다.

온도별 gref를 매번 새로 맞춘 distortion만 저장하면 전체 배율/FOV drift가 가려진다. `intrinsic_distortion`과 함께 nominal temperature의 g₀에 대한 ΔgT=gT−g₀를 저장한다. 비네팅이 MTF/spot을 개선하는 경우라도 Q와 Φexit을 같이 남겨 aperture trade-off가 보이게 한다.

### R.8 실제 구현의 수렴·교차검증·미결 사항

필수적인 검증 후보는 다음이다. 실제 엔진에서 수행하여야 하며 문서의 식만으로 통과했다고 볼 수 없다.

| 검증 | 비교량과 조건 |
|---|---|
| Lambertian 무렌즈 panel | Φs=πAsL, I(θ)=AsLcosθ, 단위 변환 |
| Lossless 직접 수송 | Raw power ledger closure; ROI 밖 power 포함 |
| Fresnel·plate | 정상 입사/경사 입사의 single-interface T; 다중 Jones와 scalar-average 차이 구분 |
| RI 적분 | 간단한 pupil의 projected solid angle; 제한된 thin/slow fixture에서만 cos⁴ |
| Forward/reverse closure | 동일 source mask·radiance·표면·loss profile에서 ray bin power/area와 reverse L integral 비교 |
| Geometric warp | affine map의 $E=Q/\lvert\det J\rvert$와 총전력 보존; finite pixel integration |
| Wave/radiometry 접속 | 동일 source illumination에서 Q와 ∫Eh의 관계, PSF window 누락 및 측도 변환 |
| Far-field | ΣIBΔΩB=Φexit, direction-cosine·spherical chart 변환 후 동일 bin power |
| RGB | 공통 좌표·lateral color 유지, source weight 중복 없는지 |
| Thermal/grid | 고정 radiance/좌표/표시 scale에서 전력과 nominal map drift 보존 |

étendue 보존은 매끄러운 4D position-angle ray bundle에서 확인하는 추가 fixture다. 단일 point의 2D pupil로 전체 4D 보존을 검증했다고 하지 않는다. 모든 데이터 sample의 첫 실행을 4D 검증으로 막을 필요도 없다.

아직 제품 호환성 또는 모델 선택으로 남아 있는 항목은 명시적으로 기록한다.

- `RI compatibility`: 제품 Rimmer 확장 알고리즘의 세부 sampling·apodization·경계 처리까지 공개식만으로 bitwise 동일하다는 근거는 없다. 물리적 reverse integral, 동일 설정에서의 수렴, 향후 reference comparison을 분리한다.
- `Physical PSF applicability`: 일반 강수차·angle/pupil caustic·높은 NA에서 scalar/local FFT가 항상 유효하지 않다. 상태는 `unsupported_by_wave_profile`일 수 있으며 ray-power 계산 가능 여부와 별도로 판정한다. 이 경우 geometric angular output을 diffraction output으로 바꾸어 표기하지 않는다.
- `Polarization`: 무편광 표면 평균과 두 편광 상태 수송은 다르다. Vector PSF와 intensity-weighted scalar PSF도 다르다. 선택을 metadata에 남긴다.
- `Absolute flux`: radiance/총광량, 실제 coating, 내부 absorption 정보가 없으면 상대 및 이상화 출력까지만 확정된다. 사용자 제공 물성이 충분하면 이 제한을 새 사용자 질문으로 만들 필요는 없다.
- `Sampling target`: 허용할 RI/각도/영상 오차는 작업 목적별 요구다. 보편적인 합격 숫자를 발명하지 않고 refinement 변화·미해결 전력·coverage를 보고한다.
- `Reciprocity`: exact ray path reversal은 유효한 검증이지만, 역방향 카메라 MTF/RI가 정방향 projector far-field와 자동으로 같지는 않다. 조명·공액·좌표·측도 및 주파수 변환이 맞는 경우만 별도의 대응 검증으로 사용한다.

---

<a id="part-6"></a>

## 6. 균일 온도 성능과 공통 재초점

### 6. 고정 초점과 공통 재초점

#### 6.1 고정 IMA 재계산

각 온도에서 `n → geometry → vertices/gaps → apertures → ray aiming → real trace → wavefront → FFT` 순서로 계산한다. 기존 온도의 chief/marginal ray, pupil mapping, OPD, diffraction reference를 값만 복사하지 않는다. 같은 물리 field·G 파장·주파수 배열·조명 프로파일을 사용한다. Diffraction-limit 비교도 그 온도의 aperture/NA/profile에서 다시 계산하고, baseline reference를 함께 표시하려면 별도 곡선으로 남긴다.

#### 6.2 Geometric RMS 기준의 공통 최적 평면

고정된 ray 집합과 가중치를 가진 균질 image space에서, 기준 평면의 교점을 $\mathbf r_{fj}$, 출사 방향을 $\mathbf d_{fj}$라 한다. Sensor가 z 방향으로 $\delta z$ 이동하면 exact straight propagation은:

$$
\mathbf r_{fj}(\delta z)=\mathbf r_{fj}+\mathbf q_{fj}\delta z,\qquad
\mathbf q_{fj}=\frac{(d_x,d_y)_{fj}}{d_{z,fj}}.
$$

이는 paraxial $\theta$ 근사를 쓰지 않는다. Field마다 normalized ray weight $\sum_jw_{fj}=1$, 지정 field weight $\eta_f\ge0$를 사용하고 $\bar{\mathbf r}_f,\bar{\mathbf q}_f$를 같은 가중 평균으로 계산한다. 각 field centroid에 대한 RMS²의 가중 합을 최소화하면:

$$
J(\delta z)=\sum_f\eta_f\sum_jw_{fj}
\left|\mathbf r_{fj}-\bar{\mathbf r}_f+
(\mathbf q_{fj}-\bar{\mathbf q}_f)\delta z\right|^2,
$$

$$
\delta z_*=-\frac{
\sum_f\eta_f\sum_jw_{fj}(\mathbf r_{fj}-\bar{\mathbf r}_f)\cdot(\mathbf q_{fj}-\bar{\mathbf q}_f)
}{\sum_f\eta_f\sum_jw_{fj}|\mathbf q_{fj}-\bar{\mathbf q}_f|^2}.
$$

분모가 0이면 이 criterion으로 유한 최적 평면을 결정할 수 없다. 이동 범위가 제한되면 이 convex quadratic에서 $\delta z_*$를 허용 구간에 제한한다. 이동 중 유효 ray 집합/가중치가 바뀌거나 다른 aperture/surface를 통과하면 이 닫힌식의 조건이 깨지므로 실제 엔진의 목적함수를 다시 평가한다.

이 focus는 **geometric centroid RMS** 기준이다. RMS wavefront, sagittal/tangential line focus, MTF 최적 초점과 같다고 label하지 않는다. 모든 field에서 별도 최적 $\delta z_f$를 계산한 뒤 합치는 것은 field-curvature 진단이지 하나의 평면 센서 성능이 아니다.

#### 6.3 FFT MTF 기준 공통 재초점

명시된 field·방향·물리 주파수 $\nu_k$의 weight로 예를 들어 다음 목적함수를 정의할 수 있다. 이는 설계자가 선택하는 merit의 예이며 보편적인 “최적 focus” 정의가 아니다.

$$
\delta z_* = \underset{\delta z\in[z_{\min},z_{\max}]}{\operatorname{argmax}}
\sum_f\eta_f\sum_{a\in\{T,S\}}\sum_k\omega_{fak}
\mathrm{MTF}_{fa}(\nu_k;\delta z).
$$

실행은 거친 공통 z 스윕으로 여러 유효 극값 후보를 찾고, 후보 bracket을 일차원 bounded search로 세분한 뒤 같은 최종 sampling으로 비교한다. MTF merit가 단봉이라고 가정하지 않는다. 후보마다 reference wavefront와 image-coordinate scale을 일관되게 다시 계산한다. IMA만 이동하는 sweep에서는 회절 장의 canonical Π와 그 위 경계장 진폭을 유지하는 `hold_boundary_field_amplitude`를 기본으로 한다. Reference가 바뀐 뒤 B=1을 다시 설정하면 조명을 재정의하는 별도 metric이므로 혼합하지 않는다. 120 lp/mm 한 점만을 최적화하면 중간 주파수 영점/contrast reversal을 놓칠 수 있으므로 전체 MTF와 complex OTF도 보존한다.

일반 FFT MTF를 z마다 반복한 sweep은 독립 구현의 정확한 정의다. OpticStudio `FFT Through Focus MTF`는 빠른 별도 경로를 사용해 일반 FFT MTF와 차이가 날 수 있으므로 그 창의 수치와 자동 동치라고 부르지 않는다. [TH03](#ref-TH03)

Afocal output에서는 관찰 평면의 z 이동만으로 angular blur가 refocus되지 않는다. 유한 source의 axial 위치, 특정 lens group 위치 등 허용 보상변수 $c$를 정하고 angular PSF/MTF 목적함수 $J(c)$를 최적화한다. 허용 보상이 없으면 `refocus_not_applicable`로 기록한다. OpticStudio의 afocal through-focus 창은 defocus 축을 diopter 단위의 optical power로 정의한다. 그 수학적 defocus sweep은 실제 source/lens 이동 보상과 별개로 label하며, 이동량 mm와 변환 없이 동일시하지 않는다. [TH03](#ref-TH03)

### 7. 온도별 출력과 thermal mapping drift

각 온도에 최소한 고정 IMA spot/FFT MTF, 선택한 common focus와 이동량, refocus 후 spot/FFT MTF를 저장한다. 본래 요청한 RI·영상/그리드·각도 성능을 해당 온도에 함께 계산하는 경우도 같은 재계산 처방을 사용한다. 온도마다 기준 G 파장·source field 좌표·정규화 분모·표시 exposure를 바꾸지 않는다.

Chief-ray mapping을 $\mathbf g_T(\mathbf q)$라 하면 baseline에 대한 위치 변화를 별도로 저장한다.

$$
\Delta\mathbf g_T(\mathbf q)=\mathbf g_T(\mathbf q)-\mathbf g_{T_0}(\mathbf q).
$$

Finite image에서는 mm, afocal에서는 명시된 angular chart의 rad로 저장한다. 온도별로 다시 fit한 paraxial/linear reference에 대한 distortion과 $\Delta\mathbf g_T$를 구분한다. 예를 들어 $g_T(q)=1.01g_{T_0}(q)$이면 온도별 배율을 다시 맞춘 distortion은 동일할 수 있지만 실제 image scale은 1% 바뀐다. 이 예는 좌표 정의의 반례이며 특정 렌즈의 열 예측값이 아니다.

권장 결과 항목:

- `temperature`, `wavelength_basis`, `index_model_id`, `geometry_scale_by_lens`, `global_vertices`, `physical_stop`.
- `fixed_ima_metrics`, `focus_criterion`, `common_focus_offset`, `refocused_metrics`.
- `chief_map_nominal_reference`, `scale_or_fov_change`, `intrinsic_distortion_at_temperature`.
- `raw_transmitted_power`, 광량 손실/미해결 장부, 표본·수렴·모델 적용성 상태.

### 8. 검증과 아직 결정해야 할 모델 선택

#### 8.1 실제 구현 경로에서 수행할 시험

| 시험 | 확인할 관계 | 한계 |
|---|---|---|
| $T=T_0$ | 기준 형상·vertex·index·결과의 재현 | fit 기준 $T_c$와 형상 기준 $T_0$를 혼동하지 않아야 함 |
| 온도 순서 변경 | 동일 온도의 결과가 sweep 순서와 무관 | 누적 갱신·잘못된 cache를 검출 |
| Index-only | 형상/STOP/좌표 고정, n만 온도 갱신 | 실제 전체 thermal 결과가 아닌 원인 진단 |
| Geometry-only | n/파장 고정, 형상·위치 정책만 갱신 | 같은 제한 |
| Coupled | 둘 모두 갱신하여 재추적 | 두 진단 변화량의 단순 합과 같을 필요 없음 |
| Sag similarity | 임의 유효 r에서 $z_T(r)=s z_0(r/s)$ | 수치식·단위·A2…A12 scaling 확인 |
| Absolute/relative 변환 | 동일 Snell ratio와 $n\ell/\lambda$ 위상 | 정확히 같은 adapter 출력끼리 비교 |
| Transmission 길이 | $\tau(\ell_1+\ell_2)=\tau(\ell_1)\tau(\ell_2)$ | 균질·일정 μ 모델에 한함 |
| 공통 geometric focus | 닫힌 RMS 최적평면과 직접 z 스윕 일치 | ray/weight 집합이 일정한 fixture |
| Nominal mapping | 단순 확대의 scale drift는 남고 재fit distortion은 동일 가능 | 왜곡률만으로 열 안정성 판단하는 오류 검출 |

독립적인 **합성 scaling 시험**에서는 모든 길이와 진공 파장을 같은 s배, n과 대응 조명을 고정했을 때 정규화 2D PSF가 $h_s(\mathbf x)=s^{-2}h_0(\mathbf x/s)$, OTF가 $\mathrm{OTF}_s(\mathbf f)=\mathrm{OTF}_0(s\mathbf f)$로 바뀌는지 검사한다. 이는 Fourier 단위와 ray–wave 연결을 확인하는 시험이다. 실제 thermal sweep은 파장이 고정되고 재료별 n/s가 다르므로 이 similarity를 열응답 공식으로 사용하지 않는다.

#### 8.2 논쟁·미확정 사항의 구분

| 항목 | 현재 판단 | 엔진에 남길 처리 |
|---|---|---|
| SCHOTT vacuum λ와 OpticStudio reference-air λ | 공개 자료의 규약이 다름; 수학적 모순이라기보다 모델/호환성 선택 | 별도 profile, 혼합 금지, 실제 제품 비교는 미검증 |
| Constant CTE vs temperature-dependent CTE | 어느 것도 무조건 우월한 보편식이 아님; 제공 데이터 정의에 달림 | 정의별 s(T) adapter |
| Lens center/front vertex/air gap 고정 | 모든 경우에 공통인 유일한 답은 없음 | 독립 위치 규칙 한 세트만 선택; 기본 center-fixed를 쓸 경우 이를 표기 |
| STOP/F-number를 온도별 재조절 | 물리 hardware sweep과 다른 실험 | 고정 hardware 결과와 분리 |
| “최적 focus” | 기준에 따라 다름 | geometric/MTF 등 목적함수와 common compensator 기록 |
| 고온에서 μ 불변 | 제공 자료가 없을 때의 명시적 근사 | throughput 결과의 모델 상태로 남김 |
| PMMA/PC의 고온 비선형·시간 의존성 | 입력 grade의 균일 평형 자료 범위 밖을 보장할 수 없음 | 물성 유효범위 확인; 범위 밖을 확정값으로 외삽하지 않음 |

제조 공차·수율, 비균일 온도, 발열원/열전도, 공간적 n 변화, 응력 복굴절, mount 응력/변형은 이번 엔진의 제외 범위다. 이들을 계산하지 않는다고 해서 명시한 균일 모델의 초점·spot·MTF 계산이 불가능해지는 것은 아니다. 해당 효과가 포함된 실물과의 일치 보장은 별도 문제다.

---

<a id="part-7"></a>

## 7. 결과·검증·미결 쟁점과 피드백 반영

### V.1 결과를 만드는 최소 단위

하나의 결과는 숫자만이 아니라 **처방 상태 + 정의 + 값 + 상태 + 재현 정보**다. 최소 실행 키는 `(design_revision, material_revision, temperature, field, wavelength, focus_state, physics_profile, numerics_profile)`다. 결과 배열의 축 이름·단위·순서와 좌표 원점을 함께 보존한다.

| 출력 | 값과 함께 저장할 정의 | 혼동하지 않을 출력 |
|---|---|---|
| G FFT MTF | 전체 f축, complex OTF, T/S 방향, 기준면, pupil/OPD profile, 유효 mask, 표본·창, 0주파수 정규화 분모 | Source pixel MTF, square-wave response, 각 field별 best-focus envelope |
| Diffraction reference | 같은 pupil mask·amplitude에서 W=0인지, 이상 원형 pupil인지 | 모든 처방의 유일한 보편 상한곡선 |
| Spot | 좌표·가중치·chief/centroid reference, RMS/GEO, 차단·미해결 비율 | 회절 PSF, 별도 focus에서 만든 spot |
| Field curvature | chief 주변 differential T/S focus, 기준 IMA·분석축·perturbation 수렴 | Full pupil의 RMS best-focus map |
| 수차 | ray fan 또는 OPD basis·indexing·normalization·제거항·잔차 | 정의 없는 coma 단일 숫자, 전체 고차 효과를 Seidel 다섯 값으로 대체 |
| EFL/BFL | first-order state·index 기준, 마지막 굴절면/IMA 원점, 파장 | 유한 공액 최적 IMA, afocal mode flag만으로 만든 무한대 값 |
| Distortion/grid | 실제 chief mapping, reference mapping, scale/FOV, 공통 IMA, 국소 blur 별도 | 기준을 매번 재피팅해서 제거된 scale 변화 |
| 일반 RI | 동일한 균일 장면, raw irradiance, field domain, peak 또는 명시한 분모 | 유한 발광면의 edge irradiance, source 전체 far-field |
| Far-field | raw bin flux, 정확한 bin solid angle, W/sr 또는 상대 intensity, 영역 안/밖 전력 | 단순 ray count, 면적 조도, angular MTF |
| Image simulation | 선형 source 값, source/Q/h 측도, field 보간, 파장·채널 대응, 노출 | RGB display code의 단순 평균, PSF shape에 throughput 이중 반영 |
| 온도 변화 | T₀ 기준 형상, n(T), STOP·anchor 정책, 고정/공통 보정 상태, 보정 이동량 | Mount FEA가 반영된 실제 조립품 전응답 |

#### V.1.1 광선과 파면의 중간 기록

최종 그래프만 저장하면 부호·기준면·광량 오류를 추적하기 어렵다. 대표 field의 원 광선 기록에는 최소 다음이 있어야 한다.

```text
ray_id, field_id, wavelength_id, launch_position, launch_direction,
launch_phase_opl, sampling_measure, sampling_weight_or_pdf,
per_surface: hit_position, incident/transmitted_direction,
             incident/transmitted_medium, path_length, cumulative_opl,
             aperture_status, transmission, cumulative_power,
stop_target, stop_residual, intersection_residual, trace_status

wave_record: boundary_plane, boundary_coordinates, mapping_jacobian,
             phase_reference, reference_point_or_plane, signed_extension,
             opd_length, amplitude_definition, illumination_profile,
             Fourier_chart, cell_measure, mask,
             output_plane, output_scale, fft_sign_and_normalization
```

모든 설계의 모든 ray를 영구 저장해야 한다는 뜻은 아니다. 검증 사례와 실패 사례에는 이를 보존하고, 대량 생성에서는 재현 가능한 표본과 집계 진단을 남길 수 있다. 원자료를 줄이는 정책도 버전으로 저장한다.

### V.2 상태와 에너지 장부

처방 전체와 metric별 상태를 구분한다. 예를 들어 afocal에서 일반 RI의 결상 가정이 성립하지 않아도 ray trace·angular spot·EFL은 유효할 수 있다.

권장 상태는 `valid`, `missing_input`, `invalid_geometry`, `no_transmitted_path`, `not_applicable`, `unsupported_profile`, `numerical_failure`, `nonconverged`, `model_applicability_unverified`, `outside_material_range`다. 일부는 배타적인 상태가 아니라 주상태와 진단 flag 조합으로 저장해도 된다.

- 물리적 aperture 차단·흡수·reflection loss·TIR에 의한 direct-path 소멸과 solver 실패를 구분한다.
- 완전히 차단된 field가 양의 전역 RI 분모를 갖는다면 RI=0은 정상이다. PSF의 총 전력이 0인데 이를 1로 정규화하는 것은 정의되지 않는다.
- `NaN`을 정상적인 0으로 바꾸거나, 미해결 ray를 단순 차단으로 집계하지 않는다.
- 명시한 reporting window 밖으로 나간 전력은 “미해결”과 별도다. 작은 PSF 창을 다시 1로 정규화하면 잘린 에너지가 숨는다.

동일한 입사 또는 방출 기준에서 직접 투과 경로의 장부를 구성한다. 시뮬레이션이 추적하지 않는 반사 경로는 해당 면에서 reflection loss로 기록하며 이후 경로를 중복 가산하지 않는다.

$$
\Phi_{\mathrm{launched}}
=\Phi_{\mathrm{reported}}+\Phi_{\mathrm{outside\ domain}}
+\Phi_{\mathrm{blocked}}+\Phi_{\mathrm{absorbed}}
+\Phi_{\mathrm{reflected/direct\ path\ rejected}}
+\Phi_{\mathrm{unresolved}}+\varepsilon_{\mathrm{balance}}.
$$

표본화되지 않은 source 방향을 사용한 importance sampling에서는 launched의 의미를 PDF로 보정된 전체 분포 추정량으로 정의해야 한다. Aim 실패 시 원인상 실제로 막힌 광선임을 입증하지 못했다면 미해결로 남긴다. 수치적 ε와 Monte Carlo 표본 오차는 같은 물리 손실이 아니다.

### V.3 필요한 검증 사례와 검증되는 범위

다음은 **후속 구현의 시험 명세**다. 본 문서 작성으로 통과한 시험 목록이 아니다. 실제 production API와 같은 교점·위상·FFT·radiometry 경로를 호출해야 한다.

| 시험 | 기대 관계·검증 대상 | 이것만으로 검증되지 않는 것 |
|---|---|---|
| 평면·구면·평행판 | 독립 해석 교점, 벡터 Snell, 매질 순서, transmitted root, OPL | 일반 비구면의 최초 교점 선택 |
| A₂와 작은 aperture | c_vertex=c+2A₂의 first-order power와 실광선의 작은 aperture 극한 | 큰 aperture 회절/수차 전체 |
| Surface gradient·Hessian | 독립 차분/자동미분 및 step 수렴, axis limit | 잘못된 물리 branch를 선택하는 문제 |
| STOP aiming | 목표 좌표 잔차, perturbation/Jacobian, 물리 clear mask | Vendor Real/Paraxial Ray Aiming의 옵션 전체 |
| 개별 ray 반전 | 원 교점·방향·OPL 복원, 동일 path attenuation | 정방향/역방향 카메라 MTF·RI 자동 동일성 |
| Phase-gradient | 매끄러운 경계면에서 ∇S=n d, reference를 뺀 OPD 미분 관계 | Caustic·clipping 불연속에서 단일 branch 근사 |
| Piston·tilt·defocus | Piston의 MTF 불변, 알려진 tilt의 PSF 이동과 OTF phase, 고정 IMA defocus 보존 | 임의 best-fit 구면을 빼도 된다는 주장 |
| 원형·비대칭 pupil | 원형 해석 MTF 및 비대칭 방향성, 물리 frequency scale | 실제 렌즈의 boundary field 정확도 |
| FFT 대 직접 적분 | 선택한 같은 scalar boundary의 적분·표본 경로 차이 | 공유한 ray/OPD 앞단 오류나 vector diffraction |
| Geometry+λ similarity | 모든 길이·λ를 s배, n 고정이면 h_s(x)=s⁻²h₀(x/s), OTF_s(f)=OTF₀(sf) | 실제 온도 응답: 실제 λ는 고정이며 재료별 팽창·n(T)가 다름 |
| Radiometry 적분 | 같은 장면에서 forward flux/area와 reverse radiance 적분의 폐합 | 서로 다른 uniform-scene와 finite-edge 모델의 동일성 |
| Étendue 선택 사례 | Regular 4D 위치·방향 bundle의 n² cosθ dA dΩ 보존 | 한 물점의 2D pupil 면적만으로 étendue 판정 |
| 차단·0·미해결 | 정의된 0/null/status, 전력 장부 | 잘린 pupil이 설계적으로 항상 나쁘다는 판정 |
| 열 identity | ΔT=0, n/CTE 변화 없음, 온도 실행 순서 불변 | 실제 mount 위치 규칙 또는 물성 유효 범위 밖 |
| 영상 에너지 | ∫E dA=∫Q dA_s, 공통 좌표 RGB, throughput 한 번 반영 | h를 단위적분한 뒤의 이 등식은 정규화상 항등식이므로 경계장 진폭의 물리 정확도 검증은 아님. 정규화 전 scalar flux와 ray flux 비교가 별도 필요 |

#### V.3.1 표본 수렴과 모델 적용성은 별도다

최소한 다음을 독립적으로 변화시킨다.

1. Ray/경계장 표본 밀도와 위상·진폭 보간 해상도.
2. Aperture edge 적분, offset grid 또는 서로 배수가 아닌 표본 격자.
3. FFT padding, PSF 물리 창 크기, 출력 grid 해상도.
4. Field 보간 밀도와 clipping 시작점 부근 표본.
5. RI 방향 격자, finite-source 면적·방향 적분 및 angular bin 크기.
6. Differential ray/aiming의 perturbation step 또는 미분 전달 정확도.

Padding만 늘려 곡선이 부드러워진 것을 원 pupil 정보의 수렴으로 세지 않는다. 표본을 늘려도 동일한 잘못된 branch나 부적절한 물리 모델에 수렴할 수 있다. 수치 변화량·보존 잔차와 scalar/FFT 적용성 판정을 따로 보고한다. 보편적인 F/#·NA·field 각도의 단일 합격선을 문서만으로 만들어 내지 않는다.

#### V.3.2 비교 목적에 따라 합격 기준을 분리한다

- **수치 계산 검증:** 선언한 수식과 알고리즘을 얼마나 정확하게 계산했는가.
- **물리 모델 검증:** 선택한 ray-derived scalar 모델이 해당 처방·영역에서 충분한가.
- **설계 적합성:** 필요한 MTF·광량·시야·온도 목표에 도달하는가.
- **제품 호환:** 특정 Zemax 버전·분석창·설정의 결과와 얼마나 일치하는가.

사용 목적의 오차 예산이 아직 없으면 수렴 변화량과 근거는 반환할 수 있지만 임의의 “정확도 99%” 또는 설계 합격값을 붙이지 않는다. 일반 FFT를 여러 IMA에서 반복하는 초점 스윕과 공식 FFT Through Focus MTF의 빠른 경로는 다를 수 있다.[G02](#ref-G02)

### V.4 현재 남는 쟁점과 구현 시 처리

논쟁적이라는 말은 알려진 물리식 자체가 미확립이라는 뜻과, 모델·convention 선택이 여러 개라는 뜻을 구분해야 한다. 다음은 구현자가 알아야 할 실제 경계다.

| ID | 쟁점 | 본 명세의 처리 | 남는 확인 / 영향 |
|---|---|---|---|
| O01 | Zemax FFT의 전체 cosine-space mapping·가중치·OPD 옵션 재현 | 공개 scalar 이론으로 닫힌 독립 FFT profile을 제공. 결과에 profile ID 기록 | 실제 같은 처방/설정의 제품 출력 대조 전에는 `Zemax exact`로 표시할 수 없음 |
| O02 | 큰 NA·강한 pupil distortion·caustic의 scalar FFT 적용 | 경계장 단일 branch, local expansion 오차, 직접 적분 비교로 적용성을 검사 | 실패 시 해당 wave metric 미지원/비수렴. 조용히 geometric MTF로 대체하지 않음 |
| O03 | 균일 lens-metric pupil과 실제 Lambertian source pupil | Geometry/OPL 공유, amplitude·measure 분리. 조건부 PSF 재사용 | Q 총량만으로 pupil 내 가중치 차이를 보정할 수 없음 |
| O04 | FFT와 Huygens 중 image kernel 방법 | 검증된 FFT가 적용되면 사용 가능. 필요한 경우 별도 scalar 적분/profile 추가 | 독립 image synthesis와 제품 Image Simulation의 모드 일치를 구분 |
| O05 | 원거리 angular 성능의 단위·좌표·대역 | Direction/solid-angle chart와 Jacobian 명시. Angular MTF 범위는 별도 입력 | 120 image lp/mm를 자동으로 같은 숫자의 cycles/degree로 전용하지 않음 |
| O06 | Afocal longitudinal/field curvature의 제품 diopter convention | 독립 local vergence 또는 요청한 보정 자유도를 명시한 결과만 반환 | 제품 부호·reference·옵션을 대조하기 전에는 같은 label로 확정하지 않음 |
| O07 | EPD/F/#와 물리 STOP 변환·제품 ray aiming | 지원한 정의에 한해 T₀ STOP resolve, 실광선 목표 맞춤 | 모든 vendor aiming 옵션·자동 semi-diameter solve는 별도 지원/대조 |
| O08 | 재료의 공기 기준과 vacuum 기준, thermal 식 인수 | Adapter별 모델 ID와 기준 λ/T/P를 보존, 공통 n_abs로 변환 | 제조사 thermal 식과 AGF 식을 파라미터 이름이 같다는 이유로 동일 취급하지 않음 |
| O09 | 균일 팽창 후 렌즈·STOP·OBJ·IMA 위치 | 명시한 anchor/spacing 정책. Fixed IMA와 마지막 gap을 중복 강제하지 않음 | 선택한 이상화 조립 조건의 결과이며 실제 mount 예측은 아님 |
| O10 | Distortion 기준과 열 scale/FOV drift | Intrinsic distortion과 nominal-referenced chief map을 모두 보존 | 제품 작은-field reference·scan 규칙의 옵션별 대조 필요 |
| O11 | Field curvature·FFT의 T/S 방향 | 분석별 축 규약을 별도 저장 | 임의 2D field에 radial/tangential 회전을 적용한 것을 vendor axis label과 혼동하지 않음 |
| O12 | 전력 손실·편광·coating phase | 명시한 scalar/unpolarized/코팅 profile만 사용 | 연속 계면에서 편광 상태 보존이나 복소 coating phase가 필요한 조건은 별도 모델 |
| O13 | 유한 source 경계와 일반 RI | Uniform-scene RI와 finite-emitter image를 별도 계산 | Edge blur, angular map fold는 단순 field별 RI 곡선으로 대표할 수 없음 |
| O14 | 선명도만 높은 설계의 평가 | 동일 field·aperture·초점, 전체 OTF와 raw flux를 함께 저장 | 최소 광량·화질 목표 없이는 임의 종합 순위/합격을 만들지 않음 |
| O15 | Canonical 경계평면과 초점 스윕의 조명 유지 | Π 선택 및 진폭 정책을 profile state로 고정. 기본은 같은 Π의 A를 유지하고 reference만 갱신 | Π를 임의 이동하거나 매번 B=1로 초기화하면 다른 조명/metric이 될 수 있음 |

O01은 공개 이론 기반 계산의 불가능을 뜻하지 않는다. **사용자가 원하는 제품 참조와 독립 알고리즘 사이의 동등성이 아직 실증되지 않았다는 뜻**이다. O08/O09는 실제 재료 데이터와 위치 정책이 제공되면 계산을 결정할 수 있는 입력·모델 문제다. O02/O12는 선언한 물리 모델의 범위 문제다. 이 세 종류를 같은 “공식이 없다”로 묶으면 안 된다.

### V.5 이전 피드백의 반영 내역

| 검토에서 나온 지적 | 본 문서 반영 |
|---|---|
| Sag와 A₂ power, thermal asphere scaling | 전체 항의 gradient/Hessian 및 c_vertex, 온도별 계수 차원 규칙 |
| Newton 수렴만으로 최초 물리 교점 보장 불가 | Branch·전방 해·형상 허용영역·접선/다중근·실패 상태 |
| Raw OPL 차이를 OPD로 사용하면 안 됨 | 공통 launch phase, 동일 출력 경계, 명시 reference, signed extension |
| P099/P100은 정책 이름만 존재 | 좌표·진폭·측도·Fourier scale을 닫는 독립 profile과 남는 vendor 동등성 구분 |
| 렌즈 PSF에 Q만 곱해 광원 영상 재사용 불가 | `lens_metric` / `source_transport`와 재사용 조건 |
| 매 field 최적 초점과 실제 센서의 혼동 | 고정 공통 IMA·공통 보정 자유도, 진단용 best-focus map 분리 |
| MTF magnitude의 spurious resolution | 전체 곡선·complex OTF·registration 저장, 임의 단조화 금지 |
| 온도별 재피팅이 scale 변화 은폐 | 기준 T₀ chief map·FOV·배율 변화와 intrinsic distortion 동시 저장 |
| Reverse camera가 자동 등가가 아님 | Ray 경로·raw transport·조건부 국소 성능 비교를 분리 |
| 일반 RI, 유한 source 영상, far-field 혼동 | 각각의 적분·측도·정규화·적용 조건을 별도로 정의 |
| 전력·진폭·n²·Jacobian 중복 | Ray packet, radiance, boundary amplitude, Q/h의 책임 분리 |
| Focal 개발만 남아 afocal 범위 축소 우려 | Afocal plane OPD·angular PSF/MTF·far-field의 별도 수식 및 조기 검증 |
| 111개 항목이 실행 schema는 아님 | 다섯 구조의 입력/출력 책임과 단위·상태·버전 규약 |
| 이론 검산을 실제 엔진 정확도로 승격 | 계획 시험과 실행 증거, 수치·모델·설계·제품 검증의 구분 |

### V.6 후속 구현에서 먼저 완성할 한 경로

첫 검증 사례는 제공된 실제 처방 하나의 고정 IMA, G 파장, 축상과 비축 field를 대상으로 한다. 입력 정규화부터 면별 ray, OPD, FFT PSF/MTF, 원 전력과 진단 기록까지 한 번에 재현해야 한다. 그 다음 finite-object collimator의 angular 계산과 유한 source far-field를 연결한다. 이 순서는 검증을 위한 것으로, field curvature·수차·image/grid·열 성능을 최종 범위에서 제외하지 않는다.

실제 처방값·재료 payload·기준 결과가 없는 곳에 설명용 예제를 채워 넣고 “사용자 광학계 결과”로 기록하지 않는다. Library/출처 URL의 개수도 시험 사례 수로 세지 않는다. 실제 엔진 release, 본 문서 revision, 입력 schema revision, 물리 profile revision을 각각 기록해야 한다.

---

<a id="part-8"></a>

## 8. 통합 참고문헌과 근거의 범위

확인일: 2026-10-18. Ansys 문서는 제품 정의, 제조사 문서는 물성식·기준 조건, 원저자 논문·교재 및 공개 source는 수학·알고리즘의 근거다. 제품 옵션의 설명은 독립 구현의 수치 검증을 뜻하지 않는다. 본문에 표시한 Jacobian, anchor, ray 미분, FFT 연결, Q/h 분해 및 실패 정책은 공개 정의를 바탕으로 한 직접 유도·독립 구현 규약이다. 여러 release의 문서를 사용했으므로 특정 한 build의 모든 기능과 동일하다는 주장을 하지 않는다.

기존 S-ID를 유지했다. W01은 공개 RS Fourier 유도, W02는 제품 표본 설명, TH는 물성·열, G는 추가 기하 분석 정의다. 동일 자료를 다른 장에서 인용한 ID가 있으면 그 관계도 보존한다.

<a id="ref-S01"></a>

### S01

**S01 — Ansys, Even Asphere**. 표면 sag·계수 정의. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v261/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Even_Asphere.html

<a id="ref-S02"></a>

### S02

**S02 — Ansys, Standard**. Standard의 plane/sphere/conic 범위. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Standard.html

<a id="ref-S03"></a>

### S03

**S03 — Michael J. Hayford, ray-optics profiles implementation**. implicit intersection·surface gradient의 공개 구현. 안정적 first-hit 규약 전체를 보장하는 근거는 아니다. https://ray-optics.readthedocs.io/en/stable/_modules/rayoptics/elem/profiles.html

<a id="ref-S04"></a>

### S04

**S04 — Ansys, How to use Ray Aiming**. 실제 STOP 목표를 반복하여 맞추는 제품 방식, paraxial/real stop sizing 구분. https://optics.ansys.com/hc/en-us/articles/42661778056083-How-to-use-Ray-Aiming

<a id="ref-S05"></a>

### S05

**S05 — Ansys, Paraxial vs. Real pupils in optical system**. pupil/chief/marginal 의미. https://optics.ansys.com/hc/en-us/articles/42661958583571-Paraxial-vs-Real-pupils-in-optical-system

<a id="ref-S06"></a>

### S06

**S06 — PBRT 4e, Specular Reflection and Transmission**. 원저자의 vector refraction·Fresnel 정의/구현. https://pbr-book.org/4ed/Reflection_Models/Specular_Reflection_and_Transmission

<a id="ref-S07"></a>

### S07

**S07 — Ansys, Effective Focal Length**. EFL과 image-space secondary focal length 구별. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Effective_Focal_Length.html

<a id="ref-S08"></a>

### S08

**S08 — Ansys, Back Focal Length**. 비공기 image space의 IMA 원점 예외. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Back_Focal_Length.html

<a id="ref-S09"></a>

### S09

**S09 — Michael J. Hayford, ray-optics first-order implementation**. paraxial power·cardinal 계산의 공개 구현. https://ray-optics.readthedocs.io/en/stable/_modules/rayoptics/parax/firstorder.html

<a id="ref-S10"></a>

### S10

**S10 — Ansys, Standard Spot Diagram**. geometric spot·reference와 표본 의존 지표. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Standard_Spot_Diagram.html

<a id="ref-S11"></a>

### S11

**S11 — Ansys, Field Curvature and Distortion**. parabasal field curvature, scan 축, small-field chief distortion reference. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Field_Curvature_and_Distortion.html

<a id="ref-S12"></a>

### S12

**S12 — Ansys, Grid Distortion**. grid chief mapping. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Grid_Distortion.html

<a id="ref-S13"></a>

### S13

**S13 — Ansys, Seidel Coefficients**. 저차 수차 descriptor의 범위. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v261/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Seidel_Coefficients.html

<a id="ref-S14"></a>

### S14

**S14 — TU Delft, Advanced Optical Imaging: Wavefield aberrations**. wavefront basis 수차 표현. https://qiweb.tudelft.nl/aoi/wavefieldaberrations/wavefieldaberrations/

예시 표의 단일 index·signed-m 표기를 변환 없이 전사하지 않는다. G9.1은 명시한 일반 다항식에서 tuple basis를 직접 생성하며 R₃¹=3ρ³−2ρ를 사용한다.

<a id="ref-S15"></a>

### S15

**[S15](#ref-S15) MIT OCW, 2.71/2.710 Optics Lecture 22, 2009.** Incoherent PSF/OTF와 pupil autocorrelation의 강의 원자료. <https://ocw.mit.edu/courses/2-71-optics-spring-2009/bebe78c6a2ae2603a0d1a3ea85e0b4ca_MIT2_71S09_lec22.pdf>

<a id="ref-S20"></a>

### S20

**S20** — Ansys, *Relative Illumination*, Release 2026 R1.02. RI 정규화·장면·적분·적용 조건의 제품 정의. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Relative_Illumination.html>

<a id="ref-S21"></a>

### S21

**S21** — Pharr, Jakob, Humphreys, *PBRT 3e, Working with Radiometric Integrals*. Radiance–irradiance, projected solid angle 및 면적-입체각 변환의 원저자 설명. <https://www.pbr-book.org/3ed-2018/Color_and_Radiometry/Working_with_Radiometric_Integrals>

<a id="ref-S23"></a>

### S23

**S23** — Ansys / Mark Nicholson, *How to simulate high resolution images*. Field별 PSF와 영상 합성 구조의 제품 튜토리얼. <https://optics.ansys.com/hc/en-us/articles/42661994652819-How-to-simulate-high-resolution-images>

<a id="ref-S27"></a>

### S27

**[S27](#ref-S27) Ansys, OpticStudio 2025 R1, TCE.** Linear TCE와 ppm/K 단위. [공식 문서](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/TCE.html)

<a id="ref-S29"></a>

### S29

**[S29](#ref-S29) SCHOTT, TIE-19: Temperature Coefficient of the Refractive Index, August 2023.** p.2 Eq.2–7: absolute thermal model, vacuum μm와 reference-index 근사; p.3 유효범위. [제조사 PDF](https://media.schott.com/api/public/content/63b62a0d73a640efa8c6fde92f2ad5c2?v=97a6a05a)

<a id="ref-S37"></a>

### S37

**S37 — Ansys, Modeling optics with realistic edge apertures**. clear aperture·chip/rim geometry 구분. https://optics.ansys.com/hc/en-us/articles/42661783548563-Modeling-optics-with-realistic-edge-apertures

<a id="ref-S39"></a>

### S39

**[S39](#ref-S39) Ansys, Afocal Image Space, 2025 R2.** Angular output/plane reference와 finite-object collimator의 제품 mode 정의. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Afocal_Image_Space.html>

<a id="ref-S40"></a>

### S40

**[S40](#ref-S40) Ansys, FFT MTF, 2025 R2.** Image-space frequency, T/S 명칭, pupil/cosine-space 적용 가정의 제품 설명. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html>

<a id="ref-S41"></a>

### S41

**[S41](#ref-S41) Ansys, Diffraction Limited, 2025 R2.** 제품 reference field/zero-OPD 기준. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Diffraction_Limited.html>

<a id="ref-S42"></a>

### S42

**[S42](#ref-S42) Ansys, FFT PSF, 2025 R1.** Chief-normal imaginary plane와 scalar/far-field 제한의 제품 설명. 본 profile의 독립 RS mapping 자체의 근거는 W01이다. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_PSF.html>

<a id="ref-S44"></a>

### S44

**S44 — Ansys, Aperture Type**. system aperture의 독립 입력과 종류. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Aperture_Type.html

<a id="ref-S47"></a>

### S47

**S47** — Ansys, *Detector Viewer*, Release 2025 R1.01. Irradiance와 radiant intensity의 단위·정의. NSC detector 기능이 본 Sequential 엔진과 동일하다는 근거는 아니다. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Detector_Viewer.html>

<a id="ref-S48"></a>

### S48

**S48 — Ansys, Field Type**. field angle/object height 및 finite ray aiming 의미. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Field_Type.html

<a id="ref-S49"></a>

### S49

**S49 — Ansys, Working F/#**. actual conjugate·4 marginal NA² 평균·aperture/vignetting 처리. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Working_F.html

<a id="ref-S51"></a>

### S51

**S51** — Ansys, *Polarization (System Explorer)*, Release 2025 R2. Unpolarized 처리와 intensity weighting/vector 분석 기능 구분. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Polarization_System_Explorer.html>

<a id="ref-S52"></a>

### S52

**S52** — Ansys, *Image Simulation*, Release 2025 R2. 제품 Diffraction/Image Simulation 및 표시/좌표 규약. <https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Image_Simulation.html>

<a id="ref-S54"></a>

### S54

**[S54](#ref-S54) Ansys, OpticStudio 2026 R1, Index of Refraction Computation.** System-air normalization, air 식, signed thermal parameter. [공식 문서](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v261/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Index_of_Refraction_Computation.html)

<a id="ref-S55"></a>

### S55

**[S55](#ref-S55) Ansys, OpticStudio 2025 R2, Wavelength Data.** System-air wavelength convention. [공식 문서](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Wavelength_Data.html)

<a id="ref-S57"></a>

### S57

**[S57](#ref-S57) Ansys, How OpticStudio calculates refractive index at arbitrary temperatures and pressures.** Reference-air λ 및 absolute/relative 변환 순서; thermal 식의 λ_rel을 명시. [공식 기술 문서](https://optics.ansys.com/hc/en-us/articles/42661799783443-How-OpticStudio-calculates-refractive-index-at-arbitrary-temperatures-and-pressures)

<a id="ref-S58"></a>

### S58

**[S58](#ref-S58) Ansys, Zemax AGF material catalog file format.** TD record와 reference temperature 등 import metadata. [공식 기술 문서](https://optics.ansys.com/hc/en-us/articles/46891219767827-Zemax-AGF-material-catalog-file-format)

<a id="ref-S59"></a>

### S59

**S59** — Pharr, Jakob, Humphreys, *PBRT 4e, Dielectric BSDF*, 특히 transmission의 radiance 변환. n² 보정과 power 수송의 구분에 대한 원저자 근거. <https://pbr-book.org/4ed/Reflection_Models/Dielectric_BSDF>

<a id="ref-S60"></a>

### S60

**S60 — Ansys, Sagittal and Tangential**. 분석별 object-local T/S 규약. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Sagittal_and_Tangential.html

<a id="ref-G01"></a>

### G01

**G01 — Ansys, Image Space F/#**. 무한 공액 paraxial EFL/EPD의 nominal 정의. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Image_Space_F.html

<a id="ref-G02"></a>

### G02

**G02 — Ansys, FFT Through Focus MTF**. ordinary FFT MTF와 다른 빠른 계산 경로의 제품 설명. https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_Through_Focus_MTF.html

<a id="ref-W01"></a>

### W01

**[W01](#ref-W01) Robert L. Lucke, “Rayleigh-Sommerfeld Fraunhofer Diffraction”, 2006, §§2–3, 식 (1), (5), (6)–(8).** 원저자 논문 전체 PDF 8쪽을 열람했다. RS obliquity를 남긴 Cartesian→direction-cosine와 역방향 Fourier 관계의 근거. 본문의 실제 ray eikonal 연결·normalization·오류상태 규약은 해당 논문이 제공한 제품 알고리즘이 아니라 명시한 자체 구현 유도다. <https://arxiv.org/pdf/physics/0604229> ; <https://doi.org/10.48550/arXiv.physics/0604229>

<a id="ref-W02"></a>

### W02

**[W02](#ref-W02) Ansys / Ken Moore, “What does the sampling correspond to in wavefront-based calculations?”** 제품의 grid centering·zero padding·OPD/PSF/OTF 관계와 다파장 scale 설명. 독립 DFT indexing을 제품과 동일하다고 주장하는 근거로 쓰지 않는다. <https://optics.ansys.com/hc/en-us/articles/42661791331859-What-does-the-sampling-correspond-to-in-wavefront-based-calculations>

<a id="ref-TH01"></a>

### TH01

**[TH01](#ref-TH01) SCHOTT, TIE-29: Refractive Index and Dispersion, August 2023.** §2.1 Abbe number, §2.2 secondary spectrum, §2.3 Sellmeier와 유효범위. [제조사 PDF](https://media.schott.com/api/public/content/aaa572afd854434fb7b3faa4bc46103f?v=c0f4fa52)

<a id="ref-TH02"></a>

### TH02

**[TH02](#ref-TH02) SCHOTT, TIE-35: Transmittance of Optical Glass, November 2020.** p.1–2 내부/외부 투과율, Beer–Lambert와 두께 환산. [제조사 PDF](https://media.schott.com/api/public/content/3d5861d824f64033aafd5353a21f7603?download=true&v=1a50ee6b)

<a id="ref-TH03"></a>

### TH03

**[TH03](#ref-TH03) Ansys, OpticStudio 2025 R1, FFT Through Focus MTF.** 일반 FFT MTF와 별도 빠른 계산 경로. [공식 문서](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_Through_Focus_MTF.html)
