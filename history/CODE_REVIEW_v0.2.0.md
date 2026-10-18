# Sequential Optics Core v0.2.0 — 후속 독립 코드 리뷰

검토일: 2026-10-18. 대상은 배포된 `sequential_optics_core_v0.2.0.zip`의 고정 사본이다. 이번 작업은 **리뷰와 재현만 수행했으며 엔진 코드를 수정하지 않았다.** 이전 R1–R9 수정 여부와 새 발견 사항을 혼동하지 않도록 새 ID를 사용한다.

## 판정

**새 구현 결함 4건과 F-θ 정의 문제 1건을 확인했다.** 특히 동일 곡면 더미 삽입으로 광선이 소실되는 문제와 detector 재표본화의 좌표 이동 불변성 위반은 데이터 생성 전에 수정할 필요가 있다. F-θ는 코드와 기존 문서가 같은 정의를 사용하지만, 일반적인 방사 방향 F-θ와 달라 이름 및 기준을 정리해야 한다.

기존 **235개 테스트를 새로 실행해 모두 통과**했다. 실행 시간은 51.515초이며 실패·오류·skip은 없다. 이 결과와 아래 반례는 양립한다. 기존 회귀 테스트가 이번 조건을 검사하지 않았기 때문이다. 반례 스크립트가 종료 코드 0으로 끝난다는 것은 결함 재현에 성공했다는 뜻이지 수정 완료를 뜻하지 않는다.

| ID | 분류 / 우선순위 | 핵심 문제 | 영향 조건 |
|---|---|---|---|
| C1 | 구현 결함 / P2, 우선 수정 | 동일 곡면의 무굴절 더미를 추가하면 유효 광선이 `escaped`로 바뀜 | 동일 위치·동일 sag의 연속 곡면, 특히 비축 |
| C2 | 수치 구현 결함 / P2, 우선 수정 | PSF와 detector를 함께 평행 이동해도 겹침 전력이 달라짐 | 작은 PSF cell을 원점에서 떨어진 detector에 resample |
| D1 | 물리량 정의 문제 / P2 | `f_theta`가 방사 방향 각도가 아닌 X/Y 개별 각도 chart를 사용 | X/Y가 함께 0이 아닌 2D 비축 field |
| C3 | 입력 경계 결함 / P2 | spectral image 합성이 잘못된 shape·음수·NaN 전력을 허용 | 외부 또는 변경된 image-result dictionary |
| C4 | 저장 계약 결함 / P2 | 허용한 NumPy 정수 설정·STOP을 JSON으로 저장하지 못함 | Python API에서 `np.int64` 사용 |

P2는 특정 조건에서 잘못된 결과 또는 정상 기능 실패를 일으키는 수정 대상이라는 의미다. 모든 처방이 실패하거나 즉시 중단해야 하는 P0/P1 결함이라고 주장하지 않는다. 입력 손상을 요구하는 C3와 정상 입력에서 발생하는 C1/C2/C4도 구별한다.

## 검토 방식과 범위

이전 구현·리뷰와 별도로 **광학 물리학자, 렌즈 설계자, 소프트웨어 개발자 역할의 새 AI 에이전트 3명**을 구성했다. 물리·소프트웨어 역할은 각각 재료/열 및 영상 재표본화 경계 검토를 보조 에이전트에 추가 분담했다. 통합 검토자가 모든 최종 재현 스크립트를 고정 소스로 다시 실행했고, 광학 물리 역할이 C1/C2의 물리적 기대값을 교차 검토했다. 사람 전문가의 인증은 아니다.

최종 통합 보고서도 세 역할이 다시 읽고 소관 수치·발생 조건·분류에 오기나 과장이 없음을 확인했다. 통합 재검사의 16개 관찰 확인은 기존 235개 unittest와 별도 기록이며 테스트 수에 합산하지 않는다.

| 관점 | 검토 대상 | 주요 독립 대조 |
|---|---|---|
| 광학 물리 | OPL/OPD, focal/afocal FFT, scalar RS, radiometry, source amplitude, n(T), CTE | 직접 autocorrelation, 원형 aperture의 해석 복소장, 두 spectral mode의 raw flux, 독립 열팽창·lensmaker 식 |
| 렌즈 설계 | 실제 교점, STOP, chief/field 정의, spot, 왜곡, T/S field curvature, 공통 focus | 무굴절 더미 삽입 불변성, 방위각 회전, 독립 Snell 미분 |
| 소프트웨어 | 입력·값 객체, 설정, 저장, 영상 합성, 보존적 resampling, 선택적 데이터 | 정상 producer 대조, 명시적 손상 입력, JSON round-trip, 좌표 이동 불변성 |

엔진의 real-ray 교점과 Snell 추적, 근축 정의량 EFL/BFL, scalar 회절 근사를 구분해 검토했다. Zemax 실행 또는 측정 데이터와의 수치 상관 검증은 이번에도 수행하지 않았다.

## C1 — 동일 곡면 더미가 가짜 광선 소실을 만든다

**위치:** `src/sequential_optics/geometry.py:276–294`, 특히 `fa==0` 및 sign-changing bracket 판정. 형상 계약은 같은 파일 153행에 coincident dummy/refracting surface를 허용한다고 명시한다.

**정상 입력:** n=1.5의 대칭 singlet, R=±50 mm, 중심 두께 3 mm, clear radius=10 mm, STOP radius=5 mm. 첫 굴절면 바로 뒤에 위치·sag가 완전히 같은 면을 추가하고 그 면 앞뒤 매질을 모두 같은 유리로 둔다. 두 처방 모두 `validate_geometry`에서 `valid`다.

이 추가 면에서는 굴절률이 변하지 않고 경로 길이가 0이다. 따라서 방향, OPL, 투과 광선 집합, spot이 보존되어야 한다.

| field (X,Y), rad | 원래 처방 | 동일 곡면 추가 후 |
|---|---|---|
| (0,0) | 217개 모두 성공 | 동일, spot 불변 |
| (0.01,0.02) | 217개 모두 성공, RMS 0.04312536724 mm | 217개 모두 `escaped`, spot 계산 불가 |
| (0,0.1) | 217개 모두 성공, RMS 0.09044853456 mm | 102개 `escaped`, 115개 성공; RMS 0.08936669979 mm, 상태 `ok` |

마지막 경우에는 quadrature 투과 비율이 1에서 **0.60828271564**로 줄고, 수치 실패를 물리적 손실처럼 포함한 spot을 반환한다. 단순 예외보다 이 조용한 결과 변화가 더 중요하다.

**원인:** 첫 교점의 양의 sag 잔차가 재현 광선에서 `5.5128124e−13 mm`다. 이는 기본 교점 tolerance `1e−10 mm`보다 작지만, 다음 동일 곡면에서는 정확한 t=0 root가 미소 음수로 이동한 것으로 보인다. 곡면 분기는 정확한 floating zero나 부호 변화만 받아들이므로 이 접속점을 놓친다. 평면 분기는 tolerance 범위의 t≈0을 허용하며 같은 더미 대조가 통과한다.

**수정 방향:** 동일 면 접속 또는 비접선 origin-on-surface를 명시적으로 처리한다. 단순히 작은 residual을 모두 root로 인정하면 기존 접선·최초 교점 문제를 되살릴 수 있으므로, 면 동일성 또는 미분/거리의 충분한 판정과 함께 처리해야 한다.

**향후 회귀 기준:** 곡면 더미 삽입 전후의 광선 방향·OPL·성공 mask·spot·투과 비율이 유지되어야 한다. 축상/비축/평면 대조와 기존 tangency/first-root 테스트도 함께 유지한다.

**증거:** `lens_design/reproduce_lens_review.py`, `lens_design/results.json`의 `LD02_coincident_curved_dummy`.

## C2 — detector overlap 면적 계산이 평행 이동에 불안정하다

**위치:** `src/sequential_optics/source_wave.py:288–290`의 `_polygon_area`, 이를 사용하는 326–339행의 보존적 재표본화.

단위 질량을 가진 0.01×0.01 mm PSF cell과 이를 완전히 포함하는 detector를 함께 이동한다. detector pitch는 0.002 mm다. 물리적 겹침은 같으므로 captured fraction은 항상 1이어야 한다.

| 같은 겹침 배치의 중심 | 기대 | 실제 |
|---|---:|---|
| (0,0) mm | 1 | 0.9999999999999999 |
| (5,5) mm | 1 | 확률 합이 1보다 크다는 `InputError` |

실제 `source_diffraction_kernel`이 생성한 결과도 검증했다. image reference=(-5,-5,20) mm, FFT cell width 약 `0.00052928734 mm`, detector 반폭 0.02 mm로 모든 FFT cell을 포함한다. 이 경우에도 오류가 발생하며, **같은 PSF·basis·상대 detector 배치를 원점으로 평행 이동한 뒤에는 0.9999999999999998**을 반환한다.

이 실제 producer의 wave 상태는 `underresolved`이며 그 사실을 증거에 보존했다. 이 결과를 검증된 광학 영상으로 주장하지 않는다. 다만 동일한 이산 cell의 overlap 면적 불변성은 phase 수렴과 독립적이며, 앞의 단일 cell 대조만으로도 수치 결함을 확인할 수 있다.

더 큰 수치 스트레스 사례인 중심 1000 mm, cell width 0.0001 mm에서는 실제 crop이 없는데 captured fraction이 **0.86147338152**가 된다. 이 13.85% 손실을 일반 카메라 조건의 오차로 확대 해석하지 않는다. 통상적인 5 mm 이동에서는 불필요한 예외까지 확인했다.

**원인:** 전역 좌표로 shoelace의 큰 두 합을 빼서 작은 다각형 면적을 구한다. 면적에 비해 위치 좌표가 크면 cancellation이 커진다. 현재 mass tolerance나 일부 분율 clamp는 이 근본 원인을 해결하지 못한다.

**수정 방향:** 각 다각형의 가까운 기준점 또는 detector-local 원점을 뺀 좌표에서 면적을 계산한다. 필요하면 안정적인 합산도 사용한다. crop 전력을 보존해야 하므로 무조건 mass를 1로 정규화하거나 허용오차만 크게 만드는 수정은 적절하지 않다.

**향후 회귀 기준:** 완전 포함·부분 crop 모두 여러 평행 이동과 cell scale에서 같은 captured mass를 유지해야 한다. 비축 projection opt-in 및 실제 producer 대조를 함께 포함한다.

**증거:** `software/boundary_probe.py`, `software/boundary_results.json`.

## D1 — 현재 F-θ는 2D 방사 방향 F-θ와 다른 양이다

**위치:** `src/sequential_optics/metrics.py:311–316`; 같은 정의를 적은 `docs/core_principles_ko.md`의 F-θ chart 설명.

현재 `Field(x,y,kind='angle')`의 방향은 `(tan x, tan y, 1)`에 비례한다. 따라서 실제 광축과의 각도와 횡방향 단위벡터는 다음과 같다.

$$
\mathbf t=(\tan x,\tan y),\qquad
\theta=\arctan\|\mathbf t\|,\qquad
\widehat{\mathbf u}=\frac{\mathbf t}{\|\mathbf t\|}.
$$

방사 방향 F-θ의 chart는 원점에서 연속적으로 0을 정의하여 다음을 사용해야 한다.

$$
\mathbf q_{F\theta}=\theta\widehat{\mathbf u}.
$$

현재 코드는 대신 `q=(x,y)`를 사용한다. 이는 X/Y 각각의 각도에 대한 chart로서는 정의할 수 있지만, 일반적인 방사 방향 F-θ와 같지 않다. 기존 문서도 현재 코드를 그대로 설명하므로 **코드-문서 불일치가 아니라 정의·명명의 문제**다. Ansys 공식 문서는 F-θ 기준을 물체 공간에서 chief ray가 이루는 각도에 연결한다.[1]

**재현:** synthetic triplet에서 실제 polar angle을 0.15 rad로 고정하고 방위각만 Y 방향에서 45° 방향으로 바꾼다. 실제 chief image 반경은 둘 다 약 `12.70884330491 mm`이며 동축 회전대칭성이 유지된다.

| 방위 | 현재 `f_theta` 왜곡 | 방사 방향 F-θ 기대값 |
|---|---:|---:|
| Y meridian | 0.68521989357% | 0.68521989357% |
| 45° diagonal | 0.30822092866% | 0.68521989357% |

차이는 약 **0.377 percentage point**다. Rectilinear 대조는 두 방향 모두 약 −0.07105439686%로 일치한다. 따라서 이 반례는 실제 chief 추적이나 회전대칭 렌즈 자체의 실패가 아니다.

**수정 방향:** Zemax의 물리적 F-θ 의미를 유지하려면 polar-angle chart와 그 역변환으로 바꾸고 원리 문서도 수정한다. 기존 componentwise chart가 필요하면 별도 이름으로 보존한다. 1D meridian에서 두 정의가 같다는 사실을 일반 2D field까지 확장하면 안 된다. 전체 Zemax 결과와의 일치까지 입증한 것은 아니다.

**증거:** `lens_design/reproduce_lens_review.py`, `lens_design/results.json`의 `LD01_f_theta_definition`.

## C3 — spectral image 합성의 물리 입력 검증이 빠져 있다

**위치:** `src/sequential_optics/imaging.py:242–250`, `registered_spectral_sum`.

grid 자체는 검사하지만 `cell_power`의 grid shape 일치, 유한성, 비음수성을 검사하지 않는다. 실제 `synthesize_image`의 정상 출력은 정상 합산된다. 그 출력 dictionary의 `cell_power`만 외부 호출자가 바꾼 경우 다음을 허용한다.

| 입력 | 실제 결과 | 필요한 처리 |
|---|---|---|
| 2×2 grid에 1×1 `cell_power=[[1]]` | cell power 합 1, irradiance×area 합 4 | broadcasting 전에 shape 오류 |
| 2×2 음수 전력 | 음수 영상 반환 | `InputError` |
| 2×2 NaN 전력 | NaN 영상 반환 | `InputError` |

이는 정상 producer가 보통 사용에서 자발적으로 음수/NaN을 만든다는 주장이 아니다. 공개 합성 API에 외부·역직렬화·변경된 결과를 넣을 때 입력 경계가 무너진다는 결함이다. 추가로 아주 큰 유한 정상 입력들의 합이 overflow할 때 출력 유한성을 검사하지 않는 점도 확인했으나, 이를 일반 광량 범위의 실패로 해석하지 않는다.

**수정 방향:** 각 `cell_power.shape==grid.shape`, finite 및 nonnegative를 검사하고, response 적용 후에도 결과 유한성을 검사한다. `wavelength_nm`과 공통 결과 계약도 검증한다. 명시적인 채널 축 이외의 NumPy broadcasting은 허용하지 않는다.

**향후 회귀 기준:** 정상 source→image→spectral 경로는 유지하고, 위 손상 입력은 도메인 `InputError`로 거부한다. 파장간 등록과 lateral color 보존도 유지한다.

**증거:** `software/reproduce_software.py`, `software/results.json`의 `SW_SPECTRAL_BOUNDARY`.

## C4 — 정상 허용 NumPy 정수가 JSON 저장에서 실패한다

**위치:** `src/sequential_optics/settings.py:100–101`, `io.py:13–18,34–35`; 허용 계약은 settings의 `Integral` 검사 및 `models.py:55`.

`SamplingSettings(pupil_rings=np.int64(2))`와 `OpticalSystem(stop_index=np.int64(0))`는 정상 생성되고 `evaluate(... first_order ...)`도 `ok`다. 그러나 `EngineSettings.to_json`과 `save_system`은 둘 다 `TypeError: Object of type int64 is not JSON serializable`로 실패한다. Python 기본 int로 같은 처방·설정을 만들면 저장·읽기 대조가 통과한다.

**원인:** 생성 단계에서 NumPy 정수를 허용하지만 저장 경로는 `asdict` 이후 표준 `json.dumps`에 NumPy scalar를 그대로 전달한다.

**수정 방향:** 값 객체 생성 시 지원 scalar를 Python 기본 scalar로 정규화하거나, 공통 JSON 변환기를 설정과 처방 저장에 일관되게 적용한다. bool/소수 index를 거부하는 기존 계약을 약화시키면 안 된다.

**향후 회귀 기준:** Python int/NumPy integer의 동일 계산 및 설정·처방 JSON round-trip을 검사한다. 계산을 끝낸 처방을 재현 가능한 형태로 저장하는 경로를 실제로 호출해야 한다.

**증거:** `software/reproduce_software.py`, `software/results.json`의 `SW_SERIALIZATION`.

## 계산 결함과 구별한 개선·명세 항목

### I1 — underresolved 상태를 최종 영상까지 보존

실제 source 파면이 grid phase step 4.2817 rad, triangle span 5.7131 rad로 `underresolved`인 사례를 만들었다. 상위 source 결과는 올바르게 그 상태를 반환하지만, `SourceKernel`과 `synthesize_image`의 최종 결과에는 상태가 없다. 전력 ledger 잔차는 0이므로 전력 보존만으로 이 품질 부족을 잡을 수 없다.

원래 source 결과를 검사하는 호출자는 이를 걸러낼 수 있고, README에도 채택 금지 지침이 있다. 따라서 이전 sampling guard가 다시 잘못된 `ok`를 반환한 문제로 세지 않았다. 다만 데이터 생성 파이프라인에는 중요하다. Kernel의 analysis status·출처·sampling 진단을 보존하고 최종 이미지까지 전달하거나, strict 모드에서 미해상 kernel을 거부하는 개선이 필요하다.

### I2 — 흡수하는 무한 물체 공간의 radiance 경계

`uniform_scene_ri`는 유한 물체에서 μℓ 감쇠를 적용하지만 `object_z=None`이면 무한 물체측 감쇠를 적용하지 않는다. μ=1/mm 사례에서 object z=−100 mm의 조도는 `1.08835e−45`, 무한 물체 설정은 0.04다.

L을 무한 물체 자체의 radiance로 정의하면 무한 경로의 극한은 0이다. 이미 입사면에 도착한 radiance로 정의하면 0.04가 성립할 수 있다. 현재 문서로는 두 경우가 구별되지 않으므로 **조건부 경계 명세 문제**로 남긴다. 투명 공기 μ=0 경로의 오류라고 세지 않는다. 적용 불가 guard 또는 명시적인 incident boundary 입력으로 정리할 수 있다.

### I3 — 값 객체의 복사·직렬화 계약

PixelGrid/SourceKernel은 직접 생성 시 bytes-backed 배열로 보호되지만 deepcopy/pickle 후 NumPy owning array로 복원되면 writable이 된다. 다만 소비 경계 `validate()` 덕분에 음수나 mass>1은 계속 거부되므로 원래 R8 전력 결함의 재발로 세지 않았다. 재구성 시 불변 배열을 재생성하는 개선이 적절하다.

FrozenDict가 포함된 settings/system의 pickle round-trip은 실패한다. ProcessPool을 공식 지원하려면 JSON 직렬화 기반 전달 또는 명시적인 pickle/reduce 계약을 정해야 한다. 현재 multiprocessing 지원을 보장한다고 문서화되어 있지 않아 개선으로 분류했다.

### I4 — 극단적인 상대 가중치의 안정성

`combine_complex_otfs`에서 finite weights `[1e308,1e308]`는 합산 overflow 때문에 정상 DC=1 대신 0을 만들 수 있다. 일반 가중치 범위의 오류가 아니라 수치 강건성 개선으로 기록한다. 양수 최대값으로 먼저 스케일한 후 정규화하면 상대비를 유지할 수 있다.

### I5 — nominal 온도 metadata

`temperature_c=None`은 재료마다 제공된 기준 n을 쓰는 명시된 convention이며 공통 실제 온도가 아니다. 현재 상위 API와 README는 이를 설명한다. 저수준 결과에도 적용 온도를 표시하면 이후 열 데이터 해석의 혼동을 줄일 수 있다. 새 계산 오류로 세지 않았다.

### I6 — 추가 검토 권고, 새 실패 미확정

강한 clipping에서는 smooth pupil용 Gauss quadrature와 별도의 경계 세분화·수렴 지표가 유용하다. 공통 초점 merit의 큰 2차 모멘트 차감은 최적 z에서 실제 centered 좌표로 다시 계산하면 더 안정적이다. 이번에는 이 두 항목에서 실용 처방의 새 오답을 확증하지 않았으므로 잠재 개선으로만 남긴다. 또한 저수준 `trace_transmission(system, ray, trace)`는 주 source 경로의 bundle provenance 검사를 거치지 않으므로, 직접 호출자의 입력 일치 책임을 명시할 필요가 있다. 정상 상위 경로에서 stale trace가 자동 생성된다는 뜻은 아니다. 상세 검토는 역할별 보고서에 있다.

## 독립 대조에서 정상 확인한 내용

| 대조 | 관찰 결과 |
|---|---:|
| 임의 비등방 complex pupil의 직접 autocorrelation과 FFT OTF | 최대 복소 오차 약 4.44×10⁻¹⁶ |
| 균일 원형 aperture의 full RS 축상 복소장 해석해 | 오차 약 7.25×10⁻¹³ |
| 두 비축 spectral mode의 정규화 전 전력 | 상대 오차 약 2.89×10⁻¹⁵ |
| n(λ,T) affine 함수의 bilinear 보간 | 최대 오차 2.22×10⁻¹⁶ |
| instantaneous CTE / reference-derivative CTE | 해석식과 부동소수점 오차 내 일치 |
| A₂–A₁₂ 전체 열팽창의 `z_T(r)=s z_0(r/s)` | 최대 sag 오차 2.78×10⁻¹⁷ mm |
| n(T)+형상 팽창을 함께 적용한 thick-lens EFL | 기대/실제 모두 50.36025862012033 mm |
| 단일 구면의 독립 Snell 미분 T/S focus, 세 방위각 | 최대 z 오차 약 5.02×10⁻⁸ mm |
| G-only 자료의 nominal 계산과 선택적 데이터 | nominal 성공, 미제공 파장·열·흡수의 해당 요청만 오류 |

이들은 선택한 모델·입력에서의 독립 불변량 또는 해석해 대조다. 모든 NA/field/비구면/온도의 정확성이나 Zemax 동일성을 보장하지 않는다. Scalar 단일 경로, sampled pupil topology, chief-normal FFT, opt-in 국소 detector 투영, 실제 코팅 위상·vector diffraction·ghost 미지원은 기존에 선언된 범위이므로 새 결함 수에 포함하지 않았다.

## 수정 권고 순서와 재검증

1. **C1, C2:** 잘못된 광선 손실과 위치 의존 전력 처리를 먼저 수정한다. 회귀 기준은 단순히 예외가 안 나는 것이 아니라 물리 불변량 보존이다.
2. **D1:** F-θ의 정의를 확정해 코드와 원리 MD/TeX를 함께 고친다. 기존 componentwise chart를 유지할 경우 이름을 분리한다.
3. **C3, C4:** 결과 저장·합성 경계의 계약을 완성한다. 정상 파이프라인과 잘못된 입력을 각각 검사한다.
4. **I1:** 최종 데이터에 품질 상태를 보존하는 채택 gate를 연결한다. `conditional`과 `underresolved`를 전력 ledger 성공과 구별한다.
5. I2–I5는 적용 환경과 병렬 실행 요구에 맞게 명세 및 강건성을 정한다.

수정 후에는 각 반례를 기대 동작 회귀 테스트로 바꾸고 기존 235개를 함께 실행해야 한다. 현재 반례 스크립트의 버그 존재 assertion을 그대로 ‘수정 검증’에 쓰면 안 된다. 이번 리뷰에서는 위 수정 작업을 수행하지 않았다.

## 재현 자료와 출처

- `physics/review.md`: 광학 물리 역할의 상세 판정과 수식 대조.
- `lens_design/review.md`: 렌즈 설계 역할의 상세 판정.
- `software/review.md`: 소프트웨어 역할의 상세 판정.
- 각 디렉터리의 `.py`, `.json`: 재현 및 실제 출력.
- `root_reproductions/`: 통합 검토자의 재실행 출력과 235개 baseline 로그.
- `run_review_probes.py`: 모든 재현을 순서대로 실행하고 핵심 관찰을 재검사하는 진입점.
- `input/engine/`: 리뷰한 변경 없는 v0.2.0 사본. 실행 시 이 src를 강제로 사용한다.
- `INPUT_IDENTITY.json`, `REVIEW_MANIFEST.json`: 입력 버전·해시 및 증거 파일의 무결성 기록.

**[1] Ansys OpticStudio 공식 도움말, Field Curvature and Distortion.** [원문](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Field_Curvature_and_Distortion.html). 2026-10-18 원문 확인. Chief-ray 기반 radial distortion, 작은 field의 reference ray, 물체 공간 chief 각도에 대한 F-θ 정의를 대조했다. 이것은 D1의 정의 근거이며 현재 엔진과 상용툴의 수치 일치 검증이 아니다.

**[2] Robert L. Lucke, Rayleigh-Sommerfeld Fraunhofer Diffraction (2006).** [원문 PDF](https://arxiv.org/pdf/physics/0604229). 2026-10-18 원문 확인. Direction-cosine Fourier 관계와 obliquity의 물리적 배경을 대조했다. 엔진의 보간·sampling guard 전체를 검증하는 출처로 확장하지 않았다.

수치 결함 C1/C2는 소스 구현과 독립 불변량으로, C3/C4는 실제 API 호출로 직접 확인했다. 확인하지 못한 문헌 수식을 새로 검증했다고 주장하지 않는다. 역할별 상세 보고서의 건수는 서로 중복될 수 있으므로 본 통합 판정의 **구현 4건 + 정의 1건**을 전체 집계로 사용한다.
