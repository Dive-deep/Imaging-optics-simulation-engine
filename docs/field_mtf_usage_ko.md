# Field별 T/S FFT MTF 사용

`evaluate`는 요청한 각 온도·파장·field마다 T와 S 두 곡선을 반환한다. 렌즈 단독 MTF이므로 LED 픽셀이나 fill factor를 입력할 필요가 없다. 파장을 생략하면 설정의 G를 쓴다.

```python
from sequential_optics import EngineSettings, Field, RunRequest, evaluate
from sequential_optics.io import load_system

system = load_system('examples/triplet_G_only.json')
settings = EngineSettings.from_json('examples/settings.json')
result = evaluate(system, RunRequest(
    fields=(Field(0, 0), Field(0, 0.03), Field(0.03, 0)),
    metrics=('mtf',), include_ray_history=False
), settings)

for state in result['states']:
    for spectral in state['wavelengths']:
        for record in spectral['fields']:
            analysis = record.get('mtf')
            if analysis is None or analysis['status'] != 'ok':
                continue  # 원래 오류 정보를 함께 보존
            wave = analysis['value']
            frequency = wave['frequencies_lp_mm']
            tangential = wave['mtf_t']
            sagittal = wave['mtf_s']
            diffraction_t = wave['diffraction_mtf_t']
            diffraction_s = wave['diffraction_mtf_s']
            quality = analysis['analysis_status']
            # underresolved는 정답 데이터로 채택하지 않는다.
            # conditional도 별도의 표본 수렴·모델 적용성 확인이 필요하다.
```

Finite object에서는 `object_z`를 정하고 `Field(x_mm,y_mm,kind='height')`를 사용한다. Afocal 요청은 `output_mode='afocal'`이며 주파수 key는 `frequencies_cycles_rad`다. 이 값은 lp/mm와 다른 단위다. 각주파수의 sampling guard가 부족하다고 판정한 값은 NaN으로 남으며 0으로 대체하면 안 된다. 투명 매질에서 명목 G 계산만 할 때 미사용 RGB·열·흡수 데이터는 필요 없다.

## T/S 이름과 실제 방향

Ansys FFT MTF의 고정 규약은 **T: object-X 방향으로 뻗은 줄무늬**, **S: object-Y 방향으로 뻗은 줄무늬**의 응답이다. 주파수 벡터는 줄무늬와 수직이므로 일반적인 X/Y 기준으로 T는 Y방향, S는 X방향의 주파수 응답이다. Field 방위각에 맞춰 T/S를 자동 회전하는 radial/meridional 규약과 구별해야 한다.[1,2]

이 독립 엔진은 chief-normal 계산면에서 `S=e1`, `T=e2`를 사용한다. `e1`은 global X를 이 면에 투영·정규화한 방향이고 `e2=chief×e1`이다. 주파수와 target 방향을 별도 metadata로 반환한다.

| 결과 key | 의미 |
|---|---|
| `axis_convention_id` | `fixed_chief_normal_xy` |
| `mtf_frequency_axes_global['s'/'t']` | 주파수가 증가하는 방향의 global 3D 단위벡터 |
| `mtf_target_line_axes_global['s'/'t']` | 해당 계산면에서 줄무늬가 뻗는 방향 |
| `output_basis` | e1, e2, chief 방향 순서의 직교 frame |
| `field_radial_axes` | False: field 방위각에 따라 자동 회전한 radial 축이 아님 |
| `otf_t`, `otf_s` | 위상을 보존하는 complex OTF |
| `diffraction_mtf_t`, `diffraction_mtf_s` | 같은 통과 pupil·진폭에서 위상 수차를 제거한 기준 |

동축 렌즈의 Y-meridional field에서는 익숙한 기하 meridional T/S와 맞지만, X-field나 대각 field에서는 같은 이름이 같은 방사 방향을 뜻하지 않는다. 90° 회전한 X/Y field에서 고정 축 T/S 곡선이 서로 교환되는 것은 정상이다. 기하 ray fan/field curvature의 radial T/S와도 별도의 축 계약이다.

현재 profile과 Zemax가 임의 2D field에서 pupil mapping과 수치까지 동일하다는 검증은 없다. 특히 projected X/Y는 일반적으로 서로 직교하지 않는데 본 엔진은 projected X로 직교 frame을 만든다. 따라서 `vendor_equivalence='unverified'`를 유지한다. Focal 결과는 chief-normal 가상 영상면의 image-space lp/mm이며 실제 global-z detector의 임의 tilt 전파와 같다고 해석하지 않는다.

## 실행 예제

```bash
python examples/run_field_mtf.py
```

`results/field_mtf_demo.json`에 축상·Y·X·대각 field 각각의 주파수, T/S, diffraction reference, complex OTF, 품질 상태, 방향 벡터를 저장한다. 설정과 status를 함께 저장하므로 배열만 떼어 성능 인증값으로 해석하지 않는다.

## 출처

1. [Ansys, FFT MTF](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html): fixed target orientation, 주파수 단위와 FFT 조건.
2. [Ansys, Sagittal and Tangential](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Sagittal_and_Tangential.html): 제품의 object-local X/Y 규약과 기하 meridional 정의의 구별.

2026-10-18 공식 문서를 확인했다. 정의 대조이며 상용툴 실행 결과와의 수치 상관 검증은 아니다.
