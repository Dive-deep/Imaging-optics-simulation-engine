# C1–C4 수정 후 독립 리뷰 — 렌즈 설계 관점

2026-10-18. 대상은 루트 `src/sequential_optics/`의 이번 수정본이다. 구현자와 별도인 렌즈 설계 역할의 AI 리뷰이며 사람 전문가 인증 또는 Zemax 수치 상관 검증은 아니다. 이 리뷰어는 엔진 소스를 수정하지 않았다. C1/C2의 광학·기하 타당성과 field별 FFT MTF 방향, D1 설명 그림을 중심으로 검토했다. C3/C4의 소프트웨어 계약은 다른 독립 리뷰에 맡긴다.

## 최종 판정

**검토 범위에서 미해소 차단 결함을 찾지 못했다.** 다만 첫 C1 수정본에서 동일 rear surface를 STOP dummy로 복제하면 aiming이 실패하는 잔여 결함을 실제로 발견했다. 구현자에게 전달한 뒤 좁은 수정이 추가되었고, 동일한 독립 반례와 광학 불변량으로 해소를 확인했다. 초기 오류를 발견하지 않았던 것처럼 검토 기록을 소급해서 바꾸지 않는다.

## C1 — 동일 곡면 접속

`geometry._coincident_successor_hit`는 작은 sag 잔차를 가진 임의 광선을 교점으로 인정하는 방식이 아니다. 직전 sequential 면을 성공적으로 만난 ray가, 동일 z/curvature/conic/A2…A12의 바로 다음 면을 만나는 경우만 이미 검증한 hit를 재사용한다. 새 면의 mechanical domain, active aperture, 법선 입사 및 Snell/TIR는 계속 적용한다. 뒤에 생략한 0 계수를 채운 뒤 처방을 비교하므로 같은 형상과 다른 직렬화 표기가 구별되지 않는다.

다음 점에서 원래 문제에 적절한 수정이다.

- 같은 점에서 추가 전파 거리 및 OPL은 정확히 0이다. 동일 매질 dummy는 방향도 보존한다.
- 독립 광선의 작은 잔차, 비동일 sag, 미소 양의 간격에는 기존 first-root solver가 계속 사용된다.
- tangent 또는 그보다 먼저 있는 미해결 root를 넘겨 뒤의 root를 고르는 완화가 없다.
- aperture 차단과 수치 실패를 서로 바꾸지 않는다.

### 후속 리뷰에서 발견하고 수정한 STOP chain 문제

독립 처방은 공기 1, 유리 1.52, 앞면 `z=0, c=.025, k=-.4, A2=.0004, A4=-1e-6`, 뒷면 `z=4, c=-.02`, clear radius 8 mm, rear STOP radius 3 mm다. 이 정상 singlet에서 `(front,rear)`를 `(front,rear,rear)`로 바꾸고 STOP을 마지막 복제 면으로 이동했다. `Field(.045,-.07)` 및 normalized pupil `(.2,-.4)`에서 원래 aiming은 성공하고, 그 ray를 복제 처방으로 추적해도 성공하지만, 복제 처방의 자체 aiming은 `first sag hit is beyond STOP target plane`으로 실패했다.

원인은 Newton의 아직 수렴하지 않은 trial ray에 대해, 이제 prefix가 된 옛 STOP 곡면을 목표 z 앞에서 반드시 만나도록 요구한 것이다. 곡면의 동일 위치 조건 하나를 중복한 탓에 정상 수렴을 먼저 차단했다.

최종 수정의 `_aiming_stop_hit`는 STOP과 정확히 같은 sag가 연속된 경우 그 체인의 **첫 면을 위치 aiming 사건으로 사용**한다. 중간 굴절은 같은 점의 위치를 바꾸지 않으므로 타당하다. 실제 추적은 여전히 모든 면의 aperture/굴절/TIR를 적용한다. 비동일 prefix의 horizon 또는 일반 first-root 규칙은 바뀌지 않았다.

독립 `postfix_lens_controls.py`는 앞면 3중·뒷면 2중 체인, 계수 0-padding과 Standard/EvenAsphere의 같은 형상, front/rear STOP, 3 fields×5 pupil×2 STOP = **30개 광선**을 비교한다. 수정 후 전부 성공했고 차이는 다음과 같다.

| 비교량 | 최대 절대 차이 |
|---|---:|
| 광선 교점 | 6.10623e-16 mm |
| 정규화 방향 성분 | 2.77556e-17 |
| OPL | 3.12250e-17 mm |

추가로 rear STOP chain의 먼저 만나는 면만 clear radius 1 mm로 줄였다. aiming 자체는 가능하지만 실제 trace는 `aperture_blocked`를 반환했다. 따라서 위치 제약을 합친 것이 실제 aperture를 없앤 것은 아니다. 동일 점의 OPL 증가는 정확히 0임도 검사했다.

STOP chain 추가 수정 후 C1 회귀 파일을 다시 실행하여 **13개 전부 통과**했다. 원래 217-ray 사례, finite object/asphere, active/mechanical aperture, net Snell, 미소 양의 간격, 다른 sag, 임의 잔차의 거부, tangent와 최초 단순 root 대조가 포함된다.

## C2 — geometry가 꼬여서 생긴 현상인가

**확인된 C2는 꼬인 렌즈 형상이나 self-intersecting polygon 때문이 아니다.** 단일 convex 사각형 PSF cell과 엄격하게 증가하는 detector edges만으로 실패할 수 있었고, 이 대조에는 렌즈 geometry 자체가 없다. 동일한 사각형과 grid를 함께 이동했을 뿐인데 겹침 면적의 합이 달라지는 floating-point cancellation이었다.

기존 global shoelace 식은 큰 좌표 곱의 차로 작은 면적을 만든다. 수정은 detector 원점을 한 번 빼고, 개별 FFT cell 중심에서도 국소 좌표로 clipping과 면적을 계산한다. 면적 공식 자체를 바꾸어 물리량을 맞추는 것이 아니라 같은 기하량을 안정적으로 평가한다. 잘린 PSF의 총 전력을 1로 다시 정규화하지 않는다.

독립 제어에서 모든 변의 회전 부호가 양수인 ordered convex rectangle을 0/5/1000 mm 이동했다. 수정 면적은 입력 float가 표현하는 실제 두 변 길이의 곱과 정확히 같았다. 예를 들어 1000 mm 위치의 면적은 `9.99999999998181e-5 mm²`인데 기존 global 식은 `9.999983012676239e-5 mm²`였다. 이미 입력 좌표를 만들 때 사라진 비트를 복구한다고 주장하지는 않는다.

별도로 C2 회귀 7개를 직접 실행해 모두 통과했다. 완전 포집, 부분 crop의 독립 1D interval-overlap 적분, 비등방·기울어진 셀의 반평면 대칭, 여러 FFT 셀의 서로 다른 전력, 정상 실제 source-wave producer 및 잘못된 grid/특이 projection 거부가 포함되어 있다. 큰 tilt/NA에서 chief-parallel detector projection이 별도의 국소 근사라는 한계는 그대로다.

## Field별 MTF와 T/S의 의미

MTF API는 요청한 각 field마다 `mtf_t`, `mtf_s`, complex `otf_t`, `otf_s`, 동일 pupil의 `diffraction_mtf_t`, `diffraction_mtf_s`를 따로 반환한다. Focal 주파수 lp/mm와 afocal cycles/rad는 별도 축이다. 신규 MTF 회귀 6개를 직접 실행해 모두 통과했다. 비축 field의 T/S가 실제로 다르고, 회전대칭 처방의 Y field를 X field로 돌리면 고정 X/Y 두 곡선이 교환되는 대조도 통과했다. 독립 직사각 pupil의 삼각형 OTF와도 비교했다.

현재 FFT profile의 방향은 `S=e1`(global/object-local X의 chief-normal 투영), `T=e2=chief×e1`이다. 이것은 field의 방사/meridional 방향을 따라 자동 회전하는 T/S가 아니다. `axis_convention_id`, `mtf_frequency_axes_global`, `mtf_target_line_axes_global`, `field_radial_axes=False`가 방향 의미를 명시한다. 주파수가 변하는 방향과 줄무늬가 뻗는 방향은 서로 직교한다.

Ansys의 공식 **Sagittal and Tangential** 설명은 OpticStudio가 field 위치와 관계없이 object-local y/x의 T/S 규약을 사용하며, Fields의 TAN이 MTF의 T/S 면을 재정의하지 않는다고 설명한다. 따라서 명명 의도를 맞추기 위해 임의 field마다 radial 축으로 조용히 교체하지 않은 것은 타당하다. 다만 이 엔진의 chief-normal 독립 scalar profile과 vendor의 실제 좌표 mapping이 수치적으로 동일함을 증명한 것은 아니다. **진정한 radial T/S와 fixed FFT T/S를 같은 열로 저장하면 안 된다.**

## D1 설명 그림

`d1_azimuth_explanation_ko.png`를 직접 확인했다. 광학계는 고정하고 방향 벡터를 z 광축 주위로 회전시키며 polar angle θ는 일정하게 둔다는 뜻이 정확하다. `(ax,ay)=(0°,10°)`와 `(7.107076°,7.107076°)`는 모두 실제 θ=10°다. 반면 두 입력 각도 성분의 norm은 후자에서 10.050923°가 된다. 이것을 F-θ의 기준 반경에 넣는 것이 D1의 원인이다.

동축 회전대칭계의 같은 θ에서 chief 영상 반경이 같은데 기준 좌표만 바뀐 사례이며, 같은 field의 tangential/sagittal 성능을 비교한 것이 아니다. 방향을 표시한 `z=1`의 좌표와 실제 lens hit를 구별한 표시도 적절하다. 설명 문서는 그림의 10°와 원래 재현의 0.15 rad를 구분한다. `metrics.py` SHA256이 v0.2.0과 동일함을 확인했으며 이번 D1은 설명만 추가했다.

## 재현과 범위

실행한 별도 script와 결과는 `postfix_lens_controls.py`, `postfix_lens_controls.json`이다. 전체 suite의 최종 집계는 root 통합 검증 기록을 따른다. 이 리뷰의 대조만으로 모든 처방의 회절 수렴, 상용툴 호환성, 실제 렌즈 성능을 인증하지 않는다.

공식 문서를 직접 읽은 근거: [Ansys OpticStudio 2025 R2 — Sagittal and Tangential](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Sagittal_and_Tangential.html). C1의 0거리 접속 및 C2의 평행이동 불변량은 이 리뷰의 독립 기하·광학 대조다.
