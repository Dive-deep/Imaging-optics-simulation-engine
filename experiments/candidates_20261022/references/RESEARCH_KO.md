# 공개 처방과 계산 정의 조사

조회일: 2026-10-22. 이 문서는 공개 자료의 조사·선별 근거다. 최종 후보의 실제 성능은 새 광선 추적 결과로 판단하며, 특허가 이번 사용자 조건을 이미 만족한다고 해석하지 않는다.

## 직접 초기값으로 사용할 수 있는 처방

[US3944337A, Three element objective lens with ellipsoidal, hyperbolic, and parabolic surfaces](https://patents.google.com/patent/US3944337A/en)는 앞쪽 두 플라스틱의 conic 비구면과 뒤쪽 구면 유리로 구성된 공기 간격 3매 처방을 공개한다. 그림을 직접 확인해 조리개가 첫 렌즈 앞에 있음을 확인했다. 수치 Example I–IV를 `US3944337A_examples.json`에 전사했다. 문헌상의 f/1.4–1.57은 이번 목표보다 느리므로 구조 출발점일 뿐이다.

원래 두 플라스틱은 PMMA와 polystyrene이다. 후자는 사용자가 허용한 PC와 같지 않으며, 유리도 사용자 FDS18W/TAFD45와 다르다. 또한 원자료의 N_D 기준은 sodium D 589.3nm이고 사용자 nd는 Fraunhofer d선이다. 따라서 문헌 숫자를 사용자 재료의 실측값으로 바꾸어 부르지 않는다. 현재 재료 모델로 바꾼 뒤 재최적화한다.

이 특허의 conic K는 명시적으로 `K=-eccentricity²`다. 현재 정의와 같은 부호를 사용한다. 원래 고차 비구면 계수는 모두0이다. Example IV의 세 번째 면은 구면이며 나머지 세 플라스틱 면은 conic이다. 12차 초과 항을 잘라야 하는 사례는 없다.

변환 절차는 다음과 같다.

1. 12°에서6.4mm인 무왜곡 기준으로 `f_target=6.4/tan(12°)=30.1096327mm`를 계산한다.
2. `s=f_target/f_patent`로 모든 반지름·축방향 길이를 곱하고, 곡률은s로 나눈다. Conic은 그대로다. 일반적인 A_m 계수는 `A_m,new=A_m,old/s^(m-1)`이나, 해당 원처방에서는 모두0이다.
3. 허용된 PMMA/PC/유리 조합으로 교체한다. 두께>2mm·공기 간격≥1mm 제약에 맞추어 수정한다. 이 순간부터 원특허 처방이 아닌 수정된 초기값이다.
4. 먼저 작은 구경에서 광선이 통과하고 초점이 형성되는지 확인한 뒤 목표 구경까지 단계적으로 키운다. 단순 확대는 f/#를 보존하므로 EPD50mm를 달성하지 않는다.
5. 첫 렌즈 opening50mm 해석과 EPD50mm 해석을 혼합하지 않고, 각 후보에 실제 stop/입사동공/렌즈 clear diameter를 기록한다.

## 검토했지만 직접 전사하지 않는 자료

[US7072123B2, Imaging lens including three lens components](https://patents.google.com/patent/US7072123B2/en)는2플라스틱+1유리와 앞조리개/내부조리개의 설계 예를 제공한다. 그러나 f/2.81–2.82이며, 일부 면에 A2 또는 홀수차 항이 있다. A4–A12만 허용한 이번 모델에 계수를 그대로 복사할 수 없으므로 구조 참고로만 사용한다. 이 자료의 K 표기 역시 새 코드에 넣기 전에 sag 식을 확인해야 한다.

[US4892398A, Triplet type objective with a diaphragm arranged behind the lens system](https://patents.google.com/patent/US4892398A/en)는2플라스틱과 높은 굴절률의 구면 유리를 사용하는 +/−/+ 계열을 뒷받침한다. 뒤쪽 조리개와 f/3.5 예가 이번 우선 조건과 멀어서 직접 초기값으로 선택하지 않았다.

본 조사에서 사용자의 모든 조건을 동시에 만족하는 공개 완성 처방은 확인하지 못했다. '검색되지 않음'을 물리적 불가능 증명으로 바꾸지 않는다. 수학적 초기값과 실제 최적화가 필요한 이유다.

## Surface 정의

[Ansys OpticStudio Even Asphere 공식 설명](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Even_Asphere.html)은 기저 conic과 짝수차 radial polynomial의 합을 사용한다. 이번 허용 형태는

`z(r)=c r²/[1+sqrt(1-(1+k)c²r²)]+A4r⁴+A6r⁶+A8r⁸+A10r¹⁰+A12r¹²`.

단위가mm일 때 A_m의 단위는mm^(1−m)이다. 이번 glass의 standard 제한은 사용자 지시에 따라 `k=0`이고 고차항도0인 구면/평면으로 제한한다. 소프트웨어의 일반 Standard surface가 conic도 지원한다는 사실로 사용자 제한을 확대하지 않는다.

## RI와 MTF의 비교 기준

[Ansys Relative Illumination 공식 설명](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v25101/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Relative_Illumination.html)은 상점에서 본 exit pupil의 유효 projected solid angle 적분으로 상대조도를 구한다. 단순 cos⁴는 얇고 느리고 수차가 없는 제한된 렌즈의 근사다. 이번 고NA 후보의 RI를 cos⁴12°라는 공통 상수로 표기하면 후보 비교가 되지 않는다. 클리핑·pupil mapping·실광선 transmission을 반영하고, 미지정 코팅은 별도 가정으로 남긴다.

[Johnson & Sasian (2020), Image distortion, pupil coma, and relative illumination](https://doi.org/10.1364/AO.392303)의 공개 초록은 distortion·pupil coma·pupil magnification이 RI에 관여함을 설명한다. 본 조사에서는 유료 본문을 읽지 않았으며 초록을 넘어서는 식을 이 논문에서 검증했다고 주장하지 않는다.

[Ansys FFT와 Huygens PSF 비교](https://optics.ansys.com/hc/en-us/articles/42661795535251-What-is-the-difference-between-the-FFT-and-Huygens-PSF)는 일반 FFT PSF에 대해 f/1.5 이상, 작은 chief 입사각 및 작은 pupil aberration 조건을 제시한다. 이번 f/~0.6·NA≥0.7 목표는 기존 FFT MTF를 무조건 인증 지표로 쓰기 어려운 영역이다. Huygens도 구현과 표본 수에 대한 검증이 필요하다.

[Ansys FFT MTF 공식 설명](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v26102/en/OpticStudio_User_Guide/OpticStudio_Help/topics/FFT_MTF.html)은 pupil phase가 급변하면 aliasing이 생기고, cosine space의 ray 분포가 불균일하면 FFT가 부정확할 수 있음을 설명한다. 큰 수차의 초기 설계에서는 낮은 공간주파수 geometric MTF를 참고할 수 있으나, 이를 diffraction FFT MTF라고 바꾸어 표시하지 않는다.

[Ansys MTF 최적화 설명](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Optimizing_for_MTF.html)은 MTF가 이미 첫 zero를 지난 상태에서는 직접 MTF local optimization이 잘 작동하지 않을 수 있어 먼저 RMS spot/wavefront를 줄이는 접근을 제안한다. 따라서 초기 수치 탐색은 다파장·다field spot/초점/개구/상고·기하 제약을 같이 줄이고, 최종 후보의0–20lp/mm 곡선을 별도로 계산한다.

## 보존한 근거

- `US3944337A_examples.json`: 전사한 네 공개 초기 처방, 원래 재료/조건 및 제한.
- `US3944337_drawing.png`: 앞조리개 배치를 직접 확인한 도면.
- `US3944337A.html/.txt`, `US7072123B2.html/.txt`, `US4892398A.html/.txt`: 조회한 특허 원자료와 검색용 텍스트.
- `ansys_relative_illumination.html`, `ansys_even_asphere.html`: 공식 정의 조회 스냅샷.
- `primary_sources.json`, `sources_sha256.json`: 출처 목록과 로컬 보존파일 해시.
