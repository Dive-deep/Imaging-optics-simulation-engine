"""Join frozen prescriptions with their hash-matched wave results and report data."""
from pathlib import Path
import collections, json
import numpy as np
from prepare_candidates import canonical_hash, eligible

HERE=Path(__file__).resolve().parent


def span(values):
    v=np.asarray(list(values),float)
    return {'minimum':float(v.min()),'maximum':float(v.max())}


def validate_curves(frequencies,tangential,sagittal,label):
    freq=np.asarray(frequencies,float);t=np.asarray(tangential,float);s=np.asarray(sagittal,float)
    if not (freq.ndim==1 and len(freq)>1 and t.shape==s.shape==freq.shape
            and np.isfinite(freq).all() and np.isfinite(t).all() and np.isfinite(s).all()
            and freq[0]==0 and freq[-1]==20 and np.all(np.diff(freq)>0)):
        raise AssertionError(f'Invalid {label} MTF array contract')
    if not (abs(t[0]-1)<1e-8 and abs(s[0]-1)<1e-8):
        raise AssertionError(f'{label} OTF normalization failure')


def main():
    package=json.loads((HERE/'results/frozen_selection.json').read_text())
    records=package['candidates'];mono=[];fft=[]
    n=len(records)
    if n<20 or len({r['source_seed_id'] for r in records})!=n or len({r['candidate_id'] for r in records})!=n:
        raise AssertionError('At least 20 distinct qualified parents are required')
    if not all(eligible(r) for r in records):
        raise AssertionError('Every candidate must pass fresh CT>2 and CT<=20 metric gates and sampled optical constraints')
    pool=package['selection_pool']
    if n!=pool['distinct_eligible_parents']:
        raise AssertionError('The final set must preserve all distinct eligible parents')
    for r in records:
        result=json.loads((HERE/'results/final_wave'/f"{r['candidate_id']}.json").read_text())
        if canonical_hash(r['design'])!=result['design_sha256']:
            raise AssertionError('Prescription changed after wave evaluation')
        r['wave']=result['wave'];r['wave_source_sha256']=result['source_sha256']
        expected_pairs={(float(field)/12.,float(wavelength)) for field in r['input']['field_angles_deg'] for wavelength in r['input']['wavelengths_nm']}
        actual_pairs=[(float(q['field_fraction']),float(q['wavelength_nm'])) for q in r['wave']['monochromatic']]
        # Round only field-key representation; this is not a physical tolerance.
        normalized=lambda pairs:{(round(f,12),w) for f,w in pairs}
        if len(actual_pairs)!=len(expected_pairs) or normalized(actual_pairs)!=normalized(expected_pairs):
            raise AssertionError('Missing requested field/wavelength evaluation')
        for q in r['wave']['monochromatic']:
            mono.append(q)
            if 'actual_fft_diagnostic' in q:fft.append(q['actual_fft_diagnostic'])
            geometric=q.get('geometric_diagnostic')
            if not isinstance(geometric,dict):raise AssertionError('Missing separately labeled RGB geometric diagnostic')
            validate_curves(q['frequencies_lp_mm'],geometric.get('mtf_t'),geometric.get('mtf_s'),'geometric diagnostic')
            if q.get('mtf_t') is None or q.get('mtf_s') is None:
                if q['status']!='not_applicable' or not q.get('geometric_diagnostic'):
                    raise AssertionError('Missing both primary wave and separately labeled geometric diagnostic')
                continue
            validate_curves(q['frequencies_lp_mm'],q['mtf_t'],q['mtf_s'],'scalar wave')
        g=r['geometric_mtf_proxy']
        r['wave']['geometric_polychromatic_diagnostic']=[
            {'field_fraction':float(field),'frequencies_lp_mm':g['frequencies_lp_mm'],
             'mtf_t':g['MTF_t'][i],'mtf_s':g['MTF_s'][i],
             'method':'geometric ray-intercept OTF, no diffraction; separate from conditional scalar wave',
             'spectral_weighting':'equal incident RGB; common centroid, preserved lateral color',
             'status':'geometric_diagnostic_only'} for i,field in enumerate(g['fields'])]
        for q in r['wave']['geometric_polychromatic_diagnostic']:
            validate_curves(q['frequencies_lp_mm'],q['mtf_t'],q['mtf_s'],'RGB geometric diagnostic')
        r['qualification']['wave_status_counts']=dict(collections.Counter(q['status'] for q in r['wave']['monochromatic']))
        r['qualification']['surrogate_training_label_ready']=False
    metrics=[r['metrics'] for r in records]
    ray_counts_4fields={r['candidate_id']:sum(r['metrics']['all_requested_native_status_counts'].values()) for r in records}
    ray_counts_13fields={r['candidate_id']:r['aperture_fit']['sampled_ray_count'] for r in records}
    summary={'selected_count':len(records),'distinct_source_seeds':len({r['source_seed_id'] for r in records}),
       'material_order_count':len({tuple(r['design']['materials']) for r in records}),
       'aperture_modes':dict(collections.Counter(r['design']['aperture_mode'] for r in records)),
       'glass_position_counts':dict(collections.Counter(str(next(i+1 for i,n in enumerate(r['design']['materials']) if n not in ('PMMA','PC'))) for r in records)),
       'EPD_mm':span(m['aperture']['EPD_mm'] for m in metrics),
       'NA_588_onaxis':span(m['NA_rgb_onaxis_conservative_complete_rim'][1] for m in metrics),
       'NA_allRGB_onaxis':span(v for m in metrics for v in m['NA_rgb_onaxis_conservative_complete_rim']),
       'vertex_BFL_mm':span(m['vertex_BFL_mm'] for m in metrics),
       'paraxial_BFL_588_mm':span(m['first_order_primary588nm']['paraxial_bfl_mm'] for m in metrics),
       'EFL_588_mm':span(m['first_order_primary588nm']['efl_mm'] for m in metrics),
       'center_thickness_mm':span(r['design']['axial_distances_mm'][j] for r in records for j in (0,2,4)),
       'air_gap_mm':span(r['design']['axial_distances_mm'][j] for r in records for j in (1,3)),
       'stop_distance_mm':span(r['design']['stop_distance_mm'] for r in records),
       'TTL_first_vertex_to_image_mm':span(sum(r['design']['axial_distances_mm']) for r in records),
       'fullfield_chief_height_mm_by_wavelength':{w:span(m['by_wavelength_nm'][w]['fields']['1.0']['chief_height_mm'] for m in metrics) for w in ('486.0','588.0','656.0')},
       'fullfield_RI_uncoated_RGB':span(m['polychromatic_fields']['1.0']['RI_uncoated_fresnel'] for m in metrics),
       'fullfield_RI_uncoated_588':span(m['by_wavelength_nm']['588.0']['fields']['1.0']['RI_uncoated_fresnel'] for m in metrics),
       'fullfield_distortion_588_percent':span(m['by_wavelength_nm']['588.0']['fields']['1.0']['distortion_percent_paraxial_EFL'] for m in metrics),
       'axis_spot_RGB_RMS_radius_mm':span(m['polychromatic_fields']['0.0']['spot_geometric_common_RGB_centroid']['rms_radius_mm'] for m in metrics),
       'fullfield_spot_RGB_RMS_radius_mm':span(m['polychromatic_fields']['1.0']['spot_geometric_common_RGB_centroid']['rms_radius_mm'] for m in metrics),
       'primary_wave_status_counts':dict(collections.Counter(q['status'] for q in mono)),
       'actual_fft_status_counts':dict(collections.Counter(q['status'] for q in fft)),
       'strict_pointwise_quadrature_pass_count':sum(q.get('quadrature_convergence_strict_relative_005',False) for q in mono),
       'primary_monochromatic_requested_pairs':len(mono),
       'primary_monochromatic_computed_curve_pairs':sum(q.get('mtf_t') is not None for q in mono),
       'polychromatic_computed_curve_pairs':sum(len(r['wave']['polychromatic']) for r in records),
       'geometric_monochromatic_curve_pairs':sum(bool(q.get('geometric_diagnostic')) for q in mono),
       'geometric_polychromatic_curve_pairs':sum(len(r['wave']['geometric_polychromatic_diagnostic']) for r in records),
       'actual_fft_requested_diagnostic_pairs':sum(sum(abs(q['wavelength_nm']-588.)<1e-9 and any(abs(q['field_fraction']-f)<1e-9 for f in (0.,.8)) for q in r['wave']['monochromatic']) for r in records),
       'actual_fft_recorded_diagnostic_pairs':len(fft),
       'sampled_ray_count_4fields_RGB_metrics_per_candidate':next(iter(ray_counts_4fields.values())) if len(set(ray_counts_4fields.values()))==1 else None,
       'sampled_ray_count_4fields_RGB_metrics_by_candidate':ray_counts_4fields,
       'sampled_ray_count_13fields_RGB_per_candidate':next(iter(ray_counts_13fields.values())) if len(set(ray_counts_13fields.values()))==1 else None,
       'sampled_ray_count_13fields_RGB_aperture_fit_by_candidate':ray_counts_13fields,
       'each_lens_CT_strictly_over2_and_at_most20_mm':all(eligible(r) for r in records),
       'all_distinct_eligible_parents_retained':n==pool['distinct_eligible_parents'],
       'presentation_order':'glass position, material order, source seed ID; no overall rank',
       'wave_seconds_sum_not_wallclock_benchmark':sum(r['wave']['seconds'] for r in records),
       'validation_scope':'sampled optical feasibility and explicitly conditional scalar curves; no manufacturing/global-optimum/Zemax/cropped-sensor certification'}
    package['summary']=summary
    (HERE/'candidates.json').write_text(json.dumps(package,indent=2,ensure_ascii=False,allow_nan=False))
    (HERE/'results/final_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False))
    aperture_description='、'.join(f"{key}: {value}개" for key,value in summary['aperture_modes'].items())
    lines=['## 10. 최종 선정·검증 요약','',
      f"새 CT 제한으로 다시 계산·검증한 합격 checkpoint {pool['eligible_checkpoints']}개를 독립 parent별로 묶었다. **서로 다른 seed의 합격 후보 {n}개를 전부 보존**했으며 상위 일부나 단일 winner를 고르지 않았다. 최종 재료 순서는 {summary['material_order_count']}종, 유리 위치별 수는 `{summary['glass_position_counts']}`다. 각 렌즈는 **2 mm < CT ≤ 20 mm**를 실제 처방과 새 수치 검사에서 모두 만족한다.",'',
      '기하 screening 진단 점수는 `0.6×평균 기하 MTF 정규화 AUC + 0.4×최악 field/axis AUC − 0.25×|RGB 무코팅 RI(1F)−1|`이다. 같은 parent의 여러 checkpoint 중 하나를 고를 때만 사용하며, 서로 다른 후보들의 전체 순위나 파동 MTF 우열을 뜻하지 않는다. 표와 그림은 유리 위치·재료 순서·원본 seed ID 순서로 배치한다. `CT20_XXXX` ID는 원본 seed에서 유도하여 후보 추가나 진단 점수 변화로 바뀌지 않는다. `provenance.generation_history`는 보존된 이전 출발점·최적화 이력과 새 `ct20_continuation`/`ct20_polish`/`ct20_rescue`/`ct20_stop_scan`, `parent_input.source_result` 연결을 담는다. STOP 스캔 이력에는 원래 반경과 변경량을 보존하며, 최종 변경 후에도 전체 광선·NA·CT 조건을 다시 검사한다.','',
      f"구경 모드 집계는 **{aperture_description}**다. `front_diameter50`은 첫 렌즈 전면 clear diameter 50 mm 대안이며 EPD 50 mm 설계라고 부르지 않는다. 축상 NA≥0.7은 세 파장 각각 검사했지만 주변 시야의 모든 방향 cone NA≥0.7을 강제한 것은 아니다. 상고6.4 mm는588 nm chief 기준이다. 다른 파장의 lateral color와 spot이 센서 가장자리를 넘어가는 양은 데이터에 남기며, detector crop·pixel response를 적용하지 않았다.",'',
      '| 집계 항목 | 최소 | 최대 |','|---|---:|---:|']
    labels={'EPD_mm':'실제 EPD [mm]','NA_allRGB_onaxis':'RGB 축상 complete-rim NA','vertex_BFL_mm':'실제 상면 BFL [mm]',
            'paraxial_BFL_588_mm':'근축 BFL,588nm [mm]','EFL_588_mm':'EFL,588nm [mm]','center_thickness_mm':'중심 두께 [mm]',
            'air_gap_mm':'축상 공기 간격 [mm]','stop_distance_mm':'앞 STOP 거리 [mm]','TTL_first_vertex_to_image_mm':'첫 꼭짓점–상면 TTL [mm]',
            'fullfield_RI_uncoated_RGB':'1F 무코팅 RGB RI','fullfield_distortion_588_percent':'1F 왜곡,588nm [%]',
            'axis_spot_RGB_RMS_radius_mm':'축상 RGB spot RMS radius [mm]','fullfield_spot_RGB_RMS_radius_mm':'1F RGB spot RMS radius [mm]'}
    for key,label in labels.items():v=summary[key];lines.append(f"| {label} | {v['minimum']:.7g} | {v['maximum']:.7g} |")
    for w,v in summary['fullfield_chief_height_mm_by_wavelength'].items():lines.append(f"| 1F chief 상고,{float(w):g}nm [mm] | {v['minimum']:.9g} | {v['maximum']:.9g} |")
    lines+=['',f"파동 평가는 전체 {n}개 후보에서 {len(mono)}개의 field/파장별 T/S 쌍을 요청했고, 주 파동 곡선을 계산할 수 있었던 것은 {summary['primary_monochromatic_computed_curve_pairs']}쌍, RGB 합성은 {summary['polychromatic_computed_curve_pairs']}field다. 적용 불가 결과는 단일 boundary ray-map의 fold/caustic 등 상세 사유와 함께 보존했다. 128→256 적분 격자 점별0.5% 기준 통과는 {summary['strict_pointwise_quadrature_pass_count']}/{len(mono)}이다. 상태 집계: `{summary['primary_wave_status_counts']}`. 1024 FFT 진단은588nm의0/.8F 두field씩 {summary['actual_fft_requested_diagnostic_pairs']}개를 계획했으나 primary 적용 조건 실패 시 실행을 생략했고, 실제 저장된 진단은 {len(fft)}개, 상태는 `{summary['actual_fft_status_counts']}`다.",'',
      f"모든 {n}개 후보의 field/RGB에 대해 기하 ray-intercept OTF/MTF 진단 곡선도 따로 제공한다. 기하 단색 T/S 쌍은 {summary['geometric_monochromatic_curve_pairs']}개, RGB 합성 T/S 쌍은 {summary['geometric_polychromatic_curve_pairs']}개다. 파동 모델 적용 불가 패널은 이를 명시하고 기하 곡선을 표시한다. 기하 MTF에는 diffraction이 없으므로 FFT나 Huygens MTF의 계산 성공으로 바꾸어 표시하지 않는다. 실패한 파동 값을 기하 값으로 덮어쓰지 않는다.",'',
      '수렴 통과도 scalar 모델·균일 angular pupil amplitude 가정 안의 적분 비교다. 최종 처방 전체의 ray-mesh·vector 고NA 물리 정확도 검증으로 확대하지 않는다. 이전 실험의 격자 probe는 이력 자료이며 CT20 최종처방의 입력hash가 다르면 새 후보 인증으로 재사용하지 않는다.','',
      '| 후보 | 시작 seed / 방법 | 재료 순서 | CT 3매 mm | EPD mm | NA588 | 실제 BFL mm | RGB RI1 | 1F RGB RMS radius μm |','|---|---|---|---|---:|---:|---:|---:|---:|']
    for r in records:
        m=r['metrics'];method='특허 변환' if int(r['source_seed_id'].split('_')[1])>=97 else '수식·ABCD'
        history=r['provenance']['generation_history']
        if any('focus_similarity_scan' in h for h in history):method+=' + 공통초점/scale'
        elif any('refinement' in h for h in history):method+=' + OTF 재최적화'
        if any('ct20_continuation' in h for h in history):method+=' + CT20 연속 재최적화'
        if any('ct20_polish' in h for h in history):method+=' + CT20 polish'
        if any('ct20_rescue' in h for h in history):method+=' + CT20 rescue'
        if any('ct20_stop_scan' in h for h in history):method+=' + STOP 반경 스캔'
        ct=' / '.join(f"{r['design']['axial_distances_mm'][j]:.4f}" for j in (0,2,4))
        lines.append(f"| {r['candidate_id']} | {r['source_seed_id']} / {method} | {' / '.join(r['design']['materials'])} | {ct} | {m['aperture']['EPD_mm']:.3f} | {m['NA_rgb_onaxis_conservative_complete_rim'][1]:.5f} | {m['vertex_BFL_mm']:.4f} | {m['polychromatic_fields']['1.0']['RI_uncoated_fresnel']:.4f} | {1000*m['polychromatic_fields']['1.0']['spot_geometric_common_RGB_centroid']['rms_radius_mm']:.2f} |")
    lines+=['','이 표와 전체 단면 그림은 합격한 독립 후보군 전체를 같은 수준으로 제시한다. 유리 위치·재료·seed ID 배치는 품질 순위가 아니다. 좋은 축상 단색 성능을 가진 후보도 주변 시야·RGB 합성은 크게 낮아질 수 있으므로 전 시야 고MTF 완성 설계로 해석하지 않는다.','']
    method=HERE/'METHOD_KO.md';before=method.read_text().split('## 10. 최종 선정·검증 요약')[0]
    method.write_text(before+'\n'.join(lines))
    print(json.dumps(summary,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
