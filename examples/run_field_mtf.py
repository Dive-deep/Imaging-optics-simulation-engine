"""Export each field's fixed-axis T/S FFT MTF and diffraction reference.

This is a synthetic lens-only, G-default example. It does not model the LED
pixel/fill response, and preserves analysis status rather than certifying data.
"""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sequential_optics import EngineSettings, Field, RunRequest, evaluate
from sequential_optics.io import save_results
from sample_systems import synthetic_triplet


def main():
    settings = EngineSettings.from_json(ROOT/'examples/settings.json')
    theta = .03
    diagonal = float(np.arctan(np.tan(theta)/np.sqrt(2)))
    fields = (Field(),Field(0.,theta),Field(theta,0.),Field(diagonal,diagonal))
    result = evaluate(synthetic_triplet(),RunRequest(fields=fields,metrics=('mtf',),
        include_ray_history=False,fail_fast=True),settings)
    spectral = result['states'][0]['wavelengths'][0]
    output = dict(engine_version=result['engine_version'],wavelength_nm=spectral['wavelength_nm'],
        settings=result['settings'],
        convention='fixed_chief_normal_xy; T=frequency y, S=frequency x; not field-radial rotation',
        vendor_equivalence='unverified',
        scope='Synthetic lens-only focal curves; conditional/underresolved are retained',fields=[])
    for record in spectral['fields']:
        analysis = record['mtf']; wave = analysis['value']
        output['fields'].append(dict(field=record['field'],status=analysis['status'],
            analysis_status=analysis['analysis_status'],frequencies_lp_mm=wave['frequencies_lp_mm'],
            mtf_t=wave['mtf_t'],mtf_s=wave['mtf_s'],otf_t=wave['otf_t'],otf_s=wave['otf_s'],
            diffraction_mtf_t=wave['diffraction_mtf_t'],diffraction_mtf_s=wave['diffraction_mtf_s'],
            axis_convention_id=wave['axis_convention_id'],
            mtf_frequency_axes_global=wave['mtf_frequency_axes_global'],
            mtf_target_line_axes_global=wave['mtf_target_line_axes_global'],
            sampling_diagnostics=wave['sampling_diagnostics']))
    target=ROOT/'results/field_mtf_demo.json'
    save_results(output,target)
    print(target)
    print([(r['field'],r['analysis_status']) for r in output['fields']])


if __name__=='__main__':
    main()
