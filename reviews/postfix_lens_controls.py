"""Independent post-fix lens controls, not an engine test-suite replacement."""
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sequential_optics import Field, OpticalSystem, Surface
from sequential_optics.geometry import aim_ray, trace_ray, validate_geometry
from sequential_optics.materials import Material
from sequential_optics.source_wave import _polygon_area


def main():
    materials = {'air': Material('air', n_constant=1.),
                 'glass': Material('glass', n_constant=1.52)}
    front = Surface(0., .025, 8., conic=-.4, coefficients=(.0004, -1e-6),
                    surface_type='even_asphere', medium_after='glass')
    rear = Surface(4., -.02, 8.)
    points = np.array([[0., 0.], [.2, -.4], [-.5, .3], [.7, .5], [-.8, -.2]])
    fields = (Field(), Field(.045, -.07), Field(-.08, .025))
    max_position = max_direction = max_opl = 0.
    count = 0
    for stop in (0, 1):
        base = OpticalSystem((front, rear), materials, stop_index=stop,
                             stop_radius=3., image_z=45.)
        # Two front dummies plus one rear dummy. Equal sag despite different
        # type/zero-padding/name; STOP can be the final duplicate in the chain.
        chain = (front, replace(front, coefficients=front.coefficients+(0.,)*4),
                 replace(front, name='front observation'), rear,
                 replace(rear, surface_type='even_asphere', coefficients=(0.,)*6))
        duplicate = replace(base, surfaces=chain, stop_index=2 if stop == 0 else 4)
        assert validate_geometry(duplicate)['status'] == 'valid'
        for field in fields:
            for pupil in points:
                ra = aim_ray(base, field, pupil)
                rb = aim_ray(duplicate, field, pupil)
                a, b = trace_ray(base, ra), trace_ray(duplicate, rb)
                assert a.success and b.success, (a.status, b.status)
                assert np.max(abs(b.positions[[1, 2]]-b.positions[0])) == 0.
                assert np.max(abs(b.positions[4]-b.positions[3])) == 0.
                assert b.opl[0] == b.opl[1] == b.opl[2]
                assert b.opl[3] == b.opl[4]
                pd = float(np.max(abs(a.positions-b.positions[[0, 3]])))
                dd = float(np.max(abs(a.directions-b.directions[[0, 3]])))
                od = float(np.max(abs(a.opl-b.opl[[0, 3]])))
                assert pd < 2e-8 and dd < 2e-9 and od < 2e-8, (pd, dd, od)
                max_position = max(max_position, pd)
                max_direction = max(max_direction, dd)
                max_opl = max(max_opl, od)
                count += 1

    # Moving the aiming endpoint to a coincident chain's first event must not
    # remove earlier physical apertures; this smaller rear opening is real.
    rear_small = replace(rear, clear_radius=1., mechanical_radius=8.)
    clipped = OpticalSystem((front, rear_small, rear), materials,
                            stop_index=2, stop_radius=3., image_z=45.)
    clipped_ray = aim_ray(clipped, Field(.045, -.07), [.8, .1])
    clipped_trace = trace_ray(clipped, clipped_ray)
    assert clipped_trace.status == 'aperture_blocked'
    assert len(clipped_trace.positions) == 1

    # Independent geometric predicate: these are ordered convex rectangles,
    # with no lens trace, surface crossing, self-intersection or projection.
    areas = []
    for shift in (0., 5., 1000.):
        p = np.array([[-.005,-.005],[.005,-.005],[.005,.005],[-.005,.005]])+shift
        edges = np.roll(p, -1, axis=0)-p
        turn = edges[:,0]*np.roll(edges,-1,axis=0)[:,1]-edges[:,1]*np.roll(edges,-1,axis=0)[:,0]
        assert np.all(turn>0.)
        exact_represented = float((p[1,0]-p[0,0])*(p[3,1]-p[0,1]))
        naive = float(abs(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))/2)
        stable = float(_polygon_area(p))
        assert stable == exact_represented
        areas.append(dict(shift_mm=shift, convex=True, represented_rectangle_area=exact_represented,
                          legacy_global_shoelace_area=naive, corrected_area=stable))
    result = {'chain_ray_comparisons': count,
              'coincident_stop_chain_physical_clip_status': clipped_trace.status,
              'max_position_difference_mm': max_position,
              'max_direction_difference': max_direction,
              'max_opl_difference_mm': max_opl,
              'convex_area_controls': areas}
    print(json.dumps(result, indent=2))
    (ROOT/'reviews'/'postfix_lens_controls.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
