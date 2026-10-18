"""D1 angular-coordinate explanation; no optics-engine dependence."""
from pathlib import Path
import json
import math
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch
import numpy as np

OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--font', type=Path, help='Korean-capable .otf or .ttf font file')
args = parser.parse_args()
if args.font:
    font_manager.fontManager.addfont(str(args.font))
    family = font_manager.FontProperties(fname=str(args.font)).get_name()
else:
    available = {entry.name for entry in font_manager.fontManager.ttflist}
    family = next((name for name in ('Noto Sans KR', 'Noto Sans CJK KR', 'NanumGothic', 'Malgun Gothic') if name in available), None)
    if family is None:
        parser.error('Provide --font with a Korean-capable .otf or .ttf file.')
plt.rcParams.update({'font.family': family,
                     'font.size': 14, 'axes.unicode_minus': False,
                     'svg.fonttype': 'path', 'mathtext.fontset': 'dejavusans'})

theta = math.radians(10)
rho = math.tan(theta)
angle_component = math.atan(rho / math.sqrt(2))
component_deg = math.degrees(angle_component)
norm_deg = math.degrees(math.sqrt(2) * angle_component)

blue, orange, ink = '#116CB0', '#D97323', '#183247'
muted, pale, border = '#536878', '#F0F5F8', '#D5E0E7'
fig = plt.figure(figsize=(13.6, 9.4), facecolor='white')
fig.text(.052, .948, 'D1: 같은 입사각에서 방위각만 바꾼다는 뜻',
         fontsize=23, color=ink, weight='bold')
fig.text(.052, .902, '광학계는 고정하고, 입사 방향을 광축 z 주위로 회전시킨 비교입니다.',
         fontsize=15, color=muted)

# The z=1 intersection of direction vectors is a slope chart, not lens geometry.
ax = fig.add_axes([.052, .348, .435, .505])
ax.set_aspect('equal')
ax.set_xlim(-.075, .245)
ax.set_ylim(-.062, .257)
ax.axis('off')
ax.set_title('① 광축을 정면에서 본 방향 투영', loc='left', fontsize=17,
             color=ink, pad=12)
ax.add_patch(Circle((0, 0), rho, fill=False, edgecolor=border, linewidth=2,
                    linestyle='--'))
for endpoint, label, pos in [((.225, 0), 'X', (.231, -.002)),
                              ((0, .231), 'Y', (-.006, .242))]:
    ax.annotate('', endpoint, (0, 0), arrowprops={'arrowstyle': '-|>', 'color': muted,
                                                'lw': 1.5, 'mutation_scale': 15})
    ax.text(*pos, label, color=muted, fontsize=13, ha='center', va='center')
a = np.array([0.0, rho])
b = np.array([rho/math.sqrt(2), rho/math.sqrt(2)])
for p, color in [(a, blue), (b, orange)]:
    ax.annotate('', p, (0, 0), arrowprops={'arrowstyle': '-|>', 'color': color,
                                         'lw': 3.2, 'mutation_scale': 20})
    ax.plot(*p, marker='o', color=color, ms=9, zorder=4)
ax.text(-.013, rho+.013, 'A', color=blue, fontsize=20, weight='bold', ha='right')
ax.text(b[0]+.012, b[1]+.004, 'B', color=orange, fontsize=20, weight='bold')
ax.text(-.012, .132, r'$\varphi=90^\circ$', color=blue, fontsize=13, ha='right')
ax.text(.13, .097, r'$\varphi=45^\circ$', color=orange, fontsize=13)
ax.add_patch(Arc((0, 0), .112, .112, theta1=0, theta2=45,
                 color=orange, linewidth=1.8))
ax.text(.063, .018, r'$\varphi$', color=orange, fontsize=15)
# A curved arrow around z makes the meaning of "rotate about z" explicit.
ax.add_patch(FancyArrowPatch((-.003, .204), (.145, .145),
                            connectionstyle='arc3,rad=-.25',
                            arrowstyle='-|>', mutation_scale=18,
                            color=ink, linewidth=1.5))
ax.text(.100, .222, '방향만 회전', color=ink, fontsize=12, ha='center')
ax.plot(0, 0, marker='o', ms=9, color=ink, zorder=5)
ax.text(-.009, -.018, '광축 z', ha='right', va='top', color=ink, fontsize=13)
ax.text(.115, -.047, r'원의 반경 $=\tan\theta$' + '\n같은 원 위 → 실제 입사각 θ가 같음',
        ha='center', va='center', color=muted, fontsize=12)
ax.text(.01, -.14, 'z=1 평면의 방향 좌표. 실제 렌즈·광선 교점 그림은 아닙니다.',
        transform=ax.transAxes, color=muted, fontsize=10.5)

# A meridional cross-section has identical angle theta for both A and B.
side = fig.add_axes([.55, .64, .40, .20])
side.set_xlim(-.03, 1.14)
side.set_ylim(-.06, .24)
side.set_aspect('equal', adjustable='box')
side.axis('off')
side.set_title('② 실제 광축과의 각은 둘 다 10°', loc='left',
               fontsize=17, color=ink, pad=10)
side.annotate('', (1.09, 0), (0, 0),
              arrowprops={'arrowstyle': '-|>', 'color': muted, 'lw': 1.5})
side.text(1.08, -.035, 'z', color=muted, fontsize=13)
side.plot([1, 1], [0, rho], '--', lw=1.4, color=border)
side.annotate('', (1, rho), (0, 0),
              arrowprops={'arrowstyle': '-|>', 'color': ink, 'lw': 2.5,
                          'mutation_scale': 18})
side.plot([.07, .39], [.07*rho, .39*rho], color=blue, lw=5)
side.plot([.43, .75], [.43*rho, .75*rho], color=orange, lw=5)
side.add_patch(Arc((0, 0), .8, .8, theta1=0, theta2=10,
                  color=ink, linewidth=1.3))
side.text(.43, .02, r'$\theta=10^\circ$', color=ink, fontsize=14)
side.text(.68, rho+.037, 'A와 B 모두', color=ink, fontsize=12)
side.text(.01, -.08, '각 입사 방향과 광축이 만드는 평면에서 본 모습',
          transform=side.transAxes, fontsize=10.5, color=muted)

# Exact numerical comparison, displayed in degrees only for human readability.
tabax = fig.add_axes([.55, .331, .40, .226])
tabax.axis('off')
tabax.set_title('③ X/Y 입력 각도를 따로 표시하면', loc='left',
                fontsize=17, color=ink, pad=12)
cells = [
    ['방위각 φ', '90°', '45°'],
    [r'$a_x$', '0°', f'{component_deg:.6f}°'],
    [r'$a_y$', '10°', f'{component_deg:.6f}°'],
    [r'$\sqrt{a_x^2+a_y^2}$', '10°', f'{norm_deg:.6f}°'],
]
table = tabax.table(cellText=cells, colLabels=['항목', 'A', 'B'],
                    colWidths=[.39,.24,.37], loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1, 1.72)
for (r,c), cell in table.get_celld().items():
    cell.set_edgecolor(border)
    cell.set_linewidth(.8)
    cell.set_facecolor(pale if r==0 else ('#FFF2E8' if r==4 else 'white'))
    cell.set_text_props(color=ink)
    if r==0 and c==1: cell.set_text_props(color=blue, weight='bold')
    if r==0 and c==2: cell.set_text_props(color=orange, weight='bold')

rect = FancyBboxPatch((.05,.058), .90,.197, boxstyle='round,pad=.012,rounding_size=.012',
                       transform=fig.transFigure, facecolor=pale, edgecolor='none', zorder=-1)
fig.add_artist(rect)
fig.text(.07, .217, 'D1의 원인: 실제 입사각 θ 대신, 두 입력 각도의 크기를 쓴 것',
         fontsize=16, color=ink, weight='bold')
fig.text(.07,.169, r'방사 방향 F-θ 기준: $r_{\rm ref}=f\theta$' +
         '  →  A와 B가 같아야 함', fontsize=15, color=ink)
fig.text(.07,.122, r'v0.2.0 성분별 기준: $r_{\rm ref}=f\sqrt{a_x^2+a_y^2}$' +
         '  →  B에서 기준 반경이 약 0.509% 커짐', fontsize=14, color=ink)
fig.text(.07,.078, 'T/S 수차의 차이를 비교한 사례가 아닙니다. 실제 chief ray의 영상 반경은 같았습니다.',
         fontsize=13, color=muted)
fig.text(.95,.025, '도식의 각도 표시는 ° · 수식 계산은 rad', ha='right',
         fontsize=10.5, color=muted)

fig.savefig(OUT / 'd1_azimuth_explanation_ko.png', dpi=170,
            facecolor='white', bbox_inches='tight')
fig.savefig(OUT / 'd1_azimuth_explanation_ko.svg', facecolor='white', bbox_inches='tight')
plt.close(fig)
(OUT / 'angle_values.json').write_text(json.dumps({
    'illustrative_polar_angle_degrees': 10.0,
    'A': {'azimuth_degrees':90.0, 'field_x_degrees':0.0, 'field_y_degrees':10.0},
    'B': {'azimuth_degrees':45.0, 'field_x_degrees':component_deg,
          'field_y_degrees':component_deg,
          'component_angle_euclidean_norm_degrees':norm_deg},
    'reference_radius_fractional_difference': norm_deg/10 - 1,
    'note':'The diagram uses theta=10 degrees for explanation. The original engine review used theta=0.15 rad.'
}, ensure_ascii=False, indent=2), encoding='utf-8')
print(OUT / 'd1_azimuth_explanation_ko.png')
