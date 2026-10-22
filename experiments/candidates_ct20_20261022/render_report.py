"""Physical ray/metric reports for initial triplet candidates.

Reads already calculated metrics and wave curves; never fabricates absent MTF.
Layouts and spot panels are freshly traced from the same stored prescription.
The offline HTML embeds all images/data and requires no external service.
"""
from pathlib import Path
import argparse, base64, html, io, json, math, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib import font_manager

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
from model import native_surfaces, sag, trace_design, first_order, geometry_checks, export_prescription
from prepare_candidates import eligible

FONT=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
if FONT.exists():
    font_manager.fontManager.addfont(str(FONT))
    plt.rcParams['font.family']=font_manager.FontProperties(fname=str(FONT)).get_name()
plt.rcParams.update({'axes.unicode_minus':False,'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white',
                     'pdf.fonttype':42,'svg.fonttype':'none'})
FIELD_COLORS={0.:'#1c6e8c',.5:'#399b60',.8:'#dc8c26',1.:'#b33a50'}
WAVE_COLORS={486.:'#3366aa',588.:'#2e8b57',656.:'#c94545'}
MATERIAL_COLORS={'PMMA':'#b9d9ef','PC':'#c3e3cf','FDS18W':'#dec6ed','TAFD45':'#efd4aa'}
FIELDS=(0.,.5,.8,1.)


def scalar(x, default=None):
    try:
        v=float(x)
        return v if math.isfinite(v) else default
    except (TypeError,ValueError):return default


def fmt(x, digits=3, suffix=''):
    x=scalar(x)
    return '—' if x is None else f'{x:.{digits}f}{suffix}'


def record_id(record):return str(record.get('candidate_id') or record['design'].get('candidate_id') or record['design'].get('seed_id','candidate'))

def normalize_wave(wave):
    if not isinstance(wave,dict):return {'monochromatic':[],'polychromatic':[]}
    if 'monochromatic' in wave:return wave
    if 'field_fraction' in wave:return {'monochromatic':[wave],'polychromatic':[]}
    return {'monochromatic':[],'polychromatic':[]}


def numeric_leaf(obj,names):
    """Read alternative schema keys, never infer a missing optical value."""
    if not isinstance(obj,dict):return None
    for name in names:
        if name in obj and scalar(obj[name]) is not None:return scalar(obj[name])
    return None


def metric_rows(metrics):
    if not isinstance(metrics,dict):return []
    if isinstance(metrics.get('by_wavelength_nm'),dict):
        rows=[]
        for wavelength,group in metrics['by_wavelength_nm'].items():
            for field,raw in group.get('fields',{}).items():
                row=dict(raw);row['wavelength_nm']=float(wavelength);row['field_fraction']=float(field)
                row['rms_radius_mm']=raw.get('spot_geometric',{}).get('rms_radius_mm')
                row['relative_illumination_uncoated']=raw.get('RI_uncoated_fresnel')
                row['relative_illumination_geometric']=raw.get('RI_geometric')
                row['distortion_percent']=raw.get('distortion_percent_paraxial_EFL')
                row['na_rim_min']=raw.get('image_space_NA',{}).get('conservative_complete_rim_NA')
                rows.append(row)
        return rows
    for key in ('per_wavelength_field','per_rgb_field','per_field_wavelength','records','per_wave','monochromatic','field_metrics'):
        value=metrics.get(key)
        if isinstance(value,list):return value
    return []


def get_metric_row(metrics,field,wavelength=588.):
    for row in metric_rows(metrics):
        f=numeric_leaf(row,('field_fraction','field'))
        w=numeric_leaf(row,('wavelength_nm','wavelength'))
        if f is not None and abs(f-field)<1e-7 and (w is None or abs(w-wavelength)<1e-7):return row
    return {}


def candidate_summary(record):
    d=record['design'];m=record.get('metrics',{});w=normalize_wave(record.get('wave',{}))
    first=first_order(d);full=get_metric_row(m,1.);axis=get_metric_row(m,0.)
    qualification=record.get('qualification') or record.get('status') or m.get('qualification') or m.get('status') or 'initial_candidate'
    if isinstance(qualification,dict):
        qualification=qualification.get('status') or ('기하 수치조건 통과 · 파동/제조 미인증' if qualification.get('geometric_sampled_constraints')=='passed' else '초기 후보 · 상세 판정 참조')
    na=numeric_leaf(axis,('na_rim_min','NA_rim_min','image_na_rim_min','na_image','NA','na','image_space_na'))
    ri=numeric_leaf(full,('relative_illumination_uncoated','RI_uncoated','ri_uncoated','relative_illumination_fresnel'))
    dist=numeric_leaf(full,('distortion_percent','rectilinear_distortion_percent','distortion_chief_slope_percent'))
    rms=numeric_leaf(full,('rms_radius_mm','spot_rms_radius_mm'))
    mtf20=None
    pp=next((r for r in w.get('polychromatic',[]) if abs(float(r['field_fraction'])-1.)<1e-7),None)
    if pp and pp.get('mtf_t') is not None and pp.get('mtf_s') is not None:mtf20=min(float(pp['mtf_t'][-1]),float(pp['mtf_s'][-1]))
    return {'candidate_id':record_id(record),'materials':' / '.join(d['materials']),
      'source_seed_id':record.get('source_seed_id',d.get('seed_id')),
      'center_thickness_mm':[float(d['axial_distances_mm'][j]) for j in (0,2,4)],
      'center_thickness_requirement':'2 mm < each CT <= 20 mm',
      'aperture_mode':d['aperture_mode'],'EPD_mm':m.get('aperture',{}).get('EPD_mm',2*d['stop_radius_mm']), 'nominal_STOP_diameter_mm':2*d['stop_radius_mm'],'first_clear_diameter_mm':2*d['clear_radii_mm'][0],
      'EFL_mm':first['efl_mm'],'BFL_vertex_mm':d['axial_distances_mm'][5],'NA_588':na,
      'RI_edge_588':ri,'distortion_edge_588_percent':dist,'spot_edge_588_rms_mm':rms,
      'poly_MTF20_worst_TS_edge':mtf20,'qualification':str(qualification)}


def plot_layout(ax,design,*,compact=False):
    s=native_surfaces(design)
    for lens in range(3):
        front,back=s[2*lens:2*lens+2]
        yf=np.linspace(-front[8],front[8],401);yb=np.linspace(-back[8],back[8],401)
        zf=front[0]+sag(front[None,:],abs(yf))[0][0]
        zb=back[0]+sag(back[None,:],abs(yb))[0][0]
        color=MATERIAL_COLORS.get(design['materials'][lens],'#d8dce1')
        ax.fill(np.r_[zf,zb[::-1]],np.r_[yf,yb[::-1]],color=color,alpha=.75,zorder=1,linewidth=0)
        ax.plot(zf,yf,color='#3a4857',lw=.8,zorder=3);ax.plot(zb,yb,color='#3a4857',lw=.8,zorder=3)
        if not compact:
            ax.text((front[0]+back[0])/2, max(front[8],back[8])+2.,design['materials'][lens],ha='center',va='bottom',fontsize=8,color='#334155')
    zstop=-float(design['stop_distance_mm']);radius=float(design['stop_radius_mm']);zimage=float(sum(design['axial_distances_mm']))
    outer=max(float(s[:,8].max()),radius)+5
    for sign in (-1,1):ax.plot([zstop,zstop],[sign*radius,sign*outer],color='#202936',lw=2.1,zorder=5)
    ax.plot([zimage,zimage],[-6.4,6.4],color='#202936',lw=2,zorder=5)
    ax.axhline(0,color='#9aa5b1',lw=.55,ls=':',zorder=0)
    ray_positions=np.array([-1.,-.67,-.33,0.,.33,.67,1.]);xy=np.c_[np.zeros(len(ray_positions)),ray_positions]
    traced=trace_design(design,fields=FIELDS,wavelengths=(588.,),xy=xy)
    for fi,field in enumerate(FIELDS):
        color=FIELD_COLORS[field]
        theta=np.radians(12*field)
        for ri,fraction in enumerate(ray_positions):
            origin=np.array([0.,fraction*radius,zstop]);source=origin.copy();source[2]-=8.;source[1]-=8*np.tan(theta)
            hist=traced['history'][0,fi,ri];hist=hist[np.all(np.isfinite(hist),axis=1)]
            points=np.vstack([source,origin,hist])
            valid=traced['status'][0,fi,ri]==0 and bool(traced['image_forward'][0,fi,ri])
            if valid:points=np.vstack([points,[*traced['image_xy'][0,fi,ri],zimage]])
            chief=fraction==0;rim=abs(fraction)==1
            ax.plot(points[:,2],points[:,1],color=color,lw=1.0 if chief else (.75 if rim else .4),
                    ls='--' if rim else '-',alpha=.95 if chief or rim else .35,zorder=4)
            if not valid:ax.plot(points[-1,2],points[-1,1],'x',color=color,ms=3,zorder=7)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xlim(zstop-9,zimage+3);ax.set_ylim(-outer-2,outer+6 if not compact else outer+1)
    ax.set_xlabel('광축 z [mm]' if not compact else 'z [mm]');ax.set_ylabel('y [mm]')
    if not compact:
        ax.text(zstop,outer+1,'STOP',ha='center',va='bottom',fontsize=7)
        ax.text(zimage,9.,'IMA\n±6.4 mm',ha='center',va='bottom',fontsize=7)
        handles=[Line2D([0],[0],color=FIELD_COLORS[f],lw=1.3,label=f'{f:g}F ({12*f:g}°)') for f in FIELDS]
        handles.extend([Line2D([0],[0],color='#555',lw=1.,label='chief'),Line2D([0],[0],color='#555',ls='--',lw=1.,label='±marginal')])
        ax.legend(handles=handles,loc='lower left',ncol=6,frameon=False,fontsize=7,bbox_to_anchor=(0,-.02))
        ax.set_title('Clear surface·실광선 추적 588 nm  |  ×: 추적 중단  |  clear 끝점 연결은 개략도',loc='left',pad=7)
    else:ax.tick_params(labelsize=6);ax.xaxis.label.set_size(6);ax.yaxis.label.set_size(6)
    return traced


def spot_panel(ax,design,field,title):
    # One common physical image plane; all wavelengths remain registered.
    raw=trace_design(design,fields=(field,),rings=8,arms=24)
    allpoints=[];allweights=[]
    for wi,wavelength in enumerate(raw['wavelengths_nm']):
        good=(raw['status'][wi,0]==0)&raw['image_forward'][wi,0]
        weights=raw['weights']*good
        if weights.sum()>0:
            p=raw['image_xy'][wi,0];allpoints.append(p);allweights.append(weights)
    if not allpoints:
        ax.text(.5,.5,'유효 광선 없음',ha='center',va='center',transform=ax.transAxes);return
    centroid=sum(np.sum(p*w[:,None],axis=0) for p,w in zip(allpoints,allweights))/sum(w.sum() for w in allweights)
    largest=0.
    for wi,wavelength in enumerate(raw['wavelengths_nm']):
        good=(raw['status'][wi,0]==0)&raw['image_forward'][wi,0]
        p=(raw['image_xy'][wi,0,good]-centroid)*1000
        if len(p):
            largest=max(largest,float(np.max(abs(p))))
            ax.scatter(p[:,0],p[:,1],s=3.,alpha=.5,color=WAVE_COLORS[float(wavelength)],rasterized=True,linewidths=0,label=f'{wavelength:g}')
    lim=max(1.,largest*1.08);ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
    ax.axhline(0,color='#bbb',lw=.4);ax.axvline(0,color='#bbb',lw=.4)
    ax.set_aspect('equal',adjustable='box');ax.set_title(title,loc='left');ax.set_xlabel('Δx [µm]');ax.set_ylabel('Δy [µm]')
    ax.text(.02,.98,'공통 RGB centroid 기준',transform=ax.transAxes,ha='left',va='top',fontsize=6,color='#666')


def ri_panel(ax,metrics):
    plotted=False
    for wavelength in (486.,588.,656.):
        values=[]
        for f in FIELDS:
            row=get_metric_row(metrics,f,wavelength)
            val=numeric_leaf(row,('relative_illumination_uncoated','RI_uncoated','ri_uncoated','relative_illumination_fresnel'))
            values.append(np.nan if val is None else val)
        if np.isfinite(values).any():ax.plot(FIELDS,values,'o-',ms=3,lw=1,color=WAVE_COLORS[wavelength],label=f'{wavelength:g} nm');plotted=True
    ax.set_title('RI 근사 — 무코팅 Fresnel·chief map',loc='left');ax.set_xlabel('정규화 field');ax.set_ylabel('E(field) / E(0)');ax.set_xlim(0,1);ax.grid(alpha=.16)
    if plotted:ax.legend(fontsize=6,frameon=False)
    else:ax.text(.5,.5,'RI 계산값 미제공',ha='center',va='center',transform=ax.transAxes)


def distortion_panel(ax,metrics):
    plotted=False
    for wavelength in (486.,588.,656.):
        values=[]
        for f in FIELDS:
            row=get_metric_row(metrics,f,wavelength)
            val=numeric_leaf(row,('distortion_percent','rectilinear_distortion_percent','distortion_chief_slope_percent'))
            values.append(np.nan if val is None else val)
        if np.isfinite(values).any():ax.plot(FIELDS,values,'o-',ms=3,lw=1,color=WAVE_COLORS[wavelength],label=f'{wavelength:g}');plotted=True
    ax.axhline(0,color='#aaa',lw=.5);ax.set_title('왜곡 — 근축 EFL·tanθ 기준',loc='left');ax.set_xlabel('정규화 field');ax.set_ylabel('distortion [%]');ax.set_xlim(0,1);ax.grid(alpha=.16)
    if not plotted:ax.text(.5,.5,'왜곡 계산값 미제공',ha='center',va='center',transform=ax.transAxes)


def mtf_panel(ax,wave,field,geometric_proxy=None):
    mono=[r for r in wave.get('monochromatic',[]) if abs(float(r['field_fraction'])-field)<1e-7]
    scalar_available=[r for r in mono if r.get('mtf_t') is not None and r.get('mtf_s') is not None]
    use_geometric=len(scalar_available)<3
    states=[str(r.get('status','unknown')) for r in mono]
    if use_geometric:
        geometric_count=0
        for r in mono:
            g=r.get('geometric_diagnostic');wavelength=float(r['wavelength_nm'])
            if not isinstance(g,dict) or g.get('mtf_t') is None or g.get('mtf_s') is None:continue
            geometric_count+=1;freq=r['frequencies_lp_mm']
            ax.plot(freq,g['mtf_t'],color=WAVE_COLORS.get(wavelength,'gray'),lw=.85,ls='-',label=f'{wavelength:g} 기하 T')
            ax.plot(freq,g['mtf_s'],color=WAVE_COLORS.get(wavelength,'gray'),lw=.85,ls='--',label=f'{wavelength:g} 기하 S')
        poly=next((r for r in wave.get('geometric_polychromatic_diagnostic',[]) if abs(float(r['field_fraction'])-field)<1e-7),None)
        if poly is None and geometric_proxy:
            ff=np.asarray(geometric_proxy.get('fields',[]),float);match=np.flatnonzero(abs(ff-field)<1e-7)
            if len(match):
                i=int(match[0]);poly={'frequencies_lp_mm':geometric_proxy['frequencies_lp_mm'],
                    'mtf_t':geometric_proxy['MTF_t'][i],'mtf_s':geometric_proxy['MTF_s'][i]}
        if poly is not None:
            ax.plot(poly['frequencies_lp_mm'],poly['mtf_t'],color='#111827',lw=1.5,label='RGB 기하 T')
            ax.plot(poly['frequencies_lp_mm'],poly['mtf_s'],color='#111827',lw=1.5,ls='--',label='RGB 기하 S')
        if not geometric_count and poly is None:
            ax.text(.5,.55,'기하·파동 계산값 미제공',ha='center',va='center',transform=ax.transAxes,fontsize=7,color='#a13131')
        flag='파동 일부 적용 불가' if scalar_available else ('파동모델 적용 불가' if mono else '파동값 미제공')
        title=f'{field:g}F ({12*field:g}°) · 기하 MTF 진단\n{flag} — 회절 MTF 아님'
    else:
        for r in mono:
            freq=r['frequencies_lp_mm'];wavelength=float(r['wavelength_nm'])
            ax.plot(freq,r['mtf_t'],color=WAVE_COLORS.get(wavelength,'gray'),lw=.9,ls='-',label=f'{wavelength:g} T')
            ax.plot(freq,r['mtf_s'],color=WAVE_COLORS.get(wavelength,'gray'),lw=.9,ls='--',label=f'{wavelength:g} S')
            fft=r.get('actual_fft_diagnostic')
            if isinstance(fft,dict) and fft.get('mtf_t') is not None:
                step=max(1,len(freq)//8)
                ax.plot(freq[::step],fft['mtf_t'][::step],ls='none',marker='x',ms=3,color='#777',alpha=.7,label='실제 FFT T 진단')
        poly=next((r for r in wave.get('polychromatic',[]) if abs(float(r['field_fraction'])-field)<1e-7),None)
        if poly is not None:
            ax.plot(poly['frequencies_lp_mm'],poly['mtf_t'],color='#111827',lw=1.5,label='RGB T')
            ax.plot(poly['frequencies_lp_mm'],poly['mtf_s'],color='#111827',lw=1.5,ls='--',label='RGB S')
        flag='조건부·격자진단 통과' if all(s=='conditional' for s in states) else '미수렴 진단곡선'
        title=f'{field:g}F ({12*field:g}°) · scalar 파동 MTF\n{flag}'
    ax.set_title(title,loc='left',fontsize=7.5,color='#8f4a28' if use_geometric else '#17354b')
    ax.set_xlim(0,20);ax.set_ylim(0,1.03);ax.set_xticks([0,5,10,15,20]);ax.set_yticks([0,.25,.5,.75,1.]);ax.grid(alpha=.18)
    ax.set_xlabel('공간주파수 [lp/mm]');ax.set_ylabel('MTF')
    if ax.lines:ax.legend(fontsize=5.2,frameon=False,ncol=2,loc='upper right')


def plot_candidate(record):
    cid=record_id(record);d=record['design'];m=record.get('metrics',{});w=normalize_wave(record.get('wave',{}));sm=candidate_summary(record)
    fig=plt.figure(figsize=(14,10.8));gs=fig.add_gridspec(3,4,height_ratios=(1.45,.95,1.15),left=.055,right=.98,top=.87,bottom=.09,hspace=.48,wspace=.37)
    fig.suptitle(f"{cid}  |  {' — '.join(d['materials'])}",x=.055,y=.98,ha='left',fontsize=15,fontweight='bold',color='#18344b')
    mode='EPD 50 mm 조건' if d['aperture_mode']=='epd50' else '첫 렌즈 clear 지름 50 mm 조건'
    firstline=f"{mode} · 실제 EPD {fmt(sm['EPD_mm'],2)} mm · EFL {fmt(sm['EFL_mm'],2)} mm · BFL {fmt(sm['BFL_vertex_mm'],2)} mm · NA(588) {fmt(sm['NA_588'],4)}"
    secondline=f"588 nm · 1F RI {fmt(sm['RI_edge_588'],4)} · 왜곡 {fmt(sm['distortion_edge_588_percent'],3,'%')} · spot RMS 반경 {fmt(None if sm['spot_edge_588_rms_mm'] is None else sm['spot_edge_588_rms_mm']*1000,2)} µm · 판정: {sm['qualification']}"
    fig.text(.055,.943,firstline,fontsize=9);fig.text(.055,.922,secondline,fontsize=9)
    ct=' / '.join(fmt(v,3) for v in sm['center_thickness_mm'])
    fig.text(.055,.901,f'CT [mm] {ct} (각 2 < CT ≤ 20) · 같은 고정 상면 · 고 NA scalar 모델 조건부 결과 · 후보 ID는 품질 순위 아님',fontsize=8,color='#8d4a25')
    plot_layout(fig.add_subplot(gs[0,:3]),d)
    pa=fig.add_subplot(gs[0,3]);pa.axis('off')
    rows=[]
    for surface in export_prescription(d):
        rows.append([str(surface['surface']), 'EA' if surface['type']=='even_asphere' else 'STD', fmt(surface['radius_mm'],2), fmt(surface['thickness_to_next_surface_mm'],2), fmt(surface['clear_semi_diameter_mm'],2)])
    table=pa.table(cellText=rows,colLabels=['S','type','R [mm]','t [mm]','CSD'],loc='upper center',bbox=[0,.52,1,.48],cellLoc='right',colLoc='right',colWidths=[.1,.18,.26,.23,.23])
    table.auto_set_font_size(False);table.set_fontsize(6.5);table.scale(1,1.5)
    for (rr,cc),cell in table.get_celld().items():
        cell.set_edgecolor('#dde5eb')
        if rr==0:cell.set_facecolor('#e8eff5');cell.set_text_props(weight='bold')
    pa.text(0,.47,'EA: even asphere / STD: spherical\nt: 다음 면까지의 축상 거리\nCSD: clear semi-diameter [mm]\n마지막 t: 상면까지의 vertex BFL\nA4–A12와 conic은 후보 JSON 참조',transform=pa.transAxes,fontsize=6.5,va='top',linespacing=1.5)
    pa.text(0,.18,'NA·RI·spot: 저장된 전체 표본 평가\n그림 spot: 가독성을 위한 8×24 표본\n단면 ray: 588 nm meridional\n기구 외곽·edge 연결은 개략도',transform=pa.transAxes,fontsize=6.5,va='top',color='#5b6874',linespacing=1.5)
    spot_panel(fig.add_subplot(gs[1,0]),d,.5,'Spot 0.5F — RGB')
    spot_panel(fig.add_subplot(gs[1,1]),d,1.,'Spot 1F — RGB')
    ri_panel(fig.add_subplot(gs[1,2]),m);distortion_panel(fig.add_subplot(gs[1,3]),m)
    for j,f in enumerate(FIELDS):mtf_panel(fig.add_subplot(gs[2,j]),w,f,record.get('geometric_mtf_proxy'))
    fig.text(.055,.035,'MTF 실선=T / 파선=S. 파동: scalar pupil autocorrelation (FFT 직접 실행과 구별). 기하 진단 패널은 회절 MTF가 아님. RGB: 공통 영상 좌표의 complex OTF 합산.',fontsize=7,color='#555')
    fig.text(.055,.019,'486 / 588 / 656 nm. RI: 균일 원거리 radiance·무코팅 Fresnel·chief-map 근사. 단면의 clear 끝점 연결은 개략도이며 기구·코팅·공차는 미설계.',fontsize=7,color='#555')
    return fig


def write_contact_sheet(records,path):
    columns=4;rows=math.ceil(len(records)/columns);fig,axes=plt.subplots(rows,columns,figsize=(16,rows*3.15),squeeze=False)
    for ax,rec in zip(axes.ravel(),records):
        sm=candidate_summary(rec);plot_layout(ax,rec['design'],compact=True)
        ct='/'.join(fmt(v,2) for v in sm['center_thickness_mm'])
        title=f"{record_id(rec)} · {' / '.join(rec['design']['materials'])}\nCT {ct} mm | EPD {fmt(sm['EPD_mm'],1)} | NA588 {fmt(sm['NA_588'],3)} | RI588 {fmt(sm['RI_edge_588'],3)}"
        ax.set_title(title,loc='left',fontsize=7,pad=5)
    for ax in axes.ravel()[len(records):]:ax.axis('off')
    fig.suptitle(f'CT20 초기 후보 {len(records)}개 전체 — 유리 위치·재료·seed ID 순서 (품질 순위 아님)\n각 2 < CT ≤ 20 mm · 실제 588 nm chief / marginal 광선',fontsize=12,y=.998)
    fig.tight_layout(rect=(0,0,1,.975));fig.savefig(path,dpi=150);plt.close(fig)


def html_report(records,images,summary,path):
    # JSON < escaped so candidate strings cannot terminate a script element.
    contact='data:image/png;base64,'+base64.b64encode((path.parent/'candidate_contact_sheet.png').read_bytes()).decode('ascii')
    data=json.dumps({'summary':summary,'images':images,'records':records,'contact_sheet':contact},ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    doc='''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>3매 광학계 초기 후보 보고서</title>
<style>
:root{--navy:#17354b;--muted:#536677;--border:#dbe3ea}*{box-sizing:border-box}body{margin:0;background:#f3f6f9;color:#172d3e;font:15px/1.55 system-ui,sans-serif}header{padding:28px 32px;background:#fff;border-bottom:1px solid var(--border)}h1{margin:0 0 9px;font-size:27px}p{margin:8px 0}.notice{background:#fff3df;border-left:4px solid #d89124;padding:10px 14px;font-size:13px;max-width:1150px}.grid{display:grid;grid-template-columns:340px minmax(0,1fr);gap:18px;padding:20px}.sidebar,.report{background:#fff;border:1px solid var(--border);border-radius:10px;overflow:hidden}.sidebar{align-self:start;position:sticky;top:12px;max-height:calc(100vh - 25px);overflow:auto}input,select{font:inherit;padding:8px;border:1px solid var(--border);border-radius:5px;width:100%;margin-bottom:8px}.filters{padding:14px}.candidate{padding:11px 14px;border-top:1px solid #edf0f3;cursor:pointer}.candidate:hover,.candidate.active{background:#e9f3fb}.candidate strong{display:block;color:var(--navy)}.candidate small{color:var(--muted);font-size:11px;line-height:1.5;display:block}.report img{display:block;width:100%;height:auto}.toolbar{display:flex;gap:10px;flex-wrap:wrap;padding:15px;border-bottom:1px solid var(--border)}button{border:1px solid #bbcbd8;background:#fff;border-radius:5px;padding:8px 12px;color:var(--navy);cursor:pointer;font:inherit}button:hover{background:#edf4f9}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:8px;padding:15px;background:#f8fbfd}.metric{font-size:11px;color:var(--muted)}.metric b{display:block;font-size:18px;color:var(--navy)}details{margin:15px}pre{font-size:11px;max-height:450px;overflow:auto;white-space:pre-wrap;background:#f7f9fa;padding:12px}footer{padding:0 32px 25px;color:#667;font-size:12px}@media(max-width:850px){.grid{grid-template-columns:1fr}.sidebar{position:static;max-height:280px}.metrics{grid-template-columns:repeat(2,1fr)}header{padding:18px}.grid{padding:12px}}@media print{.sidebar,.toolbar,header,.metrics,footer{display:none}.grid{display:block;padding:0}.report{border:0}}
</style><header><h1>CT20 · 3매 광학계 후보군 전체 보고서</h1><p>무한 물체 · ±12° · 상고 ±6.4 mm · RGB 486 / 588 / 656 nm · 플라스틱 2매 + 구면 유리 1매 · 각 2 mm &lt; CT ≤ 20 mm</p><div class="notice">합격한 독립 후보군 전체를 보존했습니다. 유리 위치·재료·원본 seed ID 순서이며 단일 winner나 품질 순위가 아닙니다. 전체 단면을 먼저 보고 목록에서 각 후보의 RI·왜곡·spot·MTF를 열 수 있습니다. 구경·NA·상고·투과·CT 제약 판정과 MTF 적용/수렴 상태를 구분합니다. 파동 MTF는 scalar pupil autocorrelation이며, 적용 불가 field는 별도 기하 MTF 진단으로 표시합니다. 기하 MTF는 회절 결과가 아니며 고 NA 벡터 회절 또는 Zemax 동치가 검증된 자료가 아닙니다.</div></header>
<main class="grid"><aside class="sidebar"><div class="filters"><input id="search" placeholder="후보 ID, 원본 seed 또는 재료 검색"><select id="mode"><option value="">모든 구경 조건</option><option value="epd50">EPD 50 mm</option><option value="front_diameter50">첫 렌즈 clear 지름 50 mm</option></select><small id="count"></small></div><div id="list"></div></aside><section class="report"><div class="toolbar"><button id="overview">전체 후보 단면</button><button id="prev">← 이전 후보</button><button id="next">다음 후보 →</button><button id="download" disabled>이 후보 JSON 받기</button><button onclick="window.print()">현재 보기 인쇄</button></div><div class="metrics" id="metrics"></div><img id="figure" alt="후보군 전체의 실제 chief 및 marginal 광선 단면"><details id="details" hidden><summary>처방 및 계산 원자료 보기</summary><pre id="raw"></pre></details></section></main><footer>외부 JavaScript·폰트·이미지 요청이 없는 오프라인 보고서입니다. 그림의 각 곡선은 저장된 계산 결과이며, 없는 성능값을 생성하지 않습니다. RI는 균일 원거리 radiance와 무코팅 표면의 국소 영상 mapping 모델입니다.</footer>
<script id="report-data" type="application/json">__DATA__</script><script>
const D=JSON.parse(document.getElementById('report-data').textContent);let selected=-1,visible=[];const $=id=>document.getElementById(id);const fmt=(v,n=3)=>v===null||v===undefined?'—':Number(v).toFixed(n);function el(tag,text){const e=document.createElement(tag);e.textContent=text;return e}
function cards(values){$('metrics').replaceChildren();values.forEach(([k,v])=>{const t=el('div',k);t.className='metric';t.append(el('b',v));$('metrics').append(t)})}
function overview(){selected=-1;$('figure').src=D.contact_sheet;$('figure').alt='합격한 전체 '+D.records.length+'개 후보의 실제 chief 및 marginal 광선 단면';$('raw').textContent='';$('details').hidden=true;$('download').disabled=true;cards([['전체 독립 후보',D.records.length+'개'],['각 렌즈 중심 두께','2 < CT ≤ 20 mm'],['배치 순서','유리 위치 · 재료 · seed'],['성능 순위','설정하지 않음']]);document.querySelectorAll('.candidate').forEach(e=>e.classList.remove('active'))}
function select(i){selected=i;const s=D.summary[i];$('figure').src=D.images[i];$('figure').alt=s.candidate_id+' 광선 단면, 상대조도, 왜곡, spot, MTF';$('raw').textContent=JSON.stringify(D.records[i],null,2);$('details').hidden=false;$('download').disabled=false;cards([['후보',s.candidate_id],['CT 3매 [mm]',s.center_thickness_mm.map(v=>fmt(v,3)).join(' / ')],['실제 EPD [mm]',fmt(s.EPD_mm,2)],['588 nm NA',fmt(s.NA_588,4)],['1F RI (588 nm)',fmt(s.RI_edge_588,4)],['EFL / BFL [mm]',fmt(s.EFL_mm,2)+' / '+fmt(s.BFL_vertex_mm,2)],['1F 왜곡 (588 nm) [%]',fmt(s.distortion_edge_588_percent)],['1F spot RMS 반경 (588 nm) [µm]',fmt(s.spot_edge_588_rms_mm===null?null:s.spot_edge_588_rms_mm*1000,2)]]);document.querySelectorAll('.candidate').forEach(e=>e.classList.toggle('active',+e.dataset.index===i));}
function filter(){const q=$('search').value.toLowerCase(),mode=$('mode').value;visible=[];$('list').replaceChildren();D.summary.forEach((s,i)=>{if((s.candidate_id+' '+s.source_seed_id+' '+s.materials).toLowerCase().includes(q)&&(!mode||s.aperture_mode===mode)){visible.push(i);const e=el('div','');e.className='candidate';e.dataset.index=i;e.append(el('strong',s.candidate_id+' · '+s.materials),el('small','CT '+s.center_thickness_mm.map(v=>fmt(v,2)).join('/')+' mm'),el('small','EPD '+fmt(s.EPD_mm,1)+' mm · NA '+fmt(s.NA_588,3)+' · RI '+fmt(s.RI_edge_588,3)),el('small',s.qualification));e.onclick=()=>select(i);$('list').append(e)}});$('count').textContent=visible.length+' / '+D.summary.length+'개 후보';$('prev').disabled=!visible.length;$('next').disabled=!visible.length;if(selected>=0&&visible.includes(selected))select(selected);else overview();}
$('search').oninput=filter;$('mode').onchange=filter;$('overview').onclick=overview;$('prev').onclick=()=>select(visible[selected<0?visible.length-1:(visible.indexOf(selected)-1+visible.length)%visible.length]);$('next').onclick=()=>select(visible[selected<0?0:(visible.indexOf(selected)+1)%visible.length]);$('download').onclick=()=>{if(selected<0)return;const blob=new Blob([JSON.stringify(D.records[selected],null,2)],{type:'application/json'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=D.summary[selected].candidate_id+'.json';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)};filter();
</script></html>'''
    path.write_text(doc.replace('__DATA__',data),encoding='utf-8')


def render(records,output):
    if not records or not all(eligible(r) for r in records):
        raise ValueError('Report accepts only candidates with all CT20 and sampled optical metric gates passed')
    output=Path(output);output.mkdir(parents=True,exist_ok=True);plots=output/'plots';plots.mkdir(exist_ok=True)
    images=[];summary=[];errors=[]
    with PdfPages(output/'candidate_layouts_and_mtf.pdf') as pdf:
        for record in records:
            cid=record_id(record);fig=plot_candidate(record)
            target=plots/(re.sub(r'[^a-zA-Z0-9_.-]','_',cid)+'.png')
            fig.savefig(target,dpi=150);pdf.savefig(fig);plt.close(fig)
            images.append('data:image/png;base64,'+base64.b64encode(target.read_bytes()).decode('ascii'));summary.append(candidate_summary(record))
    write_contact_sheet(records,output/'candidate_contact_sheet.png')
    html_report(records,images,summary,output/'candidate_report.html')
    (output/'report_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    (output/'render_manifest.json').write_text(json.dumps({'count':len(records),'candidate_ids':[record_id(r) for r in records],
      'files':['candidate_layouts_and_mtf.pdf','candidate_contact_sheet.png','candidate_report.html','report_summary.json'],
      'rays':'Fresh traces from stored design; chief and marginal kept, failed traces marked.',
      'MTF':'Saved values only; autocorrelation, actual FFT and geometric diagnostics named separately.',
      'spot_plot_sampling':'8 rings x24 arms pluschief/rim for visual readability; summaryvaluesfromfullmetrics.',
      'errors':errors},ensure_ascii=False,indent=2)+'\n')
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path,default=HERE/'reports');args=parser.parse_args()
    obj=json.loads(args.input.read_text());records=obj.get('candidates',obj.get('records',obj)) if isinstance(obj,dict) else obj
    if not isinstance(records,list):records=[records]
    render(records,args.output)
    print(json.dumps({'output':str(args.output),'candidate_count':len(records)}))
