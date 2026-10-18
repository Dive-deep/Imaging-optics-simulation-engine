"""Matched finite Lambertian source cell -> physical scalar wave -> image kernel.

This closes the radiometry/wave amplitude interface on a *regular sampled ray
map*. Successful original STOP triangles define conservative support; source
solid angle is never inferred from equal STOP quadrature weights. Each triangle
uses a constant physical boundary amplitude and a mean transmission. Refinement
is required. This is scalar direct-path diffraction, not vector Maxwell optics.
"""
import numpy as np
from .errors import InputError, MissingDataError, ModelApplicabilityError
from .models import Field,Ray
from .settings import EngineSettings
from .wave import (make_boundary_state,_mapped_interpolation,fft_psf,afocal_fft,
                   sample_axis_otf,angular_otf,focal_direct_reference,VACUUM_IMPEDANCE_OHM,PROFILE_ID)
from .radiometry import trace_transmission,_absorption,spherical_bin_solid_angles
from .imaging import SourceKernel


def spherical_triangle_measures(directions):
    """Signed exact solid angle and projected solid angle of unit great-circle triangles.

    directions has shape (...,3,3), ordered vertices. The z projection is
    integral cos(theta_source) dOmega for a +z Lambertian source. The vector
    area identity is 1/2 integral r cross dr around the spherical boundary.
    This describes geodesic triangle edges, which approximate curved ray-map
    edges as the original STOP mesh is refined.
    """
    d=np.asarray(directions,float)
    if d.shape[-2:]!=(3,3) or not np.all(np.isfinite(d)):
        raise InputError('Direction triangles must have trailing shape (3,3)')
    norm=np.linalg.norm(d,axis=-1)
    if np.any(norm==0): raise InputError('Direction vertices must be nonzero')
    d=d/norm[...,None]
    if np.any(d[...,2]<=0): raise ModelApplicabilityError('Source directions must lie in forward hemisphere')
    a,b,c=d[...,0,:],d[...,1,:],d[...,2,:]
    solid=2*np.arctan2(np.sum(a*np.cross(b,c),axis=-1),
                     1+np.sum(a*b+b*c+c*a,axis=-1))
    vector=np.zeros(d.shape[:-2]+(3,))
    for p,q in ((a,b),(b,c),(c,a)):
        cross=np.cross(p,q);sine=np.linalg.norm(cross,axis=-1)
        angle=np.arctan2(sine,np.sum(p*q,axis=-1))
        factor=np.divide(angle,sine,out=np.ones_like(angle),where=sine>0)
        vector+=cross*factor[...,None]/2
    return solid,vector[...,2]


def _triangle_area(uv,tri):
    p=uv[tri];a=p[:,1]-p[:,0];b=p[:,2]-p[:,0]
    return abs(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2


def _raster(state,mode,image_z,N,grid_margin,triangle_amplitude,direction_cosines):
    uv=state.uv_mm;good=state.valid_ray_mask;n=state.output_index
    if mode=='focal':
        chief=state.chief_position;d=state.chief_direction
        q0=chief+(image_z-chief[2])/d[2]*d
        L=float(np.dot(q0-state.origin,d))
        if L<=0: raise ModelApplicabilityError('Image plane must be downstream of physical wave boundary')
        R=np.sqrt(L*L+np.sum(uv*uv,axis=1));vertices=uv/R[:,None]
        opd=state.opl_mm-state.chief_opl_mm+n*np.sum(uv*uv,axis=1)/(R+L)
    else:
        L=None;q0=None;vertices=uv;opd=state.opl_mm-state.chief_opl_mm
    ext=np.max(np.abs(vertices[good]),axis=0)*grid_margin
    if np.any(ext<=0): raise ModelApplicabilityError('Degenerate boundary pupil')
    spacing=2*ext/N
    x=(np.arange(N)-N//2)*spacing[0];y=(np.arange(N)-N//2)*spacing[1]
    cosine=np.full(len(uv),np.nan)
    # Actual ray normal cosines are stored separately by the bridge.
    cosine[good]=direction_cosines[good]
    values=np.column_stack([opd,cosine])
    vals,mask,owner=_mapped_interpolation(vertices,state.triangles,values,x,y)
    amplitude=np.zeros(mask.shape);amplitude[mask]=triangle_amplitude[owner[mask]]
    xx,yy=np.meshgrid(x,y)
    if mode=='focal':
        g2=1-xx*xx-yy*yy
        if np.any(mask&(g2<=0)): raise ModelApplicabilityError('Focal pupil outside reference hemisphere')
        g2=np.where(mask,g2,1.)
        jac=L*L/g2**2
        fourier_amplitude=amplitude*L/g2
    else:
        jac=np.ones(mask.shape);fourier_amplitude=amplitude
    phase=2*np.pi*vals[:,:,0]/(state.wavelength_nm*1e-6)
    field=np.where(mask,fourier_amplitude*np.exp(1j*phase),0j)
    boundary_flux=float(np.sum(amplitude**2*vals[:,:,1]*jac)*spacing.prod()*state.output_index/(2*VACUUM_IMPEDANCE_OHM)*1e-6)
    max_step=0.
    for axis in (0,1):
        m=np.take(mask,range(N-1),axis=axis)&np.take(mask,range(1,N),axis=axis)
        if np.any(m): max_step=max(max_step,float(np.max(abs(np.diff(phase,axis=axis)[m]))))
    triangle_span=float(np.max(np.ptp(opd[state.triangles],axis=1))*2*np.pi/(state.wavelength_nm*1e-6))
    return dict(field=field,mask=mask,owner=owner,x=x,y=y,spacing=spacing,
                opd_mm=vals[:,:,0],boundary_amplitude_v_m=amplitude,
                boundary_grid_flux_w=boundary_flux,max_grid_phase_step_rad=max_step,
                max_triangle_phase_span_rad=triangle_span,distance_mm=L,image_reference_point_mm=q0)


def source_diffraction_kernel(system,source_xy_mm,source_area_mm2,source_radiance=None,*,
                              fill_factor=1.,bundle=None,wavelength_nm=None,settings=None,
                              mode='focal',image_z_mm=None,image_grid=None,
                              grid_size=None,padding=None,frequencies=None,
                              transmission_mode='unpolarized_jones_power',missing_absorption='transparent',
                              raw_power_grid_size=None,grid_margin=1.05,
                              phase_step_limit_rad=np.pi,allow_local_detector_projection=False,
                              direct_check_points=None):
    """Create matched Q and scalar diffraction for one finite source macro-cell.

    ``source_area_mm2`` is the actual quadrature cell area, not the emitter's
    whole area. Radiance is active-area W/(mm2 sr); fill is applied exactly once.
    With radiance=None, amplitudes/powers use a documented relative L=1 scale.
    Q is collected power per source area. Optional ``image_grid`` produces a
    cropped SourceKernel on global detector coordinates without renormalizing.

    Off-axis chief-normal -> z-detector conversion is an additional local
    parallel-projection approximation and requires explicit opt-in. For exact
    detector propagation at large tilt/NA use a dedicated scalar/vector solver.
    """
    from .geometry import trace_bundle,gradient,refract
    settings=settings or EngineSettings();wl=settings.wavelength(wavelength_nm)
    xy=np.asarray(source_xy_mm,float)
    if xy.shape!=(2,) or not np.all(np.isfinite(xy)): raise InputError('source_xy_mm must be finite x,y')
    if system.object_z is None: raise MissingDataError('Finite physical object_z required by source diffraction')
    if not np.isfinite(source_area_mm2) or source_area_mm2<=0: raise InputError('Positive source cell area required')
    Lsource=1. if source_radiance is None else float(source_radiance)
    if not np.isfinite(Lsource) or Lsource<=0 or not np.isfinite(fill_factor) or not 0<fill_factor<=1:
        raise InputError('Source diffraction requires positive radiance and fill in (0,1]')
    if mode not in ('focal','afocal'): raise InputError('mode must be focal or afocal')
    if mode!='focal' and (image_grid is not None or direct_check_points is not None):
        raise InputError('Focal detector grid and near-focus RS points require focal mode')
    N=settings.sampling.wave_grid_size if grid_size is None else grid_size
    pad=settings.sampling.fft_padding if padding is None else padding
    Nraw=N if raw_power_grid_size is None else raw_power_grid_size
    if any(int(v)!=v or v<8 or v%2 for v in (N,Nraw)) or grid_margin<=1:
        raise InputError('Wave grids must be even integers >=8 and grid_margin >1')
    z=system.image_z if image_z_mm is None else image_z_mm
    if mode=='focal' and z is None: raise MissingDataError('Focal source diffraction needs image_z')
    if bundle is None:
        bundle=trace_bundle(system,Field(float(xy[0]),float(xy[1]),'height'),wl,settings=settings)
    elif (bundle.field.kind!='height' or not np.allclose([bundle.field.x,bundle.field.y],xy,rtol=0,atol=settings.ray.intersection_tol_mm)
          or bundle.wavelength_nm!=wl):
        raise InputError('Source wave bundle must match finite source point and requested wavelength')
    state=make_boundary_state(system,bundle,settings,image_z=z,mode=mode)
    good=state.valid_ray_mask;tri=state.triangles
    source_point=np.array([xy[0],xy[1],system.object_z])
    directions=np.full((len(bundle.traces),3),np.nan);cosine=np.full(len(bundle.traces),np.nan)
    tau=np.full(len(bundle.traces),np.nan);assumptions=set()
    if transmission_mode!='lossless_scalar':
        mu_out=_absorption(system,system.surfaces[-1].medium_after,wl,missing_absorption,assumptions)
        if mu_out!=0:
            raise ModelApplicabilityError('Scalar source-wave propagation requires nonabsorbing output medium; complex-k propagation is not implemented')
    for i,t in enumerate(bundle.traces):
        if not good[i]: continue
        initial=t.positions[0]-source_point;initial/=np.linalg.norm(initial)
        first=system.surfaces[0];g=gradient(first,*t.positions[0,:2])
        normal=np.array([-g[0],-g[1],1.]);normal/=np.linalg.norm(normal)
        n0=system.index(system.ambient_material,wl);n1=system.index(first.medium_after,wl)
        expected_direction=refract(initial,normal,n0,n1)
        expected_opl=n0*np.linalg.norm(t.positions[0]-source_point)
        tolerance=max(settings.ray.intersection_tol_mm*n0,
                      np.finfo(float).eps*max(1.,abs(expected_opl)))*10
        if (np.linalg.norm(expected_direction-t.directions[0])>settings.ray.direction_tol*10
                or abs(expected_opl-t.opl[0])>tolerance):
            raise InputError('Provided trace does not originate at the declared physical source point with its finite-source OPL')
        directions[i]=initial;cosine[i]=np.dot(t.last_direction,state.chief_direction)
        transport=trace_transmission(system,Ray(source_point,initial),t,wl,settings,
                                     transmission_mode,missing_absorption)
        tau[i]=transport['transmission'];assumptions.update(transport['assumptions'])
    omega,projected=spherical_triangle_measures(directions[tri])
    signs=np.sign(projected)
    if np.any(signs==0) or np.any(signs!=signs[0]) or np.any(np.sign(omega)!=signs):
        raise ModelApplicabilityError('Launch direction map folds or degenerates on sampled STOP triangles')
    projected=np.abs(projected);omega=np.abs(omega)
    if projected.sum()>np.pi*(1+1e-10):
        raise ModelApplicabilityError('Sampled source map covers more than the Lambertian hemisphere')
    q_triangle=Lsource*fill_factor*projected*np.mean(tau[tri],axis=1)
    area=_triangle_area(state.uv_mm,tri)
    normal_cos=np.mean(cosine[tri],axis=1)
    if np.any(area<=0) or np.any(normal_cos<=0):
        raise ModelApplicabilityError('Degenerate or grazing physical boundary triangles')
    p_triangle=q_triangle*source_area_mm2/area
    amplitude=np.sqrt(2*VACUUM_IMPEDANCE_OHM*1e6*p_triangle/(state.output_index*normal_cos))
    raster=_raster(state,mode,z,int(N),grid_margin,amplitude,cosine)
    lam=wl*1e-6;n=state.output_index
    if mode=='focal':
        wave=fft_psf(raster['field'],raster['spacing'],lam,index=n,padding=pad)
        dl=fft_psf(abs(raster['field']),raster['spacing'],lam,index=n,padding=pad)
        f=np.linspace(0,settings.sampling.max_frequency_lp_mm,settings.sampling.frequency_samples) if frequencies is None else np.asarray(frequencies,float)
        iy,ix=np.where(raster['mask']);cutoffs=np.array([ix.max()-ix.min()+1,iy.max()-iy.min()+1])*raster['spacing']/(lam/n)
        wave['frequencies_lp_mm']=f
        for axis,label,cutoff in [('x','s',cutoffs[0]),('y','t',cutoffs[1])]:
            wave['otf_'+label]=sample_axis_otf(wave,f,axis,cutoff)
            wave['mtf_'+label]=abs(wave['otf_'+label])
            wave['diffraction_mtf_'+label]=abs(sample_axis_otf(dl,f,axis,cutoff))
        wave['same_pupil_diffraction_psf']=dl['psf']
        cc,ee=np.meshgrid(raster['x'],raster['y']);mask=raster['mask']
        L=raster['distance_mm'];gamma=np.sqrt(1-cc[mask]**2-ee[mask]**2)
        u=L*cc[mask]/gamma;v=L*ee[mask]/gamma;R=L/gamma
        phase_residual=0.;amplitude_residual=0.
        for qx in (wave['x_mm'][0],wave['x_mm'][-1]):
            for qy in (wave['y_mm'][0],wave['y_mm'][-1]):
                rq=np.sqrt((u-qx)**2+(v-qy)**2+L*L)
                phase_residual=max(phase_residual,float(np.max(abs(rq-R+cc[mask]*qx+ee[mask]*qy)))*2*np.pi*n/lam)
                amplitude_residual=max(amplitude_residual,float(np.max(abs(R*R/(rq*rq)-1))))
        wave['local_expansion_diagnostics']={'full_window_phase_residual_rad':phase_residual,
                'full_window_relative_amplitude_residual':amplitude_residual,'not_an_mtf_error_bound':True}
        if direct_check_points is not None:
            points=np.atleast_2d(np.asarray(direct_check_points,float))
            exact=focal_direct_reference(raster['field'],raster['x'],raster['y'],L,lam,points,index=n)
            local=np.array([np.sum(raster['field']*np.exp(-2j*np.pi*n/lam*(qx*cc+qy*ee)))*raster['spacing'].prod()/(1j*(lam/n)) for qx,qy in points])
            wave['direct_scalar_check']={'points_xy_mm':points,'exact_intensity':abs(exact)**2,
                    'local_fft_integral_intensity':abs(local)**2,'shared_boundary_input':True}
    else:
        wave=afocal_fft(raster['field'],raster['spacing'],lam,index=n,padding=pad)
        dl=afocal_fft(abs(raster['field']),raster['spacing'],lam,index=n,padding=pad)
        wave['same_pupil_diffraction_h_solid_angle']=dl['h_solid_angle']
        if frequencies is not None:
            f=np.asarray(frequencies,float);wave['frequencies_cycles_rad']=f
            for axis,label in [('x','s'),('y','t')]:
                wave['otf_'+label]=angular_otf(wave,f,axis=axis);wave['mtf_'+label]=abs(wave['otf_'+label])
    # Raw wave/boundary/ray comparison uses a separate uniform PHYSICAL-plane
    # field and its angular spectrum even when the requested output is focal.
    # Focal reference phase removal cannot be reused for this power diagnostic.
    raw_raster=raster if mode=='afocal' and Nraw==N else _raster(state,'afocal',None,int(Nraw),grid_margin,amplitude,cosine)
    raw_wave=afocal_fft(raw_raster['field'],raw_raster['spacing'],lam,index=n,padding=pad)
    raw_resolved=raw_raster['max_grid_phase_step_rad']<phase_step_limit_rad and raw_raster['max_triangle_phase_span_rad']<phase_step_limit_rad
    q=float(q_triangle.sum());ray_power=q*source_area_mm2
    raw_wave_power=raw_wave['electric_field_wave_power_w']
    phase_ok=max(raster['max_grid_phase_step_rad'],raster['max_triangle_phase_span_rad'])<phase_step_limit_rad
    wave.update(profile_id=PROFILE_ID+'_source_triangle',mode=mode,
                status='conditional' if phase_ok and wave.get('linear_autocorrelation_supported',True) else 'underresolved',
                vendor_equivalence='unverified',illumination='physical_source_triangle_boundary_amplitude',
                output_basis=state.basis,boundary_origin_mm=state.origin,
                image_reference_point_mm=raster['image_reference_point_mm'],distance_mm=raster['distance_mm'],
                pupil=raster['field'],pupil_x=raster['x'],pupil_y=raster['y'],pupil_mask=raster['mask'],
                pupil_opd_mm=raster['opd_mm'],pupil_spacing=raster['spacing'],
                axis_convention='S=projected object X; T=e3 cross e1; not radial/meridional',
                removed_phase_terms=['piston'],
                model_validity='scalar sampled single-branch boundary; convergence and near-focus approximation unverified',
                sampling_diagnostics={'max_grid_phase_step_rad':raster['max_grid_phase_step_rad'],
                                      'max_triangle_phase_span_rad':raster['max_triangle_phase_span_rad'],
                                      'support':'conservative_transmitted_STOP_triangles',
                                      'amplitude_model':'constant physical electric-field amplitude per mapped triangle',
                                      'transmission_quadrature':'vertex_mean within spherical source triangle',
                                      'source_map_edge_model':'great_circle chords; refine original STOP mesh'})
    result=dict(wave=wave,collected_power_density=q,source_cell_area_mm2=float(source_area_mm2),
                source_cell_collected_power=ray_power,source_xy_mm=xy,wavelength_nm=wl,
                source_hemisphere_power_density=np.pi*Lsource*fill_factor,
                source_projected_solid_angle_sr=float(projected.sum()),
                boundary_triangle_amplitude_v_m=amplitude,boundary_triangle_area_mm2=area,
                boundary_triangle_power_w=q_triangle*source_area_mm2,
                triangle_mean_direction_cosine=normal_cos,
                raw_power_diagnostics={'ray_triangle_power_w':ray_power,
                     'boundary_grid_power_w':raw_raster['boundary_grid_flux_w'],
                     'boundary_grid_to_ray_ratio':raw_raster['boundary_grid_flux_w']/ray_power,
                     'wave_power_w_unvalidated':raw_wave_power,
                     'wave_to_boundary_ratio_unvalidated':raw_wave_power/raw_raster['boundary_grid_flux_w'],
                     'raw_wave_phase_sampling_status':'necessary_check_passed' if raw_resolved else 'underresolved',
                     'raw_max_grid_phase_step_rad':raw_raster['max_grid_phase_step_rad'],
                     'raw_max_triangle_phase_span_rad':raw_raster['max_triangle_phase_span_rad'],
                     'wave_window_edge_power_fraction':raw_wave['edge_energy_fraction'],
                     'note':'Ratios precede PSF normalization; triangle flux itself is by construction. Passing phase check is necessary, not sufficient.'},
                assumptions=sorted(assumptions),units='W, W/mm2_source' if source_radiance is not None else 'relative source scale L=1',
                status=wave['status'],boundary_state_geometry=state)
    if image_grid is not None:
        if mode!='focal': raise InputError('Image PixelGrid is a focal detector; use bin_source_angular_response for afocal')
        kernel,projection=wave_to_source_kernel(result,image_grid,
                              allow_local_detector_projection=allow_local_detector_projection)
        result['kernel']=kernel;result['detector_projection']=projection
    return result


def _clip_polygon(poly,axis,value,greater):
    if not len(poly): return poly
    out=[];prev=poly[-1];inside_prev=prev[axis]>=value if greater else prev[axis]<=value
    for point in poly:
        inside=point[axis]>=value if greater else point[axis]<=value
        if inside!=inside_prev:
            t=(value-prev[axis])/(point[axis]-prev[axis]);out.append(prev+t*(point-prev))
        if inside: out.append(point)
        prev=point;inside_prev=inside
    return np.asarray(out)


def _polygon_area(poly):
    if len(poly)<3: return 0.
    return abs(np.dot(poly[:,0],np.roll(poly[:,1],-1))-np.dot(poly[:,1],np.roll(poly[:,0],-1)))/2


def wave_to_source_kernel(source_result,image_grid,*,allow_local_detector_projection=False,tilt_tolerance=1e-12):
    """Power-conservative overlap of FFT cells onto a global z detector grid.

    At nonzero chief tilt the affine projection is an extra *local parallel-ray*
    approximation, not exact diffraction propagation to a tilted detector.
    Cropped tails remain missing probability mass. No source Q, fill or tau is
    multiplied here; SourceKernel mass will later be multiplied by matched Q.
    """
    w=source_result['wave']
    if w['mode']!='focal' or w['illumination']!='physical_source_triangle_boundary_amplitude':
        raise InputError('A matched physical-source focal wave result is required')
    basis=w['output_basis'];tilted=np.linalg.norm(basis[2,:2])>tilt_tolerance
    if tilted and not allow_local_detector_projection:
        raise ModelApplicabilityError('Off-axis detector conversion needs explicit allow_local_detector_projection=True or exact detector propagation')
    J=basis[:2,:2]
    if abs(np.linalg.det(J))<tilt_tolerance:
        raise ModelApplicabilityError('Singular chief-parallel detector projection')
    inverse=np.linalg.inv(J);origin=w['image_reference_point_mm'][:2]
    dx=w['dx_mm'];dy=w['dy_mm'];cell_area=dx*dy/abs(np.linalg.det(J))
    mass=w['psf']*dx*dy
    output=np.zeros(image_grid.shape)
    xe=image_grid.x_edges_mm;ye=image_grid.y_edges_mm
    offsets=np.array([[-dx/2,-dy/2],[dx/2,-dy/2],[dx/2,dy/2],[-dx/2,dy/2]])@inverse.T
    for iy,y in enumerate(w['y_mm']):
        for ix,x in enumerate(w['x_mm']):
            weight=mass[iy,ix]
            if weight==0: continue
            center=origin+inverse@np.array([x,y]);poly=center+offsets
            lo=poly.min(axis=0);hi=poly.max(axis=0)
            if hi[0]<=xe[0] or lo[0]>=xe[-1] or hi[1]<=ye[0] or lo[1]>=ye[-1]: continue
            i0=max(0,np.searchsorted(xe,lo[0],side='right')-1);i1=min(len(xe)-2,np.searchsorted(xe,hi[0],side='left'))
            j0=max(0,np.searchsorted(ye,lo[1],side='right')-1);j1=min(len(ye)-2,np.searchsorted(ye,hi[1],side='left'))
            for j in range(j0,j1+1):
                for i in range(i0,i1+1):
                    if lo[0]>=xe[i] and hi[0]<=xe[i+1] and lo[1]>=ye[j] and hi[1]<=ye[j+1]:
                        fraction=1.
                    else:
                        clipped=poly
                        for axis,value,greater in ((0,xe[i],True),(0,xe[i+1],False),(1,ye[j],True),(1,ye[j+1],False)):
                            clipped=_clip_polygon(clipped,axis,value,greater)
                            if len(clipped)==0: break
                        fraction=min(1.,max(0.,_polygon_area(clipped)/cell_area))
                    output[j,i]+=weight*fraction
    # Round-off can exceed one by a few ulps; preserve crop power and only
    # remove a demonstrably numerical positive excess, never a physical tail.
    if output.sum()>1 and output.sum()<1+1e-10: output/=output.sum()
    kernel=SourceKernel(image_grid,output,source_result['wavelength_nm'],
                        model='scalar_source_fft_chief_parallel_detector_projection' if tilted else 'scalar_source_fft_on_axis_detector')
    return kernel,dict(model='chief_parallel_local_affine' if tilted else 'coincident_normal_planes',
                       approximation=bool(tilted),jacobian_virtual_from_detector=J,
                       captured_fraction=kernel.captured_fraction,
                       resampling='cell_overlap_power_conserving_no_crop_renormalization')


def bin_source_angular_response(source_result,theta_edges,phi_edges):
    """Rotate local angular PSF directions into global spherical bins, preserving power."""
    w=source_result['wave']
    if w['mode']!='afocal': raise InputError('An afocal source response is required')
    l,m=np.meshgrid(w['l'],w['m']);valid=w['propagating_mask']
    local=np.column_stack([l[valid],m[valid],w['gamma'][valid]])
    global_direction=local@w['output_basis']
    theta=np.arccos(np.clip(global_direction[:,2],-1,1))
    phi=(np.arctan2(global_direction[:,1],global_direction[:,0])-phi_edges[0])%(2*np.pi)+phi_edges[0]
    powers=w['cell_power'][valid]*source_result['source_cell_collected_power']
    solid=spherical_bin_solid_angles(theta_edges,phi_edges)
    histogram=np.histogram2d(theta,phi,bins=(theta_edges,phi_edges),weights=powers)[0]
    total=source_result['source_cell_collected_power']
    return dict(bin_power=histogram,solid_angle_sr=solid,intensity=histogram/solid,
                source_cell_collected_power=total,outside_roi_power=max(0.,total-float(histogram.sum())),
                frame='global_system',registration='full_unit_direction_rotation',
                model='scalar_source_diffraction_cell',status=source_result['status'])
