"""Independent original-v0.2.1 intersection/Snell comparison, not a speed test."""
from pathlib import Path
import json,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(HERE.parent))
from native.raytrace import trace
from sequential_optics.models import Surface,Ray
from sequential_optics.geometry import intersect,gradient,refract
from sequential_optics.errors import GeometryError,NumericalError,TotalInternalReflection

def reference(surfaces,o,d,indices):
    ray=Ray(o,d);p=ray.origin.copy();d=ray.direction.copy();opl=0;points=[];dirs=[];opls=[]
    for j,row in enumerate(surfaces):
        surface=Surface(z=row[0],curvature=row[1],conic=row[2],coefficients=(0.,*row[3:8]),clear_radius=row[8],surface_type='even_asphere')
        try:
            point,distance=intersect(surface,Ray(p,d));g=gradient(surface,*point[:2])
            nd=refract(d,np.array([-g[0],-g[1],1]),indices[j],indices[j+1]) if indices[j]!=indices[j+1] else d.copy()
            if nd[2]<=1e-13:return 4,j,np.array(points),np.array(dirs),np.array(opls)
        except TotalInternalReflection:return 2,j,np.array(points),np.array(dirs),np.array(opls)
        except (GeometryError,NumericalError) as ex:
            code=1 if str(ex).startswith('escaped:') else 3
            return code,j,np.array(points),np.array(dirs),np.array(opls)
        opl+=indices[j]*distance;p=point;d=nd;points.append(p);dirs.append(d);opls.append(opl)
    return 0,-1,np.array(points),np.array(dirs),np.array(opls)

def main():
    rng=np.random.default_rng(20261022);cases=[]
    for curvature,conic,coeff in [(0,0,[0]*5),(.03,0,[0]*5),(-.028,-1,[0]*5),(.038,-1.5,[-2e-6,3e-10,-2e-14,1e-18,-2e-23]),(0,0,[3e-6,-4e-9,2e-12,-1e-15,2e-19])]:
        s=np.array([[0,curvature,conic,*coeff,20.]])
        for _ in range(60):
            o=np.r_[rng.uniform(-15,15,2),-40.];d=np.r_[rng.uniform(-.25,.25,2),1.]
            cases.append((s,o,d,np.array([1.,1.492])))
    for angle in [0,10,30,40,41.8103149,45,50,70]:
        a=np.deg2rad(angle);cases.append((np.array([[0,0,0,0,0,0,0,0,100.]]),np.array([0,0,-1.]),np.array([np.sin(a),0,np.cos(a)]),np.array([1.5,1.])))
    s=np.zeros((6,9));s[:,0]=[0,12,15,23,26,38];s[:,1]=[1/60,-1/60,-1/80,1/80,1/45,-1/45];s[:,8]=20
    s[0,2]=-1;s[0,3]=-3e-7;s[5,2]=-1;s[5,3]=2e-7
    for _ in range(100):cases.append((s,np.r_[rng.uniform(-12,12,2),-30.],np.r_[rng.uniform(-.12,.12,2),1.],np.array([1.,1.492,1.,1.94594,1.,1.586,1.])))
    max_pos=max_dir=max_opl=0.;mismatches=[];codes={};started=time.perf_counter()
    for index,(s,o,d,ns) in enumerate(cases):
        actual=trace(s,o[None],d[None],ns);ref=reference(s,o,d,ns)
        status=int(actual['status'][0]);failed=int(actual['failed_surface'][0]);codes[status]=codes.get(status,0)+1
        if (status,failed)!=ref[:2]:mismatches.append({'case':index,'actual':[status,failed],'reference':ref[:2]});continue
        length=len(ref[2])
        if length:
            max_pos=max(max_pos,float(np.max(np.abs(actual['history'][0,:length]-ref[2]))))
            max_dir=max(max_dir,float(np.max(np.abs(actual['direction_history'][0,:length]-ref[3]))))
            max_opl=max(max_opl,float(np.max(np.abs(actual['opl_history'][0,:length]-ref[4]))))
    report={'cases':len(cases),'status_counts':codes,'status_mismatches':mismatches,'max_position_difference_mm':max_pos,'max_direction_difference':max_dir,'max_opl_difference_mm':max_opl,'seconds':time.perf_counter()-started,'passed':not mismatches and max_pos<1e-8 and max_dir<1e-9 and max_opl<1e-8}
    (HERE/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    assert report['passed']

if __name__=='__main__':main()
