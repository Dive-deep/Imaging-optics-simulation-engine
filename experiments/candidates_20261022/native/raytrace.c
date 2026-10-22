/* Candidate-search sequential real-ray kernel. All optical arithmetic is double.
 * Surface row: z,c,k,A4,A6,A8,A10,A12,clear_radius.
 * First-hit isolation uses interval sag/derivative bounds, not nearest Newton root.
 * No fast-math; uncertain or tangent first intervals are numerical failures.
 */
#include <math.h>
#include <stddef.h>
#include <float.h>

#define TOL 1e-10
#define MAX_INTERVALS 4096
#define MAX_DEPTH 48

static double mn(double a,double b){return a<b?a:b;}
static double mx(double a,double b){return a>b?a:b;}
static double norm3(const double *x){return hypot(hypot(x[0],x[1]),x[2]);}

static int sag_grad(const double *s,double x,double y,double *z,double *b){
    double q=x*x+y*y,g=(1+s[2])*s[1]*s[1],disc=1-g*q;
    if (!(disc>0) || !isfinite(q))return 0;
    double root=sqrt(disc),poly=s[7],der=6*s[7];
    for(int j=6;j>=3;j--){poly=poly*q+s[j];der=der*q+(j-1)*s[j];}
    *z=s[1]*q/(1+root)+q*q*poly;
    *b=s[1]/root+2*q*der;
    return isfinite(*z)&&isfinite(*b);
}

static int bounds(const double *s,double q0,double q1,double *sl,double *sh,double *bl,double *bh){
    double g=(1+s[2])*s[1]*s[1],d0=1-g*q0,d1=1-g*q1;
    if(!(d0>0&&d1>0))return 0;
    double r0=sqrt(d0),r1=sqrt(d1);
    double a=s[1]*q0/(1+r0),b=s[1]*q1/(1+r1);
    *sl=mn(a,b);*sh=mx(a,b);
    a=s[1]/r0;b=s[1]/r1;*bl=mn(a,b);*bh=mx(a,b);
    double p0=q0,p1=q1;
    for(int j=3;j<=7;j++){
        int order=j-1;
        a=2*order*s[j]*p0;b=2*order*s[j]*p1;
        *bl+=mn(a,b);*bh+=mx(a,b);
        p0*=q0;p1*=q1;a=s[j]*p0;b=s[j]*p1;
        *sl+=mn(a,b);*sh+=mx(a,b);
    }
    *sl=nextafter(*sl,-INFINITY);*sh=nextafter(*sh,INFINITY);
    *bl=nextafter(*bl,-INFINITY);*bh=nextafter(*bh,INFINITY);
    return isfinite(*sl)&&isfinite(*sh)&&isfinite(*bl)&&isfinite(*bh);
}

static int f_df(const double *s,const double *o,const double *d,double t,double *f,double *df){
    double x=o[0]+t*d[0],y=o[1]+t*d[1],sg,b;
    if(!sag_grad(s,x,y,&sg,&b))return 0;
    *f=o[2]+t*d[2]-s[0]-sg;
    *df=d[2]-b*(x*d[0]+y*d[1]);
    return isfinite(*f)&&isfinite(*df);
}

static int interval_bounds(const double*s,const double*o,const double*d,double lo,double hi,double *fl,double *fh,double *dl,double *dh){
    double xa=o[0]+lo*d[0],ya=o[1]+lo*d[1],xb=o[0]+hi*d[0],yb=o[1]+hi*d[1];
    double qa=xa*xa+ya*ya,qb=xb*xb+yb*yb,q0=mn(qa,qb),q1=mx(qa,qb);
    double dd=d[0]*d[0]+d[1]*d[1];
    if(dd>0){double tv=-(o[0]*d[0]+o[1]*d[1])/dd;if(tv>lo&&tv<hi){double x=o[0]+tv*d[0],y=o[1]+tv*d[1];q0=x*x+y*y;}}
    double sl,sh,bl,bh;
    if(!bounds(s,mx(0,nextafter(q0,-INFINITY)),nextafter(q1,INFINITY),&sl,&sh,&bl,&bh))return 0;
    double za=o[2]+lo*d[2]-s[0],zb=o[2]+hi*d[2]-s[0];
    *fl=nextafter(mn(za,zb)-sh,-INFINITY);*fh=nextafter(mx(za,zb)-sl,INFINITY);
    double ha=xa*d[0]+ya*d[1],hb=xb*d[0]+yb*d[1];
    double p0=bl*ha,p1=bl*hb,p2=bh*ha,p3=bh*hb;
    *dl=nextafter(d[2]-mx(mx(p0,p1),mx(p2,p3)),-INFINITY);
    *dh=nextafter(d[2]-mn(mn(p0,p1),mn(p2,p3)),INFINITY);
    return 1;
}

/* Return 0 hit,1 escaped physical clear aperture,3 unresolved/invalid. */
static int hit(const double *s,const double *o,const double *d,double *t_hit){
    double rad=s[8]+TOL,g=(1+s[2])*s[1]*s[1];
    if(!(rad>0)||g*s[8]*s[8]>=1)return 3;
    if(g>0)rad=mn(rad,nextafter(sqrt(1/g),0));
    int polynomial=0;for(int j=3;j<8;j++)polynomial|=(s[j]!=0);
    if(s[1]==0&&!polynomial){
        if(fabs(d[2])<1e-14)return 3;
        double t=(s[0]-o[2])/d[2];
        if(t < -TOL)return 1;
        t=mx(t,0);double x=o[0]+t*d[0],y=o[1]+t*d[1];
        if(hypot(x,y)>rad)return 1;
        *t_hit=t;return 0;
    }
    double lo=0,hi=INFINITY,aa=d[0]*d[0]+d[1]*d[1],bb=2*(o[0]*d[0]+o[1]*d[1]),qq=o[0]*o[0]+o[1]*o[1],cc=qq-rad*rad;
    if(aa==0){if(cc>TOL*mx(1,rad))return 1;}
    else{
        double disc=bb*bb-4*aa*cc,tol=32*DBL_EPSILON*(bb*bb+4*aa*(qq+rad*rad));
        if(disc < -tol)return 1;
        double root=sqrt(mx(0,disc)),q=-.5*(bb+copysign(root,bb));
        double t0=q/aa,t1=q!=0?cc/q:0;
        lo=mx(lo,nextafter(mn(t0,t1),-INFINITY));hi=mn(hi,nextafter(mx(t0,t1),INFINITY));
    }
    double sl,sh,bl,bh;
    if(!bounds(s,0,rad*rad,&sl,&sh,&bl,&bh))return 3;
    if(fabs(d[2])>1e-14){
        double a=(s[0]+sl-o[2])/d[2],b=(s[0]+sh-o[2])/d[2];
        lo=mx(lo,mn(a,b)-TOL);hi=mn(hi,mx(a,b)+TOL);
    }else if(o[2]<s[0]+sl-TOL||o[2]>s[0]+sh+TOL)return 1;
    if(hi<lo-TOL||hi < -TOL)return 1;
    if(!isfinite(hi))return 3;
    lo=mx(0,lo);hi=mx(0,hi);
    struct interval{double a,b;int depth;} stack[MAX_DEPTH+3];
    int top=1,visits=0;stack[0]=(struct interval){lo,hi,0};
    while(top){
        struct interval cur=stack[--top];double a=cur.a,b=cur.b,fl,fh,dl,dh;
        if(++visits>MAX_INTERVALS)return 3;
        if(!interval_bounds(s,o,d,a,b,&fl,&fh,&dl,&dh))return 3;
        if(fl>TOL||fh < -TOL)continue;
        double fa,fb,dfa,dfb;
        if(!f_df(s,o,d,a,&fa,&dfa)||!f_df(s,o,d,b,&fb,&dfb))return 3;
        if(dl>0||dh<0){
            if(fa==0){*t_hit=a;return fabs(dfa)>1e-13?0:3;}
            if(fa*fb<0){
                double t=(s[0]-o[2])/d[2];if(!(t>a&&t<b))t=(a+b)/2;
                for(int iteration=0;iteration<72;iteration++){
                    double f,df;if(!f_df(s,o,d,t,&f,&df))return 3;
                    if(fabs(f)<TOL*.1 && fabs(df)>1e-13){*t_hit=t;return 0;}
                    if((fa<0)==(f<0)){a=t;fa=f;}else{b=t;fb=f;}
                    double nt=t-f/df;
                    if(!(nt>a&&nt<b)||!isfinite(nt)||iteration%8==7)nt=(a+b)/2;
                    if(nt==t || b-a<TOL*.01){
                        if(fabs(f)<=8*TOL&&fabs(df)>1e-13){*t_hit=t;return 0;}return 3;
                    }
                    t=nt;
                }
                return 3;
            }
            if(fb==0){*t_hit=b;return fabs(dfb)>1e-13?0:3;}
            continue;
        }
        if(b-a<=TOL||cur.depth>=MAX_DEPTH)return 3;
        double mid=(a+b)/2;if(mid==a||mid==b)return 3;
        stack[top++]=(struct interval){mid,b,cur.depth+1};stack[top++]=(struct interval){a,mid,cur.depth+1};
    }
    return 1;
}

void trace_rays(int N,int S,const double*surfaces,const double*origins,const double*directions,const double*indices,
                double*positions,double*outdirs,double*opls,int*statuses,int*failed_surface,double*history,double*dirhistory,double*oplhistory){
    for(int i=0;i<N;i++){
        double o[3],d[3],opl=0;for(int k=0;k<3;k++){o[k]=origins[3*i+k];d[k]=directions[3*i+k];}
        double dn=norm3(d);int status=0,failed=-1;
        if(!(dn>0)||!isfinite(dn)){status=3;failed=0;}else for(int k=0;k<3;k++)d[k]/=dn;
        for(int j=0;j<S&&!status;j++){
            const double *s=surfaces+9*j;double t,n1=indices[i*(S+1)+j],n2=indices[i*(S+1)+j+1];
            if(!(n1>0&&n2>0)||!isfinite(n1+n2)){status=3;failed=j;break;}
            status=hit(s,o,d,&t);if(status){failed=j;break;}
            double p[3];for(int k=0;k<3;k++)p[k]=o[k]+t*d[k];
            double sg,b;if(!sag_grad(s,p[0],p[1],&sg,&b)){status=3;failed=j;break;}
            double normal[3]={-b*p[0],-b*p[1],1},nn=norm3(normal);
            for(int k=0;k<3;k++)normal[k]/=nn;
            double dot=d[0]*normal[0]+d[1]*normal[1]+d[2]*normal[2];
            if(dot>0){dot=-dot;for(int k=0;k<3;k++)normal[k]=-normal[k];}
            double ci=-dot,eta=n1/n2,disc=1-eta*eta*(1-ci*ci),roundoff=32*DBL_EPSILON*mx(1,eta*eta);
            if(disc < -roundoff){status=2;failed=j;break;}
            if(disc<=roundoff){status=3;failed=j;break;}
            double factor=eta*ci-sqrt(disc),nd[3];for(int k=0;k<3;k++)nd[k]=n1==n2?d[k]:eta*d[k]+factor*normal[k];
            double norm=norm3(nd);if(!(norm>0)||!isfinite(norm)){status=3;failed=j;break;}
            if(nd[2]/norm<=1e-13){status=4;failed=j;break;}
            opl+=n1*t;for(int k=0;k<3;k++){
                o[k]=p[k];d[k]=nd[k]/norm;
                history[(i*S+j)*3+k]=o[k];dirhistory[(i*S+j)*3+k]=d[k];
            }
            oplhistory[i*S+j]=opl;
        }
        for(int k=0;k<3;k++){positions[3*i+k]=o[k];outdirs[3*i+k]=d[k];}
        opls[i]=opl;statuses[i]=status;failed_surface[i]=failed;
    }
}
