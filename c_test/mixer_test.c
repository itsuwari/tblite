#include "../c_src/scf/mixer.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void ref_lineq(int m, double *a, double *c){
    for(int k=0;k<m;k++){
        double pivot=a[k*m+k];
        for(int j=k+1;j<m;j++) a[k*m+j]/=pivot;
        c[k]/=pivot;
        for(int i=k+1;i<m;i++){
            double factor=a[i*m+k];
            for(int j=k+1;j<m;j++) a[i*m+j]-=factor*a[k*m+j];
            c[i]-=factor*c[k];
        }
    }
    for(int i=m-1;i>=0;i--){
        for(int j=i+1;j<m;j++) c[i]-=a[i*m+j]*c[j];
    }
}

static void ref_broyden(int n,int mem,double alpha,double *q,double *qlast,
                        double *dq,double *dqlast,double *df,double *u,
                        double *a,double *omega,int iter){
    int itn=iter-1;
    int it1=(itn-1)%mem;
    if(iter==1){
        memcpy(dqlast,dq,sizeof(double)*n);
        memcpy(qlast,q,sizeof(double)*n);
        for(int i=0;i<n;i++) q[i]+=alpha*dq[i];
        return;
    }
    int limit = itn < mem ? itn : mem;
    double wfac=0.01, minw=1.0, maxw=100000.0;
    double w=0.0;
    for(int i=0;i<n;i++) w+=dq[i]*dq[i];
    w=sqrt(w);
    if(w > (wfac/maxw)) w = wfac/w; else w = maxw;
    if(w < minw) w = minw;
    omega[it1]=w;
    int it1p=it1;
    for(int i=0;i<n;i++) df[it1p*n+i]=dq[i]-dqlast[i];
    double inv=0.0;
    for(int i=0;i<n;i++) inv+=df[it1p*n+i]*df[it1p*n+i];
    inv = 1.0/sqrt(inv>1e-300?inv:1e-300);
    for(int i=0;i<n;i++) df[it1p*n+i]*=inv;
    double *beta=malloc(sizeof(double)*limit*limit);
    double *cvec=malloc(sizeof(double)*limit);
    for(int j=itn-limit+1;j<=itn;j++){
        int jp=(j-1)%mem;
        double dot=0.0;
        for(int i=0;i<n;i++) dot+=df[jp*n+i]*df[it1p*n+i];
        a[jp*mem+it1p]=dot;
        a[it1p*mem+jp]=dot;
        double dot2=0.0;
        for(int i=0;i<n;i++) dot2+=df[jp*n+i]*dq[i];
        int cidx = j - (itn - limit + 1);
        cvec[cidx] = omega[jp]*dot2;
    }
    for(int j=itn-limit+1;j<=itn;j++){
        int jp=(j-1)%mem;
        int col = j - (itn - limit + 1);
        for(int i=0;i<limit;i++){
            int ip=(itn - limit + i)%mem;
            beta[i*limit+col]=omega[ip]*omega[jp]*a[ip*mem+jp];
        }
        beta[col*limit+col]+=0.0001;
    }
    ref_lineq(limit,beta,cvec);
    for(int i=0;i<n;i++)
        u[it1p*n+i]=alpha*df[it1p*n+i]+inv*(q[i]-qlast[i]);
    memcpy(dqlast,dq,sizeof(double)*n);
    memcpy(qlast,q,sizeof(double)*n);
    for(int i=0;i<n;i++) q[i]+=alpha*dq[i];
    for(int j=itn-limit+1;j<=itn;j++){
        int jp=(j-1)%mem;
        int col = j - (itn - limit + 1);
        double coeff=omega[jp]*cvec[col];
        for(int i=0;i<n;i++) q[i]-=coeff*u[jp*n+i];
    }
    free(beta); free(cvec);
}

int main(void){
    srand(0);
    int n=3, mem=5; double damp=0.5;
    mixer_type *mix;
    new_mixer(&mix, mem, n, damp);
    double q0[3], q1[3], q2[3];
    for(int i=0;i<n;i++){ q0[i]=(double)rand()/RAND_MAX; q1[i]=(double)rand()/RAND_MAX; q2[i]=(double)rand()/RAND_MAX; }
    mix->set(mix, q0, n);
    mix->diff(mix, q1, n);
    mix->next(mix);
    double out1[3]; mix->get(mix,out1,n);
    mix->set(mix, out1, n);
    mix->diff(mix, q2, n);
    mix->next(mix);
    double out2[3]; mix->get(mix,out2,n);
    double err = mix->get_error(mix);

    double qref[3]; memcpy(qref,q0,sizeof(q0));
    double qlast[3]={0};
    double dqlast[3]={0};
    double *df=calloc(n*mem,sizeof(double));
    double *u=calloc(n*mem,sizeof(double));
    double *a=calloc(mem*mem,sizeof(double));
    double *omega=calloc(mem,sizeof(double));
    double dq1[3]; for(int i=0;i<n;i++) dq1[i]=q1[i]-q0[i];
    ref_broyden(n,mem,damp,qref,qlast,dq1,dqlast,df,u,a,omega,1);
    for(int i=0;i<n;i++) assert(fabs(out1[i]-qref[i])<1e-12);
    double dq2[3]; for(int i=0;i<n;i++) dq2[i]=q2[i]-qref[i];
    ref_broyden(n,mem,damp,qref,qlast,dq2,dqlast,df,u,a,omega,2);
    for(int i=0;i<n;i++) assert(fabs(out2[i]-qref[i])<1e-12);
    double ref_err=0.0; for(int i=0;i<n;i++){ double d=q2[i]-out1[i]; ref_err+=d*d/n; }
    ref_err=sqrt(ref_err);
    assert(fabs(err-ref_err)<1e-12);
    free(df); free(u); free(a); free(omega);
    mix->del(mix); free(mix);
    printf("mixer tests passed\n");
    return 0;
}
