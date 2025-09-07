#include "mixer_broyden.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void set_1d(mixer_type *base, const double *qvec, int n){
    broyden_mixer *self = (broyden_mixer*)base;
    memcpy(self->q_in + self->iset, qvec, sizeof(double)*n);
    self->iset += n;
}

static void diff_1d(mixer_type *base, const double *qvec, int n){
    broyden_mixer *self = (broyden_mixer*)base;
    for(int i=0;i<n;i++){
        int idx = self->idif + i;
        self->dq[idx] = qvec[i] - self->q_in[idx];
    }
    self->idif += n;
}

static void get_1d(mixer_type *base, double *qvec, int n){
    broyden_mixer *self = (broyden_mixer*)base;
    memcpy(qvec, self->q_in + self->iget, sizeof(double)*n);
    self->iget += n;
}

static double get_error_fn(const mixer_type *base){
    const broyden_mixer *self = (const broyden_mixer*)base;
    double err=0.0;
    for(int i=0;i<self->ndim;i++) err += self->dq[i]*self->dq[i]/self->ndim;
    return sqrt(err);
}

static void lineq(int m, double *a, double *c){
    for(int k=0;k<m;k++){
        double pivot = a[k*m+k];
        for(int j=k+1;j<m;j++) a[k*m+j] /= pivot;
        c[k] /= pivot;
        for(int i=k+1;i<m;i++){
            double factor = a[i*m+k];
            for(int j=k+1;j<m;j++) a[i*m+j] -= factor*a[k*m+j];
            c[i] -= factor*c[k];
        }
    }
    for(int i=m-1;i>=0;i--){
        for(int j=i+1;j<m;j++) c[i] -= a[i*m+j]*c[j];
    }
}

static void broyden_step(broyden_mixer *self, int info[1]){
    int n=self->ndim, mem=self->memory;
    int iter=self->iter;
    int itn=iter-1;
    int it1=(itn-1)%mem;
    info[0]=0;
    if(iter==1){
        memcpy(self->dqlast, self->dq, sizeof(double)*n);
        memcpy(self->qlast_in, self->q_in, sizeof(double)*n);
        for(int i=0;i<n;i++) self->q_in[i] += self->damp*self->dq[i];
        return;
    }
    int limit = itn < mem ? itn : mem;
    double *beta = malloc(sizeof(double)*limit*limit);
    double *c = malloc(sizeof(double)*limit);
    int it1p = it1;
    for(int i=0;i<n;i++){
        self->df[it1p*n+i] = self->dq[i] - self->dqlast[i];
    }
    double inv_norm=0.0;
    for(int i=0;i<n;i++) inv_norm += self->df[it1p*n+i]*self->df[it1p*n+i];
    inv_norm = 1.0/sqrt(inv_norm>1e-300?inv_norm:1e-300);
    for(int i=0;i<n;i++) self->df[it1p*n+i]*=inv_norm;
    for(int j=itn-limit+1;j<=itn;j++){
        int jp=(j-1)%mem;
        double dot=0.0;
        for(int i=0;i<n;i++) dot += self->df[jp*n+i]*self->df[it1p*n+i];
        self->a[jp*mem+it1p]=dot;
        self->a[it1p*mem+jp]=dot;
        double dot2=0.0;
        for(int i=0;i<n;i++) dot2+= self->df[jp*n+i]*self->dq[i];
        c[jp-(itn-limit+1)] = self->omega[jp]*dot2;
    }
    for(int j=itn-limit+1;j<=itn;j++){
        int jp=(j-1)%mem;
        for(int i=0;i<limit;i++){
            int ip=(itn-limit+1+i-1)%mem;
            beta[i*limit+(jp-(itn-limit+1))] = self->omega[ip]*self->omega[jp]*self->a[ip*mem+jp];
        }
        int idx=jp-(itn-limit+1);
        beta[idx*limit+idx] += 0.0001; /* omega0^2 */
    }
    lineq(limit, beta, c);
    for(int i=0;i<n;i++){
        self->u[it1p*n+i] = self->damp*self->df[it1p*n+i] + inv_norm*(self->q_in[i]-self->qlast_in[i]);
    }
    memcpy(self->dqlast, self->dq, sizeof(double)*n);
    memcpy(self->qlast_in, self->q_in, sizeof(double)*n);
    for(int i=0;i<n;i++) self->q_in[i] += self->damp*self->dq[i];
    for(int j=itn-limit+1;j<=itn;j++){
        int jp=(j-1)%mem;
        double coeff = self->omega[jp]*c[jp-(itn-limit+1)];
        for(int i=0;i<n;i++) self->q_in[i] -= coeff*self->u[jp*n+i];
    }
    free(beta); free(c);
}

static void next_fn(mixer_type *base){
    broyden_mixer *self = (broyden_mixer*)base;
    self->iset=0; self->idif=0; self->iget=0;
    self->iter++;
    int info[1];
    broyden_step(self, info);
    (void)info;
}

static void delete_fn(mixer_type *base){
    broyden_mixer *self = (broyden_mixer*)base;
    free(self->df); free(self->u); free(self->a); free(self->dq);
    free(self->dqlast); free(self->qlast_in); free(self->omega); free(self->q_in);
}

void new_broyden(broyden_mixer *self, int ndim, broyden_input input){
    self->base.set = set_1d;
    self->base.diff = diff_1d;
    self->base.get = get_1d;
    self->base.next = next_fn;
    self->base.get_error = get_error_fn;
    self->base.del = delete_fn;
    self->ndim=ndim;
    self->memory=input.memory;
    self->iter=0; self->iset=self->idif=self->iget=0;
    self->damp=input.damp;
    self->df = calloc(ndim*input.memory, sizeof(double));
    self->u = calloc(ndim*input.memory, sizeof(double));
    self->a = calloc(input.memory*input.memory, sizeof(double));
    self->dq = calloc(ndim, sizeof(double));
    self->dqlast = calloc(ndim, sizeof(double));
    self->qlast_in = calloc(ndim, sizeof(double));
    self->omega = calloc(input.memory, sizeof(double));
    self->q_in = calloc(ndim, sizeof(double));
}
