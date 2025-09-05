#include "../c_src/scf/diag.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define IDX2(i,j,n1) ((j)*(n1)+(i))
#define IDX3(i,j,k,n1,n2) ((k)*(n1)*(n2) + (j)*(n1) + (i))

double rand_double(void){
    return (double)rand()/RAND_MAX*2.0 - 1.0;
}

void identity_solve(diag_solver_type *self, double *hmat, const double *smat,
                    double *eval, int norb){
    (void)self; (void)hmat; (void)smat; (void)eval; (void)norb;
}

/* Reference Fermi filling for tests */
static void ref_aufbau(double nel, int n, int *homo, double *occ){
    for(int i=0;i<n;i++) occ[i]=0.0;
    *homo = (int)floor(nel);
    for(int i=0;i<*homo && i<n;i++) occ[i]=1.0;
    if(*homo < n) occ[*homo] = fmod(nel,1.0);
    if(fmod(nel,1.0)>0.5 && *homo<n) (*homo)++;
}

static void ref_fermi_fill(int homo, double kt, const double *emo, int n,
                           double *occ, double *efermi){
    const int max_cycle=200; const double thr=sqrt(DBL_EPSILON);
    const double sqrttiny=sqrt(DBL_MIN);
    double e=0.5*(emo[homo>0?homo-1:0]+emo[(homo<n?homo:n-1)]);
    double occt=homo;
    for(int cycle=1;cycle<=max_cycle;++cycle){
        double total=0.0,total_d=0.0;
        for(int i=0;i<n;i++){
            double ferm=0.0,df=0.0; double x=(emo[i]-e)/kt;
            if(x<50.0){ double ex=exp(x); ferm=1.0/(ex+1.0); df=ex/(kt*(ex+1.0)*(ex+1.0)); }
            occ[i]=ferm; total+=ferm; total_d+=df;
        }
        double change=0.0; if(total_d>sqrttiny) change=(occt-total)/total_d;
        e+=change; if(fabs(occt-total)<=thr) break;
    }
    *efermi=e;
}

static void ref_get_fermi(double nel, double kt, const double *emo, int n,
                          double *occ){
    int homo; double ef;
    ref_aufbau(nel,n,&homo,occ);
    if(homo>0){ ref_fermi_fill(homo,kt,emo,n,occ,&ef); }
}

int main(void){
    srand(0);
    int norb=4, nspin=2;
    double kt=0.1;
    double nel[2]={3.0,2.0};
    diag_solver_type solver; new_diag_solver(&solver, kt, nel, identity_solve);

    double *smat=malloc(sizeof(double)*norb*norb);
    double *hmat=malloc(sizeof(double)*norb*norb*nspin);
    double *eval=malloc(sizeof(double)*norb*nspin);
    double *focc=malloc(sizeof(double)*norb*nspin);
    double *density=malloc(sizeof(double)*norb*norb*nspin);

    for(int i=0;i<norb*norb;i++) smat[i]=rand_double();
    for(int spin=0;spin<nspin;spin++){
        for(int i=0;i<norb;i++){
            for(int j=0;j<norb;j++){
                hmat[IDX3(j,i,spin,norb,norb)] = (i==j)?0.5:0.0;
            }
        }
    }
    for(int i=0;i<norb*nspin;i++) eval[i]=rand_double();
    memset(focc,0,sizeof(double)*norb*nspin);
    memset(density,0,sizeof(double)*norb*norb*nspin);

    /* Reference */
    double *ref_focc=malloc(sizeof(double)*norb*nspin);
    double *ref_density=malloc(sizeof(double)*norb*norb*nspin);
    for(int spin=0;spin<nspin;spin++){
        ref_get_fermi(nel[spin],kt,&eval[spin*norb],norb,&ref_focc[spin*norb]);
        for(int i=0;i<norb;i++){
            for(int j=0;j<norb;j++){
                ref_density[IDX3(j,i,spin,norb,norb)] =
                    (i==j)?ref_focc[IDX2(i,spin,norb)]:0.0;
            }
        }
    }

    get_density(&solver, smat, hmat, eval, focc, density, norb, nspin);

    for(int i=0;i<norb*nspin;i++) assert(fabs(focc[i]-ref_focc[i])<1e-12);
    for(int i=0;i<norb*norb*nspin;i++) assert(fabs(density[i]-ref_density[i])<1e-12);

    /* Test weighted density */
    for(int spin=0;spin<nspin;spin++){
        for(int i=0;i<norb;i++){
            for(int j=0;j<norb;j++){
                hmat[IDX3(j,i,spin,norb,norb)] = (i==j)?1.0:0.0;
            }
        }
    }
    for(int i=0;i<norb*nspin;i++){ focc[i]=rand_double(); eval[i]=rand_double(); }
    memset(density,0,sizeof(double)*norb*norb*nspin);
    for(int spin=0;spin<nspin;spin++){
        for(int i=0;i<norb;i++){
            for(int j=0;j<norb;j++){
                ref_density[IDX3(j,i,spin,norb,norb)] =
                    (i==j)?focc[IDX2(i,spin,norb)]*eval[IDX2(i,spin,norb)]:0.0;
            }
        }
    }
    get_wdensity(&solver, smat, hmat, eval, focc, density, norb, nspin);
    for(int i=0;i<norb*norb*nspin;i++) assert(fabs(density[i]-ref_density[i])<1e-12);

    free(smat); free(hmat); free(eval); free(focc); free(density);
    free(ref_focc); free(ref_density);
    delete_diag_solver(&solver);
    printf("diag tests passed\n");
    return 0;
}
