#include "diag.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define IDX2(i,j,n1) ((j)*(n1)+(i))
#define IDX3(i,j,k,n1,n2) ((k)*(n1)*(n2) + (j)*(n1) + (i))

void new_diag_solver(diag_solver_type *s, double kt, const double nel[2], solve_fn solve){
    s->kt = kt;
    s->nel[0] = nel[0];
    s->nel[1] = nel[1];
    s->solve = solve;
}

void delete_diag_solver(diag_solver_type *s){
    (void)s;
}

static void get_aufbau_filling(double nel, int n, int *homo, double *occ){
    for(int i=0;i<n;i++) occ[i]=0.0;
    *homo = (int)floor(nel);
    for(int i=0;i<*homo && i<n;i++) occ[i]=1.0;
    if(*homo < n) occ[*homo] = fmod(nel,1.0);
    if(fmod(nel,1.0)>0.5 && *homo < n) (*homo)++;
}

static void get_fermi_filling_internal(int homo, double kt, const double *emo, int n,
                                       double *occ, double *efermi){
    const int max_cycle = 200;
    const double thr = sqrt(DBL_EPSILON);
    const double sqrttiny = sqrt(DBL_MIN);
    double e = 0.5*(emo[homo>0? homo-1:0] + emo[(homo<n? homo: n-1)]);
    double occt = homo;
    for(int cycle=1; cycle<=max_cycle; ++cycle){
        double total=0.0, total_d=0.0;
        for(int iao=0; iao<n; ++iao){
            double fermif=0.0, dfermi=0.0;
            double x = (emo[iao]-e)/kt;
            if(x < 50.0){
                double ex = exp(x);
                fermif = 1.0/(ex+1.0);
                dfermi = ex/ (kt*(ex+1.0)*(ex+1.0));
            }
            occ[iao] = fermif;
            total += fermif;
            total_d += dfermi;
        }
        double change = 0.0;
        if(total_d > sqrttiny) change = (occt-total)/total_d;
        e += change;
        if(fabs(occt-total) <= thr) break;
    }
    *efermi = e;
}

void get_fermi_filling(double nel, double kt, const double *emo, int n,
                       int *homo, double *focc, double *efermi){
    double etmp = 0.0;
    get_aufbau_filling(nel, n, homo, focc);
    if(*homo > 0){
        get_fermi_filling_internal(*homo, kt, emo, n, focc, &etmp);
        *efermi = 0.5*etmp;
    }else{
        *efermi = 0.0;
    }
}

void get_density_matrix(const double *focc, const double *coeff,
                        double *pmat, int n){
    double *scratch = malloc(sizeof(double)*n*n);
    for(int iao=0; iao<n; ++iao){
        for(int jao=0; jao<n; ++jao){
            scratch[IDX2(jao, iao, n)] = coeff[IDX2(jao, iao, n)] * focc[iao];
        }
    }
    for(int iao=0; iao<n; ++iao){
        for(int jao=0; jao<n; ++jao){
            double sum=0.0;
            for(int k=0;k<n;k++){
                sum += scratch[IDX2(jao, k, n)] * coeff[IDX2(iao, k, n)];
            }
            pmat[IDX2(jao, iao, n)] = sum;
        }
    }
    free(scratch);
}

void get_density(diag_solver_type *self, const double *smat, double *hmat,
                 double *eval, double *focc, double *density,
                 int norb, int nspin){
    (void)smat;
    int n = nspin;
    if(n == 2){
        for(int i=0;i<norb*norb*n;i++) hmat[i] *= 2.0;
        for(int spin=0; spin<n; ++spin){
            self->solve(self, &hmat[spin*norb*norb], smat, &eval[spin*norb], norb);
            int homo;
            double efermi;
            get_fermi_filling(self->nel[spin], self->kt, &eval[spin*norb], norb,
                              &homo, &focc[spin*norb], &efermi);
            get_density_matrix(&focc[spin*norb], &hmat[spin*norb*norb],
                               &density[spin*norb*norb], norb);
        }
    } else {
        self->solve(self, &hmat[0], smat, &eval[0], norb);
        for(int i=0;i<norb*2;i++) focc[i]=0.0;
        for(int spin=0; spin<2; ++spin){
            int homo; double efermi;
            get_fermi_filling(self->nel[spin], self->kt, &eval[0], norb,
                              &homo, &focc[spin*norb], &efermi);
        }
        for(int i=0;i<norb;i++) focc[i] += focc[norb+i];
        get_density_matrix(&focc[0], &hmat[0], &density[0], norb);
        for(int i=0;i<norb;i++) focc[i] -= focc[norb+i];
    }
}

void get_wdensity(diag_solver_type *self, const double *smat, double *hmat,
                  double *eval, double *focc, double *density,
                  int norb, int nspin){
    (void)self; (void)smat;
    if(nspin==1){
        // assume focc has two columns
        for(int i=0;i<norb;i++) focc[i] += focc[norb+i];
        double *tmp = malloc(sizeof(double)*norb);
        for(int i=0;i<norb;i++) tmp[i] = focc[i]*eval[i];
        get_density_matrix(tmp, hmat, density, norb);
        free(tmp);
        for(int i=0;i<norb;i++) focc[i] -= focc[norb+i];
    }else{
        double *tmp = malloc(sizeof(double)*norb);
        for(int spin=0; spin<nspin; ++spin){
            for(int i=0;i<norb;i++) tmp[i] = focc[IDX2(i,spin,norb)] * eval[IDX2(i,spin,norb)];
            get_density_matrix(tmp, &hmat[spin*norb*norb], &density[spin*norb*norb], norb);
        }
        free(tmp);
    }
}
