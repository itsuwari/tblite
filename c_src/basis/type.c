#include "type.h"
#include <stdlib.h>
#include <math.h>

static double clip(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

static double integral_cutoff(double acc) {
    const double min_intcut = 5.0;
    const double max_intcut = 25.0;
    const double max_acc = 1.0e-4;
    const double min_acc = 1.0e+3;
    double val = max_intcut - 10.0 * log10(clip(acc, min_acc, max_acc));
    if (val < min_intcut) val = min_intcut;
    if (val > max_intcut) val = max_intcut;
    return val;
}

void new_basis(basis_type *bas, int nat, int nspecies,
               const int *id, const int *nshell,
               cgto_type **cgto, double acc) {
    bas->nsh_id = (int*)malloc(sizeof(int)*nspecies);
    for (int i=0;i<nspecies;i++) bas->nsh_id[i] = nshell[i];
    bas->intcut = integral_cutoff(acc);

    bas->nsh_at = (int*)malloc(sizeof(int)*nat);
    for (int i=0;i<nat;i++) bas->nsh_at[i] = nshell[id[i]-1];

    bas->nsh = 0;
    for (int i=0;i<nat;i++) bas->nsh += bas->nsh_at[i];

    bas->ish_at = (int*)malloc(sizeof(int)*nat);
    bas->sh2at = (int*)malloc(sizeof(int)*bas->nsh);
    int ii=0;
    for (int i=0;i<nat;i++) {
        bas->ish_at[i] = ii;
        for (int ish=0; ish<bas->nsh_at[i]; ++ish) {
            bas->sh2at[ii+ish] = i;
        }
        ii += bas->nsh_at[i];
    }

    bas->nao_sh = (int*)malloc(sizeof(int)*bas->nsh);
    for (int i=0;i<nat;i++) {
        int isp = id[i]-1;
        ii = bas->ish_at[i];
        for (int ish=0; ish<bas->nsh_at[i]; ++ish) {
            bas->nao_sh[ii+ish] = 2*cgto[isp][ish].ang + 1;
        }
    }

    bas->nao = 0;
    for (int ish=0; ish<bas->nsh; ++ish) bas->nao += bas->nao_sh[ish];

    bas->iao_sh = (int*)malloc(sizeof(int)*bas->nsh);
    bas->ao2sh = (int*)malloc(sizeof(int)*bas->nao);
    bas->ao2at = (int*)malloc(sizeof(int)*bas->nao);
    ii = 0;
    for (int ish=0; ish<bas->nsh; ++ish) {
        bas->iao_sh[ish] = ii;
        for (int iao=0; iao<bas->nao_sh[ish]; ++iao) {
            bas->ao2sh[ii+iao] = ish;
            bas->ao2at[ii+iao] = bas->sh2at[ish];
        }
        ii += bas->nao_sh[ish];
    }

    bas->maxl = 0;
    bas->min_alpha = 1e308;
    for (int isp=0; isp<nspecies; ++isp) {
        for (int ish=0; ish<nshell[isp]; ++ish) {
            if (cgto[isp][ish].ang > bas->maxl)
                bas->maxl = cgto[isp][ish].ang;
            for (int ip=0; ip<cgto[isp][ish].nprim; ++ip) {
                if (cgto[isp][ish].alpha[ip] < bas->min_alpha)
                    bas->min_alpha = cgto[isp][ish].alpha[ip];
            }
        }
    }
}

double basis_get_cutoff(const basis_type *bas, double acc) {
    const double max_cutoff = 40.0;
    double intcut = bas->intcut;
    if (acc > 0.0) intcut = integral_cutoff(acc);
    double cutoff = sqrt(2.0*intcut/bas->min_alpha);
    if (cutoff > max_cutoff) cutoff = max_cutoff;
    return cutoff;
}

void free_basis(basis_type *bas) {
    free(bas->nsh_id);
    free(bas->nsh_at);
    free(bas->nao_sh);
    free(bas->iao_sh);
    free(bas->ish_at);
    free(bas->ao2at);
    free(bas->ao2sh);
    free(bas->sh2at);
}

