#include "potential.h"
#include <string.h>
#include <stdio.h>

static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s);
    if (!p) {
        fprintf(stderr, "Allocation failed\n");
        exit(1);
    }
    return p;
}

void new_potential(potential_type *self, const structure_type *mol,
                   const basis_type *bas, int nspin, int grad) {
    self->grad = grad ? 1 : 0;
    self->nat = mol->nat;
    self->nsh = bas->nsh;
    self->nao = bas->nao;
    self->nspin = nspin;

    self->vat = xcalloc((size_t)self->nat * nspin, sizeof(double));
    self->vsh = xcalloc((size_t)self->nsh * nspin, sizeof(double));
    self->vao = xcalloc((size_t)self->nao * nspin, sizeof(double));
    self->vdp = xcalloc((size_t)3 * self->nat * nspin, sizeof(double));
    self->vqp = xcalloc((size_t)6 * self->nat * nspin, sizeof(double));

    if (self->grad) {
        self->dvatdr = xcalloc((size_t)3 * self->nat * self->nat * nspin, sizeof(double));
        self->dvatdL = xcalloc((size_t)3 * 3 * self->nat * nspin, sizeof(double));
        self->dvshdr = xcalloc((size_t)3 * self->nat * self->nsh * nspin, sizeof(double));
        self->dvshdL = xcalloc((size_t)3 * 3 * self->nsh * nspin, sizeof(double));
    } else {
        self->dvatdr = NULL;
        self->dvatdL = NULL;
        self->dvshdr = NULL;
        self->dvshdL = NULL;
    }
}

void reset_potential(potential_type *self) {
    size_t nv = (size_t)self->nat * self->nspin;
    memset(self->vat, 0, nv * sizeof(double));
    nv = (size_t)self->nsh * self->nspin;
    memset(self->vsh, 0, nv * sizeof(double));
    nv = (size_t)self->nao * self->nspin;
    memset(self->vao, 0, nv * sizeof(double));
    nv = (size_t)3 * self->nat * self->nspin;
    memset(self->vdp, 0, nv * sizeof(double));
    nv = (size_t)6 * self->nat * self->nspin;
    memset(self->vqp, 0, nv * sizeof(double));
    if (self->grad) {
        nv = (size_t)3 * self->nat * self->nat * self->nspin;
        memset(self->dvatdr, 0, nv * sizeof(double));
        nv = (size_t)9 * self->nat * self->nspin;
        memset(self->dvatdL, 0, nv * sizeof(double));
        nv = (size_t)3 * self->nat * self->nsh * self->nspin;
        memset(self->dvshdr, 0, nv * sizeof(double));
        nv = (size_t)9 * self->nsh * self->nspin;
        memset(self->dvshdL, 0, nv * sizeof(double));
    }
}

void free_potential(potential_type *self) {
    free(self->vat);
    free(self->vsh);
    free(self->vao);
    free(self->vdp);
    free(self->vqp);
    free(self->dvatdr);
    free(self->dvatdL);
    free(self->dvshdr);
    free(self->dvshdL);
    memset(self, 0, sizeof(*self));
}

/* Helper macros for indexing */
#define IDX2(i,j,n1) ((j)*(n1) + (i))
#define IDX3(i,j,k,n1,n2) ((k)*(n1)*(n2) + (j)*(n1) + (i))

static void add_vat_to_vsh(const basis_type *bas, const potential_type *pot) {
    int nat = pot->nat;
    int nsh = pot->nsh;
    int nspin = pot->nspin;
    (void)nsh;
    for (int spin = 0; spin < nspin; ++spin) {
        for (int iat = 0; iat < nat; ++iat) {
            int ii = bas->ish_at[iat];
            double val = pot->vat[IDX2(iat, spin, nat)];
            for (int ish = 0; ish < bas->nsh_at[iat]; ++ish) {
                int idx = IDX2(ii + ish, spin, nsh);
                pot->vsh[idx] += val;
            }
        }
    }
}

static void add_vsh_to_vao(const basis_type *bas, const potential_type *pot) {
    int nsh = pot->nsh;
    int nao = pot->nao;
    int nspin = pot->nspin;
    (void)nao;
    for (int spin = 0; spin < nspin; ++spin) {
        for (int ish = 0; ish < nsh; ++ish) {
            int ii = bas->iao_sh[ish];
            double val = pot->vsh[IDX2(ish, spin, nsh)];
            for (int iao = 0; iao < bas->nao_sh[ish]; ++iao) {
                int idx = IDX2(ii + iao, spin, nao);
                pot->vao[idx] += val;
            }
        }
    }
}

static void add_vao_to_h1(const basis_type *bas, const integral_type *ints,
                          const potential_type *pot, double *h1) {
    int nao = bas->nao;
    int nspin = pot->nspin;
    for (int spin = 0; spin < nspin; ++spin) {
        for (int iao = 0; iao < nao; ++iao) {
            for (int jao = 0; jao < nao; ++jao) {
                double s = ints->overlap[IDX2(jao, iao, nao)];
                double v = pot->vao[IDX2(jao, spin, nao)] + pot->vao[IDX2(iao, spin, nao)];
                h1[IDX3(jao, iao, spin, nao, nao)] -= 0.5 * s * v;
            }
        }
    }
}

static void add_vmp_to_h1(const basis_type *bas, const double *mpint,
                          const double *vmp, int nmp, const potential_type *pot,
                          double *h1) {
    int nao = bas->nao;
    int nat = bas->nat;
    int nspin = pot->nspin;
    (void)nat;
    for (int spin = 0; spin < nspin; ++spin) {
        for (int iao = 0; iao < nao; ++iao) {
            for (int jao = 0; jao < nao; ++jao) {
                int iat = bas->ao2at[iao];
                int jat = bas->ao2at[jao];
                double dot1 = 0.0, dot2 = 0.0;
                for (int k = 0; k < nmp; ++k) {
                    double m1 = mpint[IDX3(k, jao, iao, nmp, nao)];
                    double m2 = mpint[IDX3(k, iao, jao, nmp, nao)];
                    double v1 = vmp[IDX3(k, iat, spin, nmp, nat)];
                    double v2 = vmp[IDX3(k, jat, spin, nmp, nat)];
                    dot1 += m1 * v1;
                    dot2 += m2 * v2;
                }
                h1[IDX3(jao, iao, spin, nao, nao)] -= 0.5 * (dot1 + dot2);
            }
        }
    }
}

static void magnet_to_updown(double *h1, int nao, int nspin) {
    (void)h1; (void)nao; (void)nspin; /* no-op placeholder */
}

void add_pot_to_h1(const basis_type *bas, const integral_type *ints,
                   potential_type *pot, double *h1) {
    int nao = bas->nao;
    int nspin = pot->nspin;
    size_t m = (size_t)nao * nao;
    memcpy(h1, ints->hamiltonian, m * sizeof(double));
    for (int spin = 1; spin < nspin; ++spin) {
        memset(&h1[spin * m], 0, m * sizeof(double));
    }

    add_vat_to_vsh(bas, pot);
    add_vsh_to_vao(bas, pot);
    add_vao_to_h1(bas, ints, pot, h1);
    add_vmp_to_h1(bas, ints->dipole, pot->vdp, 3, pot, h1);
    add_vmp_to_h1(bas, ints->quadrupole, pot->vqp, 6, pot, h1);
    magnet_to_updown(h1, nao, nspin);
}
