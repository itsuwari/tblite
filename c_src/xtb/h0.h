#ifndef TBLITE_XTB_H0_H
#define TBLITE_XTB_H0_H

#include <stddef.h>

/* Simplified Hamiltonian container */
typedef struct {
    int mshell; /* maximum number of shells per element */
    int nid;    /* number of distinct elements (Z index) */
    double *selfenergy; /* size mshell*nid */
    double *kcn;        /* size mshell*nid */
    double *kq1;        /* size mshell*nid */
    double *kq2;        /* size mshell*nid */
    double *refocc;     /* size mshell*nid */
} tb_hamiltonian;

void new_hamiltonian(tb_hamiltonian *h0, int mshell, int nid);
void free_hamiltonian(tb_hamiltonian *h0);

void get_selfenergy(const tb_hamiltonian *h0, int nat,
                    const int *id, const int *ish_at,
                    const int *nshell, const double *cn,
                    const double *qat, double *selfenergy,
                    double *dsedcn, double *dsedq);

void get_occupation(const tb_hamiltonian *h0, int nat, int nsh_tot,
                    const int *id, const int *ish_at,
                    const int *nshell, double charge,
                    double *nocc, double *n0at, double *n0sh);

void shift_operator(const double vec[3], double s,
                    const double di[3], const double qi[6],
                    double dj[3], double qj[6]);

#endif /* TBLITE_XTB_H0_H */
