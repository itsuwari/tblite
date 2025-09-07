#ifndef TBLITE_BASIS_TYPE_H
#define TBLITE_BASIS_TYPE_H

#include <stddef.h>

typedef struct {
    int ang;        /* angular momentum */
    int nprim;      /* number of primitives */
    double alpha[12];
    double coeff[12];
} cgto_type;

typedef struct {
    int maxl;
    int nsh;
    int nao;
    double intcut;
    double min_alpha;
    int *nsh_id;   /* size: nspecies */
    int *nsh_at;   /* size: nat */
    int *nao_sh;   /* size: nsh */
    int *iao_sh;   /* size: nsh */
    int *ish_at;   /* size: nat */
    int *ao2at;    /* size: nao */
    int *ao2sh;    /* size: nao */
    int *sh2at;    /* size: nsh */
} basis_type;

void new_basis(basis_type *bas, int nat, int nspecies,
               const int *id, const int *nshell,
               cgto_type **cgto, double acc);

double basis_get_cutoff(const basis_type *bas, double acc);
void free_basis(basis_type *bas);

#endif /* TBLITE_BASIS_TYPE_H */
