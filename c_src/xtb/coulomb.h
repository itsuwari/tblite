#ifndef TBLITE_XTB_COULOMB_H
#define TBLITE_XTB_COULOMB_H

#include "../coulomb/charge.h"
#include "../coulomb/multipole.h"

/* Container for isotropic and multipole Coulomb terms */
typedef struct {
    int nat;
    double *gamma;         /* size nat*nat */
    double *amat_sd;       /* size 3*nat*nat */
    double *amat_dd;       /* size 3*nat*3*nat */
    double *amat_sq;       /* size 6*nat*nat */
} tb_coulomb;

void tb_coulomb_new(tb_coulomb *c, int nat);
void tb_coulomb_free(tb_coulomb *c);
void tb_coulomb_energy(const tb_coulomb *c, const double *qat,
                       const double *dpat, const double *qpat,
                       double *energy);
void tb_coulomb_potential(const tb_coulomb *c, const double *qat,
                          const double *dpat, const double *qpat,
                          double *vat, double *vdp, double *vqp);

#endif /* TBLITE_XTB_COULOMB_H */
