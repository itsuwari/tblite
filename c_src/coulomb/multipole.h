#ifndef TBLITE_COULOMB_MULTIPOLE_H
#define TBLITE_COULOMB_MULTIPOLE_H

void multipole_mrad(int nat, const int *id, const double *cn,
                    const double *rad, const double *vcn,
                    double shift, double kexp, double rmax,
                    double *mrad, double *dmrdcn);

void multipole_matrix_0d(int nat, const double *xyz, const double *mrad,
                         double kdmp3, double kdmp5,
                         double *amat_sd, double *amat_dd, double *amat_sq);

void multipole_energy(int nat, const double *amat_sd, const double *amat_dd,
                      const double *amat_sq, const double *qat,
                      const double *dpat, const double *qpat, double *energy);

void multipole_potential(int nat, const double *amat_sd, const double *amat_dd,
                         const double *amat_sq, const double *qat,
                         const double *dpat, const double *qpat,
                         double *vdp, double *vat, double *vqp);

#endif /* TBLITE_COULOMB_MULTIPOLE_H */
