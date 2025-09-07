#ifndef D4_H
#define D4_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct d4_param {
    double s6;
    double s8;
    double a1;
    double a2;
};

void d4_get_dispersion_matrix(int nat, int mref, int zmax,
                              const int *ref, const int *id,
                              const double *xyz, const struct d4_param *param,
                              const double *r4r2, const double *c6,
                              const double *trans, int ntrans,
                              double cutoff, double *dispmat);

void d4_self_energy(int nat, int mref, const double *dispmat,
                    const double *gw, double *eat);

void d4_self_potential(int nat, int mref, const double *dispmat,
                       const double *gw, const double *dgwdq,
                       double *vat);

#ifdef __cplusplus
}
#endif

#endif
