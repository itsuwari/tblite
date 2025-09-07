#ifndef TBLITE_REPULSION_EFFECTIVE_H
#define TBLITE_REPULSION_EFFECTIVE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compute screened Coulomb repulsion energy per atom. */
void effective_repulsion_energy(int nat, const int *id, const double *xyz,
                                int nid, const double *alpha,
                                const double *zeff, const double *kexp,
                                const double *rexp, double cutoff,
                                double *energies);

#ifdef __cplusplus
}
#endif

#endif /* TBLITE_REPULSION_EFFECTIVE_H */
