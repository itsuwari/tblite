#include "effective.h"
#include <math.h>
#include <stddef.h>

/* Helper to index pair parameters stored in row-major nid x nid matrices */
static inline size_t idx(int i, int j, int nid) { return (size_t)j * nid + i; }

void effective_repulsion_energy(int nat, const int *id, const double *xyz,
                                int nid, const double *alpha,
                                const double *zeff, const double *kexp,
                                const double *rexp, double cutoff,
                                double *energies) {
    double cutoff2 = cutoff * cutoff;
    for (int i = 0; i < nat; ++i) energies[i] = 0.0;
    for (int i = 0; i < nat; ++i) {
        int izp = id[i];
        const double *xi = xyz + 3 * i;
        for (int j = 0; j <= i; ++j) {
            int jzp = id[j];
            const double *xj = xyz + 3 * j;
            double rij0 = xi[0] - xj[0];
            double rij1 = xi[1] - xj[1];
            double rij2 = xi[2] - xj[2];
            double r2 = rij0 * rij0 + rij1 * rij1 + rij2 * rij2;
            if (r2 > cutoff2 || r2 < 1.0e-12) continue;
            double r = sqrt(r2);
            size_t p = idx(izp, jzp, nid);
            double k = kexp[p];
            double r1k = pow(r, k);
            double exa = exp(-alpha[p] * r1k);
            double r1r = pow(r, rexp[p]);
            double dE = zeff[p] * exa / r1r;
            energies[i] += 0.5 * dE;
            if (i != j) energies[j] += 0.5 * dE;
        }
    }
}
