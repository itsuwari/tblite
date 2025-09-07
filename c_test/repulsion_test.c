#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "repulsion/effective.h"

static double drand() { return rand() / (double)RAND_MAX; }

int main(void) {
    srand(1);
    for (int t = 0; t < 10; ++t) {
        int nat = 1 + rand() % 4;
        int nid = 1 + rand() % 3;
        int *id = malloc(sizeof(int) * nat);
        double *xyz = malloc(sizeof(double) * 3 * nat);
        for (int i = 0; i < nat; ++i) {
            id[i] = rand() % nid;
            xyz[3*i+0] = drand();
            xyz[3*i+1] = drand();
            xyz[3*i+2] = drand();
        }
        double *alpha = malloc(sizeof(double) * nid * nid);
        double *zeff = malloc(sizeof(double) * nid * nid);
        double *kexp = malloc(sizeof(double) * nid * nid);
        double *rexp = malloc(sizeof(double) * nid * nid);
        for (int i = 0; i < nid * nid; ++i) {
            alpha[i] = 0.5 + drand();
            zeff[i] = drand();
            kexp[i] = 1.0 + drand();
            rexp[i] = 1.0 + drand();
        }
        double cutoff = 10.0;
        double *ener_ref = calloc(nat, sizeof(double));
        for (int i = 0; i < nat; ++i) {
            for (int j = 0; j <= i; ++j) {
                int izp = id[i];
                int jzp = id[j];
                double dx = xyz[3*i+0] - xyz[3*j+0];
                double dy = xyz[3*i+1] - xyz[3*j+1];
                double dz = xyz[3*i+2] - xyz[3*j+2];
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 > cutoff*cutoff || r2 < 1e-12) continue;
                double r = sqrt(r2);
                int p = jzp * nid + izp;
                double r1k = pow(r, kexp[p]);
                double exa = exp(-alpha[p] * r1k);
                double r1r = pow(r, rexp[p]);
                double dE = zeff[p] * exa / r1r;
                ener_ref[i] += 0.5 * dE;
                if (i != j) ener_ref[j] += 0.5 * dE;
            }
        }
        double *ener_test = malloc(sizeof(double) * nat);
        effective_repulsion_energy(nat, id, xyz, nid, alpha, zeff, kexp,
                                   rexp, cutoff, ener_test);
        for (int i = 0; i < nat; ++i) {
            if (fabs(ener_ref[i] - ener_test[i]) > 1e-12) {
                fprintf(stderr, "energy mismatch at atom %d\n", i);
                return 1;
            }
        }
        free(id); free(xyz); free(alpha); free(zeff); free(kexp); free(rexp);
        free(ener_ref); free(ener_test);
    }
    printf("OK\n");
    return 0;
}
