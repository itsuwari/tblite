#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ncoord/gfn.h"

static double drand() { return rand() / (double)RAND_MAX; }

int main(void) {
    srand(0);
    for (int t = 0; t < 10; ++t) {
        int nid = 1 + rand() % 5;
        double *rcov = malloc(sizeof(double) * nid);
        for (int i = 0; i < nid; ++i) rcov[i] = 0.5 + drand();
        gfn_ncoord_type nc;
        new_gfn_ncoord(&nc, nid, rcov, 25.0);
        for (int k = 0; k < 20; ++k) {
            size_t izp = rand() % nid;
            size_t jzp = rand() % nid;
            double r = 0.1 + 5.0 * drand();
            double rc = rcov[izp] + rcov[jzp];
            double c1 = 1.0 / (1.0 + exp(-10.0 * (rc / r - 1.0)));
            double c2 = 1.0 / (1.0 + exp(-20.0 * ((rc + 2.0) / r - 1.0)));
            double expect = c1 * c2;
            double got = ncoord_count(&nc, izp, jzp, r);
            if (fabs(expect - got) > 1e-12) {
                fprintf(stderr, "count mismatch\n");
                return 1;
            }
            double d1_exp = exp(-10.0 * (rc / r - 1.0));
            double d1 = (-10.0 * rc * d1_exp) / (r * r * pow(d1_exp + 1.0, 2));
            double d2_exp = exp(-20.0 * ((rc + 2.0) / r - 1.0));
            double d2 = (-20.0 * (rc + 2.0) * d2_exp) / (r * r * pow(d2_exp + 1.0, 2));
            double dexpect = d1 * c2 + c1 * d2;
            double dgot = ncoord_dcount(&nc, izp, jzp, r);
            if (fabs(dexpect - dgot) > 1e-12) {
                fprintf(stderr, "dcount mismatch\n");
                return 1;
            }
        }
        free(rcov);
        free(nc.rcov);
    }
    printf("OK\n");
    return 0;
}
